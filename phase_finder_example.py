#!/usr/bin/env python3
"""
phase_finder_example.py ────────────────────────────────────────────────────────
A **teaching‑grade** re‑implementation of the core Phase Finder workflow,
_now with line‑by‑line commentary_.

Goal
────
• Show **exactly** how ref‑based inversion calling works without hidden magic.
• Keep a single readable file you can tweak in a Jupyter cell or embed in a
  Snakemake rule.

Dependencies (install via conda or pip + system packages):
  • Python ≥3.8, Biopython, pysam
  • EMBOSS (for `einverted`), Bowtie 2, SAMtools

Typical usage
-------------
```bash
phase_finder_example.py \
  --reference genome.fna \
  --reads reads_R1.fastq.gz reads_R2.fastq.gz \
  --outdir phase_finder_out --threads 8
```
Output folder structure
```
phase_finder_out/
  ├─ IR_pairs.tsv            # coordinates of every inverted‑repeat pair
  ├─ candidates.fna          # F/R sequences ready for mapping
  ├─ bt2_index.*             # Bowtie 2 index files
  ├─ alignments.bam          # read mappings (sorted + indexed)
  └─ orientation_counts.tsv  # ON/OFF ratios per locus
```

IMPORTANT  This script is **didactic**: error trapping is minimal and no clever
optimisations are applied.  For production work, parallelise the alignment step
and add QC (coverage filters, MAPQ cut‑offs, etc.).
"""

# ──────────────────────────────────────────────────────────────────────────────
# Imports & type aliases
# ──────────────────────────────────────────────────────────────────────────────

# ─ py stdlib ─────────────────────────────────────────────────────────────────
import argparse               # parse CLI flags
import os
import re                     # regex for einverted parsing
import subprocess             # run external tools
from pathlib import Path      # nicer path handling than os.path
from typing import List, Tuple

# ─ 3rd‑party ─────────────────────────────────────────────────────────────────
from Bio import SeqIO         # reading/writing FASTA/FASTQ
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam                  # lightweight SAM/BAM wrapper

# ──────────────────────────────────────────────────────────────────────────────
# Generic helper: run a shell command verbosely and abort on failure
# ──────────────────────────────────────────────────────────────────────────────

def run_cmd(cmd: List[str], msg: str) -> None:
    """Pretty‑print *cmd*, run it, and raise if the exit code ≠ 0."""
    print(f"[+] {msg}")
    print("    $", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)

# ──────────────────────────────────────────────────────────────────────────────
# STEP 1 ─ Locate inverted repeats with EMBOSS *einverted*
# ──────────────────────────────────────────────────────────────────────────────


def find_inverted_repeats(
    ref_fasta: Path,
    outfile: Path,
    min_repeat: int = 8,
    max_repeat: int = 100,
    max_mismatch: int = 4,
    gap: int = 12,
) -> None:
    """Call *einverted*, then convert its verbose output into a neat TSV.

    Parameters
    ----------
    ref_fasta  : Path   Reference genome (one or many contigs) in FASTA.
    outfile    : Path   Destination .tsv path (5 cols: contig / A0 A1 B0 B1).
    min_repeat : int    Shortest IR motif to keep.
    max_repeat : int    Longest IR motif to keep.
    max_mismatch: int   Max mismatches allowed in IR comparison.
    gap        : int    Max permitted loop (gap) length inside the IR arms.
    """

    # *einverted* spits human‑readable text; store it temporarily then parse.
    einv_raw = outfile.with_suffix(".einv.txt")

    # Build the command – see EMBOSS docs for flag meanings.
    run_cmd(
        [
            "einverted",
            "-sequence", str(ref_fasta),      # input genome
            "-gap", str(gap),                # IR loop ≤ gap bp
            "-match", "3", "-mismatch", "-4", # scoring matrix
            "-minrepeat", str(min_repeat),
            "-maxrepeat", str(max_repeat),
            "-threshold", "30",              # min alignment score to keep
            "-outfile", str(einv_raw),
        ],
        "Scanning genome for inverted repeats (IRs)",
    )

    # Regex captures 4 coordinate integers per hit line in *einverted* output.
    hit_line = re.compile(r"^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)")
    records: List[Tuple[str, int, int, int, int]] = []
    contig = None  # will be updated when a line starts with '>'

    with open(einv_raw) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                # Header line: store current contig ID (strip '>').
                contig = line[1:].split()[0]
                continue
            m = hit_line.match(line)
            if m and contig:
                a0, a1, b0, b1 = map(int, m.groups())
                # Convert 1‑based inclusive → 0‑based half‑open coordinates.
                # Phase Finder only needs bounding boxes of the two IR arms.
                records.append((contig, a0 - 1, a1, b0 - 1, b1))

    # Dump to a clean tabular file.
    with open(outfile, "w") as out:
        out.write("contig\tIR_A_start\tIR_A_end\tIR_B_start\tIR_B_end\n")
        for rec in records:
            out.write("\t".join(map(str, rec)) + "\n")

# ──────────────────────────────────────────────────────────────────────────────
# STEP 2 ─ For each IR pair, build forward (F) and flipped (R) sequences
# ──────────────────────────────────────────────────────────────────────────────


def load_reference(ref_fasta: Path):
    """Return a dictionary keyed by contig ID → Bio.SeqRecord."""
    return {rec.id: rec for rec in SeqIO.parse(ref_fasta, "fasta")}


def make_candidates_fasta(
    ref_fasta: Path,
    ir_tsv: Path,
    flank: int,
    out_fasta: Path,
) -> None:
    """Generate the two orientations (F, R) for every candidate locus.

    • *Flank* is extra context so reads map uniquely (default 200 bp).
    • The flipped version is produced by reverse‑complementing only the IR
      *core* (between the repeat arms), leaving flanks unchanged so read
      aligners can distinguish the two orientations.
    """

    ref = load_reference(ref_fasta)
    out_records: List[SeqRecord] = []

    for idx, line in enumerate(open(ir_tsv)):
        if line.startswith("contig"):
            continue  # skip header
        contig, a0, a1, b0, b1 = line.rstrip().split("\t")
        a0, a1, b0, b1 = map(int, (a0, a1, b0, b1))

        # ─ Region coordinates including flanking context ──────────────
        start = max(0, a0 - flank)
        end = min(len(ref[contig]), b1 + flank)

        # Forward orientation = as in reference.
        forward_seg = ref[contig].seq[start:end]

        # Construct the flipped orientation:
        #   flanks stay the same; middle block (a0:b1) is reverse‑complemented.
        flipped_core = ref[contig].seq[a0:b1].reverse_complement()
        flipped_seg = ref[contig].seq[start:a0] + flipped_core + ref[contig].seq[b1:end]

        locus_id = f"locus{idx + 1:04d}_{contig}_{start + 1}_{end}"
        out_records.append(SeqRecord(forward_seg, id=f"{locus_id}|F", description="forward"))
        out_records.append(SeqRecord(flipped_seg, id=f"{locus_id}|R", description="flipped"))

    # Write all loci to a single FASTA – Bowtie 2 will index the lot.
    SeqIO.write(out_records, out_fasta, "fasta")

# ──────────────────────────────────────────────────────────────────────────────
# STEP 3 ─ Build Bowtie 2 index and align reads
# ──────────────────────────────────────────────────────────────────────────────


def build_bowtie2_index(fasta: Path, prefix: Path, threads: int) -> None:
    """Create an FM‑index of the F/R sequences using Bowtie 2‑build."""
    run_cmd([
        "bowtie2-build",
        "--threads", str(threads),
        str(fasta),        # input FASTA
        str(prefix),       # output prefix (generates 6 files)
    ], "Building Bowtie 2 index")


def align_reads_bt2(
    index_prefix: Path,
    r1: Path,
    r2: Path,
    bam_out: Path,
    threads: int,
) -> None:
    """Map paired‑end reads to the F/R index –> sorted, indexed BAM."""

    # 1. Alignment → SAM
    sam_tmp = bam_out.with_suffix(".sam")
    run_cmd(
        [
            "bowtie2",
            "-x", str(index_prefix),  # index prefix produced above
            "-1", str(r1), "-2", str(r2),
            "--very-sensitive",       # preset with relaxed seed settings
            "--threads", str(threads),
            "-S", str(sam_tmp),
        ],
        "Aligning reads to F/R database",
    )

    # 2. Convert SAM → BAM, sort by coordinate (required for indexing).
    run_cmd([
        "samtools", "sort",
        "-@", str(threads),
        "-o", str(bam_out),
        str(sam_tmp),
    ], "Sorting alignments")

    # 3. Build .bai index for random access.
    run_cmd(["samtools", "index", str(bam_out)], "Indexing BAM")

    # Housekeeping: remove large intermediate SAM.
    sam_tmp.unlink(missing_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# STEP 4 ─ Count how many reads back F vs R at each locus
# ──────────────────────────────────────────────────────────────────────────────


def count_orientation(bam: Path, tsv_out: Path) -> None:
    """Loop over every alignment, bucket counts by (locus, orientation).

    Orientation is encoded in the reference name of each BAM record, e.g.
    "locus0003_chr1_12345_13000|R".  We tally *all* reads (proper + orphan
    pairs).  For production you might require:
      • `read.is_proper_pair` AND NOT `read.is_secondary` AND NOT `read.is_qcfail`
      • MAPQ ≥ 20
    """

    bamfile = pysam.AlignmentFile(bam, "rb")
    counts = {}  # {locus: {"F": n, "R": n}}

    for aln in bamfile:
        if aln.is_unmapped:
            continue
        ref_name = bamfile.get_reference_name(aln.reference_id)
        locus, orient = ref_name.rsplit("|", 1)  # fast split from the right
        buckets = counts.setdefault(locus, {"F": 0, "R": 0})
        buckets[orient] += 1

    with open(tsv_out, "w") as out:
        out.write("locus\tF_reads\tR_reads\tR_fraction\n")
        for locus, bucket in sorted(counts.items()):
            f = bucket["F"]
            r = bucket["R"]
            ratio = r / (f + r) if (f + r) else 0.0
            out.write(f"{locus}\t{f}\t{r}\t{ratio:.4f}\n")

# ──────────────────────────────────────────────────────────────────────────────
# CLI glue – parse arguments, run steps in order
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser("Minimal, fully‑commented Phase Finder workflow")
    p.add_argument("--reference", required=True, type=Path, help="Reference genome FASTA")
    p.add_argument("--reads", nargs=2, required=True, type=Path, metavar=("R1", "R2"),
                   help="Paired‑end FASTQ files")
    p.add_argument("--outdir", required=True, type=Path, help="Output directory")
    p.add_argument("--threads", type=int, default=4, help="CPU cores for Bowtie 2 / SAMtools")
    p.add_argument("--flank", type=int, default=200, help="Extra bp to add either side of IR core")
    args = p.parse_args()

    # Create output dir first – abort early if path invalid.
    args.outdir.mkdir(parents=True, exist_ok=True)

    # 1. Locate IRs
    ir_tsv = args.outdir / "IR_pairs.tsv"
    find_inverted_repeats(args.reference, ir_tsv)

    # 2. Build F/R FASTA
    candidates_fna = args.outdir / "candidates.fna"
    make_candidates_fasta(args.reference, ir_tsv, args.flank, candidates_fna)

    # 3. Index & align
    index_prefix = args.outdir / "bt2_index"
    build_bowtie2_index(candidates_fna, index_prefix, args.threads)

    bam_out = args.outdir / "alignments.bam"
    align_reads_bt2(index_prefix, args.reads[0], args.reads[1], bam_out, args.threads)

    # 4. Count F vs R support
    orientation_tsv = args.outdir / "orientation_counts.tsv"
    count_orientation(bam_out, orientation_tsv)

    print(f"[✓] Complete!  Per‑locus ON/OFF ratios: {orientation_tsv}")

# Only run *main()* when the file is executed directly, not imported.
if __name__ == "__main__":
    main()
