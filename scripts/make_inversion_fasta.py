#This code takes the coordinates of the inversion sites from the bed file
# and puts them in a fasta file with a flank of 500 nt.

from Bio import SeqIO
from Bio.Seq import Seq

# Load genome FASTA as a dictionary
genome = SeqIO.to_dict(SeqIO.parse("contigs.fasta", "fasta"))

# Parameters
flank = 500  # How much upstream/downstream to include

# Open BED file and read coordinates of inversion sites
inversions = []
with open("putative_inversions.bed") as bed:
    for line in bed:
        parts = line.strip().split('\t')
        contig = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        name = parts[3] if len(parts) > 3 else f"inv_{start}"
        inversions.append((contig, start, end, name))

# Write sequences in fasta file
with open("inversion_candidates.fasta", "w") as out_fasta:
    for i, (contig, start, end, name) in enumerate(inversions, 1):
        if contig not in genome:
            print(f"Warning: {contig} not in genome. Skipping.")
            continue

        seq = genome[contig].seq
        flank_L = seq[max(0, start - flank):start]
        inv_seq = seq[start:end]
        flank_R = seq[end:end + flank]

        full_F = flank_L + inv_seq + flank_R
        full_R = flank_L + inv_seq.reverse_complement() + flank_R

        out_fasta.write(f">{name}_F\n{full_F}\n")
        out_fasta.write(f">{name}_R\n{full_R}\n")

print(f"✅ Created inversion_candidates.fasta with {len(inversions)*2} sequences.")
