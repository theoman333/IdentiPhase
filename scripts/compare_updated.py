import pandas as pd
from Bio import SeqIO

# 1. Load your BED with inv_id
bed = pd.read_csv("beds/putative_inversions.bed", sep="\t", header=None,
                  names=["Chrom", "Start", "End", "inv_id"])

# 2. Load ratios.tsv and filter for non-zero R_to_total
ratios = pd.read_csv("ratios.tsv", sep="\t")
ratios_filtered = ratios[ratios["R_to_total"] != 0]

# 3. Filter BED to only keep inv_id present in filtered ratios
bed_filtered = bed[bed["inv_id"].isin(ratios_filtered["inv_id"])]

# 4. Load your genome FASTA file
genome_fasta = "contigs.fasta"  # replace with your actual fasta file path
genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# 5. Extract sequences for filtered inversions
with open("filtered_inversion_sequences.fa", "w") as out_f:
    for _, row in bed_filtered.iterrows():
        chrom = row["Chrom"]
        start = row["Start"]
        end = row["End"]
        inv_id = row["inv_id"]

        # Extract sequence (assuming 0-based start, 1-based end as in BED)
        seq = genome[chrom].seq[start:end]

        # Write to fasta file
        out_f.write(f">{inv_id}|{chrom}:{start}-{end}\n{seq}\n")

print(f"Extracted {len(bed_filtered)} inversion sequences to filtered_inversion_sequences.fa")
