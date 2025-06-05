from collections import defaultdict

input_file = "inversion_site_alignment.bed"  # your BED file
output_file = "inversion_read_ratios.txt"    # output file

# Store counts for each site
counts = defaultdict(lambda: {'F': 0, 'R': 0})

# Parse BED file
with open(input_file) as f:
    for line in f:
        chrom = line.strip().split('\t')[0]  # e.g., inv_actual_1011_F
        parts = chrom.rsplit('_', 1)
        site, orientation = parts[0], parts[1]
        if orientation in ['F', 'R']:
            counts[site][orientation] += 1

# Write output with F / (F + R) ratio
with open(output_file, 'w') as out:
    out.write("Site\tForward_Reads\tReverse_Reads\tF/(F+R)_Ratio\n")
    for site, vals in sorted(counts.items()):
        f_reads = vals['F']
        r_reads = vals['R']
        total_reads = f_reads + r_reads
        ratio = f_reads / total_reads if total_reads > 0 else "NA"
        out.write(f"{site}\t{f_reads}\t{r_reads}\t{ratio:.4f}\n")

print(f"✅ Ratios written to {output_file}")
