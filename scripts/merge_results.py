import pandas as pd

# Load the ratios file and the bed file
ratios_df = pd.read_csv('ratios.tsv', sep='\t', header=None, names=['inversion', 'ratio'])
bed_df = pd.read_csv('beds/putative_inversions.bed', sep='\t', header=None, names=['contig', 'start', 'end', 'inversion'])

# Filter the bed file to keep only rows with inversion IDs that exist in the ratios file
filtered_bed_df = bed_df[bed_df['inversion'].isin(ratios_df['inversion'])]

# Merge the filtered bed file with the ratios file
merged_df = pd.merge(filtered_bed_df, ratios_df, on='inversion', how='left')

# Save the result to a new file
merged_df.to_csv('beds/merged_results.bed', sep='\t', header=False, index=False)

