import pandas as pd

# Read your bed file (inv_actuals)
inv_df = pd.read_csv('beds/filtered_merged_results.bed', sep='\t', header=None, names=['contig', 'start', 'end', 'inv_id', 'ratio'])

# Read your PhaseFinder results
phase_df = pd.read_csv('PhaseFinderResults.out.txt.txt', sep='\t')

# Extract start1 and end2 from the PhaseFinder 'ID' column
def extract_start1_end2(phase_id):
    parts = phase_id.split(':')[1]
    nums = list(map(int, parts.split('-')))
    start1 = nums[0]
    end2 = nums[3]
    return pd.Series({'phase_start1': start1, 'phase_end2': end2})

phase_df[['phase_start1', 'phase_end2']] = phase_df['ID'].apply(extract_start1_end2)

# Now match: for each inv_actual, find phaseFinder IDs where start1 or end2 is close (±500)
matches = []

for idx, inv_row in inv_df.iterrows():
    inv_start = inv_row['start']
    inv_end = inv_row['end']
    inv_id = inv_row['inv_id']
    ratio = inv_row['ratio']

    # Find matches
    nearby = phase_df[
        ((abs(phase_df['phase_start1'] - inv_start) <= 500) | 
         (abs(phase_df['phase_end2'] - inv_start) <= 500) |
         (abs(phase_df['phase_start1'] - inv_end) <= 500) |
         (abs(phase_df['phase_end2'] - inv_end) <= 500))
    ]

    if not nearby.empty:
        for _, phase_row in nearby.iterrows():
            matches.append({
                'inv_id': inv_id,
                'inv_start': inv_start,
                'inv_end': inv_end,
                'inv_ratio': ratio,
                'matched_phase_id': phase_row['ID'],
                'phase_start1': phase_row['phase_start1'],
                'phase_end2': phase_row['phase_end2'],
                'Pe_ratio': phase_row['Pe_ratio']
            })
    else:
        matches.append({
            'inv_id': inv_id,
            'inv_start': inv_start,
            'inv_end': inv_end,
            'inv_ratio': ratio,
            'matched_phase_id': None,
            'phase_start1': None,
            'phase_end2': None,
            'Pe_ratio': None
        })

# Save or display the results
matches_df = pd.DataFrame(matches)
matches_df.to_csv('matched_phasefinder_results.csv', index=False)
print(matches_df)
