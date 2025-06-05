#This is the algorithm for finding putative regions. It takes the BAM file as input and outputs a
#BED file with the putative regions.

import pysam
from collections import defaultdict

bam = pysam.AlignmentFile("aligned_reads_sorted.bam", "rb")

window_size = 1000
step_size = 500
min_discordant = 5

# Step 1: Count discordant reads in windows
discordant_windows = defaultdict(list)  # key: (chrom, window_start), value: list of positions

for read in bam:
    if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped:
        # Check for flipped orientation
        if read.is_reverse == read.mate_is_reverse:
            pos = min(read.reference_start, read.next_reference_start)
            window_start = (pos // step_size) * step_size
            discordant_windows[(read.reference_name, window_start)].append(pos)

# Step 2: Write results — both window and narrow region
with open("putative_inversions.bed", "w") as out:
    for i, ((chrom, window_start), positions) in enumerate(discordant_windows.items(), 1):
        if len(positions) >= min_discordant:
            window_end = window_start + window_size
            actual_start = min(positions)
            actual_end = max(positions)
            #out.write(f"{chrom}\t{window_start}\t{window_end}\tinv_window_{i}\n")
            out.write(f"{chrom}\t{actual_start}\t{actual_end}\tinv_actual_{i}\n")

bam.close()
print("Done. Wrote window and actual inversion regions to putative_inversions.bed")
