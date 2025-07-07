Workflow
run_spades.sh

Step 1) Use spades to build makeshift reference. 
This outputs contigs.fasta


IdentiPhase.sh
Step 2) Use bowtie to create genome index 
Step 3) Align the segments to the genome using Bowtie2 and create SAM file (which details the allignment)
output: aligned_reads.sam
Step 3a) Convert SAM file to BAM file (smaller file, more efficient. )
output: aligned_reads_sorted_bam

Step 4)
putative.py - the algorithm for finding putative regions. 
outputs the regions to putative_inversions.bed

Step 5)
make_inversion_fasta.py
Takes the coordinates of the inversion sites from the bed file  and puts them in a fasta file with a flank of 500 nt/
We do this so we can use part of the PhaseFinder code to find the ratios of F/R/
output: inversion_canditates.fasta

Step 6)
align_reads_to_sites.sh
This script takes inversion_candidates.fasta and uses bowtie to align the original samples to the sites.
This is how we caclulate the ratio of alligned forward and reverese reads, similar to PhaseFinder/
Output: beds/inversion_site_alignment.bed
This bed file tells me the site that the read was aligned to and 
the start and end position of the alignment on the site.
Now I know how many were aligned to forward or reverse.
 
Step 7)
count_inversions_ratio.py
output: inversion_read_ratios.txt
Step 8)
Use PhaseFioinder code to filter results.
->ratios.tsv

where i am now:
need to compare my results with PhaseFinder results. my results are saved in ratios.tsv
use einverted to validate results. 

Need to find a way to organize this. Run e-inverted on suspected finds!
Locate my putatate inversion fastas and run einverted on them.


merge_results.py:
merges the bed file with fasta so we know the genome locations.


log:
Identiphase found area 626-747:
PhaseFinder knew this area, found 88 in the forward but none in the reverse. 

*padded putatate_inversions.bed using bedtools slop, output -> padded_putative_inversions
*created contigs.genome file using samtools faidx on contigs.fasta

then createds regions.fa using bedtools getfasta on padded_putative_inversions.bed
to extract the actual sequences of these regions.

-then ran einverted on this
 ------------------------
 this was not so good. so i combined ratios.tsv with putate_inversions.bed, then got the fasta sequence of the inverted regions
 and saved to filtered_inverstion_sequences.fa.
 Then i ran einverted and saved output to ir.results.txt but this tool a long time, need to run again and see why it's not working well.
 what i need to do is run einverted on the sites i want with only the reads in that area.
 need to find a way to validate results.
 need to run einverted
 -------------------------------------
Workflow
run_spades.sh

Step 1) Use spades to build makeshift reference. 
This outputs contigs.fasta


IdentiPhase.sh
Step 2) Use bowtie to create genome index 
Step 3) Align the segments to the genome using Bowtie2 and create SAM file (which details the allignment)
output: aligned_reads.sam
Step 3a) Convert SAM file to BAM file (smaller file, more efficient. )
output: aligned_reads_sorted_bam

Step 4)
putative.py - the algorithm for finding putative regions. 
outputs the regions to putative_inversions.bed

Step 5)
make_inversion_fasta.py
Takes the coordinates of the inversion sites from the bed file  and puts them in a fasta file with a flank of 500 nt/
We do this so we can use part of the PhaseFinder code to find the ratios of F/R/
output: inversion_canditates.fasta

Step 6)
align_reads_to_sites.sh
This script takes inversion_candidates.fasta and uses bowtie to align the original samples to the sites.
This is how we caclulate the ratio of alligned forward and reverese reads, similar to PhaseFinder/
Output: beds/inversion_site_alignment.bed
This bed file tells me the site that the read was aligned to and 
the start and end position of the alignment on the site.
Now I know how many were aligned to forward or reverse.
 
Step 7)
count_inversions_ratio.py
output: inversion_read_ratios.txt
Step 8)
Use PhaseFioinder code to filter results.
->ratios.tsv

where i am now:
need to compare my results with PhaseFinder results. my results are saved in ratios.tsv
use einverted to validate results. 

Need to find a way to organize this. Run e-inverted on suspected finds!
Locate my putatate inversion fastas and run einverted on them.


merge_results.py:
merges the bed file with fasta so we know the genome locations.


log:
Identiphase found area 626-747:
PhaseFinder knew this area, found 88 in the forward but none in the reverse. 

*padded putatate_inversions.bed using bedtools slop, output -> padded_putative_inversions
*created contigs.genome file using samtools faidx on contigs.fasta

then createds regions.fa using bedtools getfasta on padded_putative_inversions.bed
to extract the actual sequences of these regions.

-then ran einverted on this
 ------------------------
 this was not so good. so i combined ratios.tsv with putate_inversions.bed, then got the fasta sequence of the inverted regions
 and saved to filtered_inverstion_sequences.fa.
 Then i ran einverted and saved output to ir.results.txt but this tool a long time, need to run again and see why it's not working well.
 what i need to do is run einverted on the sites i want with only the reads in that area.
 need to find a way to validate results.
 need to run einverted
 -------------------------------------
 **next step: validation. use the sites found in ratios.tsv. 
 allign reads to site and calculate percetange of inverted reads.**

