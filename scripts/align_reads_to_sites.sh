#!/bin/bash

#align the original reads to the site canditates.

bowtie -p 4 -a --best --strata inversion_candidates -1 SRR6837539_1.fastq.gz -2 SRR6837539_2.fastq.gz -S | \
samtools view -@ 4 -F 4 -h | \
sam2bed -d | \
sortBed | \
cut -f 1-4,7 > beds/inversion_site_alignment.bed