#!/bin/bash


#update the code so it takes input here and then uses it on run_spades.sh
#begin with SPAdes assemble
#will not run for now because takes a long time
#This is the "main" of the code 

#Step 1, using spades to assemble the makeshift reference from the samples.
# ./run_spades.sh


#now allign reads to contigs using bowtie

# Step 1: Build the genome index using Bowtie2 (only needs to be done once)
echo "Building Bowtie2 index for the genome..."
#need to change path name eventually.
bowtie2-build contigs.fasta genome_index

#   Step 2: Align the segments to the genome using Bowtie2 (for each segment!) outputs sam file
# echo Aligning segments to the genome...
bowtie2 -x genome_index -1 spades_results/C_PD_TP7_1p1_sampled3.5_R1.fastq.gz -2 spades_results/C_PD_TP7_1p1_sampled3.5_R2.fastq.gz -S aligned_reads.sam
#We use spade_results because it is already in contigs format,  more useful than using samples. (not sure about this acutally.)

# use samtools to convert to Bam
# samtools view -bS aligned_reads.sam | samtools sort -o aligned_reads_sorted.bam
#samtools index aligned_reads_sorted.bam


# #here continue by finnding dosconcorant reads '
#look at grok, chatGPT
