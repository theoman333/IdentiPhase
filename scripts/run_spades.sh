#!/bin/bash

#begin with SPAdes assemble
#This script takes the paired-end data and uses spades to assemble the genome.

spades.py --isolate -1 SRR6837539_1.fastq.gz -2 SRR6837539_2.fastq.gz -o home\omri\Identiphase --memory64


ls
# Check if the correct number of arguments is passed
if [ $# -ne 3 ]; then
  echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> <output_directory>"
  exit 1
fi

# Assign command-line arguments to variables
R1=$1
R2=$2
output_dir=$3

# Check if the files exist
if [[ ! -f "$R1" ]]; then
  echo "Error: File $R1 does not exist."
  exit 1
fi

if [[ ! -f "$R2" ]]; then
  echo "Error: File $R2 does not exist."
  exit 1
fi

# Run SPAdes with the provided files
spades.py --isolate -1 "$R1" -2 "$R2" -o "$output_dir" --memory 64

echo "SPAdes assembly complete. Output stored in $output_dir"

#use contigs.fasta file to llok for inverstions 