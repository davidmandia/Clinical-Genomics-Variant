#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.bam output.sam"
    exit 1
fi

# Assign arguments to variables
input_bam=$1
output_sam=$2

# Convert BAM to SAM
samtools view -h -o $output_sam $input_bam

# Check if the conversion was successful
if [ $? -eq 0 ]; then
    echo "Conversion successful: $output_sam created."
else
    echo "Error during conversion."
fi





