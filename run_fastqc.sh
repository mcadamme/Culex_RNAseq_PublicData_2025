#! /bin/bash

#This is the script I used to run fastqc on the RNAseq data from the SRA
#MLF 05262025


#!/bin/bash

# Set the base directory
BASE_DIR="/home/megan/Desktop/Public_mosquito_rnaseq"

# Loop through the first and second level directories
for first_level in "$BASE_DIR"/SRR*/; do
    for second_level in "$first_level"*/; do
        # Run fastqc on each file in each directory
            /home/megan/src/FastQC/fastqc *.fastq.gz

    done
done



# Setting the base directory
BASE_DIR="/home/megan/Desktop/Public_mosquito_rnaseq"

for first_level in "$BASE_DIR"/SRR*/; do
    for second_level in "$first_level"*/; do
        
        # Expand both files
        files=("${second_level}"*.fastq.gz)
        
        # Run fastqc on each file in this directory
            /home/megan/src/FastQC/fastqc "${files}"
    
    done
done

