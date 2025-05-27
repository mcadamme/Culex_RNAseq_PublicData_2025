#! /bin/bash

#This is the script I used to run fastqc on the RNAseq data from the SRA
#MLF 05262025

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

