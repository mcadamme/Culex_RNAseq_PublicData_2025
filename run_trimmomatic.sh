#! /bin/bash

#This is the script I used to run Trimmomatic on the RNAseq data from the SRA
#MLF 05262025

# Setting the base directory
BASE_DIR="/home/megan/Desktop/Public_mosquito_rnaseq"

for first_level in "$BASE_DIR"/SRR*/; do

    for second_level in "$first_level"*/; do
    
    	#getting SRR name
	SUBDIR_NAME=$(basename "$(dirname "${second_level}")")
        
        # Expand both files
        files=("${second_level}"*.fastq.gz)
        
        # Define output file names
    	OUT_P1="${SUBDIR_NAME}"_trimmed_R1_paired.fastq.gz
    	OUT_UP1="${SUBDIR_NAME}"_trimmed_R1_unpaired.fastq.gz
    	OUT_P2="${SUBDIR_NAME}"_trimmed_R2_paired.fastq.gz
    	OUT_UP2="${SUBDIR_NAME}"_trimmed_R2_unpaired.fastq.gz
        
        # Run Trimmomatic on each file in this directory
        java -jar /home/megan/src/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 "${files[0]}" "${files[1]}" "$OUT_P1" "$OUT_UP1" "$OUT_P2" "$OUT_UP2" ILLUMINACLIP:/home/megan/src/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50
    
    done
done

