#! /bin/bash

#This is the script I used for sam to bam conversion
#MF 05292025

#Setting the base directory
BASE_DIR="/home/megan/Desktop/Public_mosquito_rnaseq"

for first_level in "$BASE_DIR"/SRR*/; do
     
     for second_level in "$first_level"*/; do
     
     cd "${second_level}"
    
    
     #getting SRR name
     SUBDIR_NAME=$(basename "$(dirname "${second_level}")")
     
     for sample in "${SUBDIR_NAME}".sam; do
     echo "$sample"
     
     #convert file from SAM to BAM format
     samtools view -bS $sample -o ${SUBDIR_NAME}.uns.bam
     
     #Sort BAM file
     samtools sort ${SUBDIR_NAME}.uns.bam -o ${SUBDIR_NAME}.bam
     
     #Index BAM file
     samtools index ${SUBDIR_NAME}.bam
     
     #Remove intermediate files
     rm ${SUBDIR_NAME}.uns.bam
     
     cd ~
	   
	   done
     done
done
