#! /bin/bash

#This is the script I used to run fastqc on the RNAseq data from the SRA
#MLF 05262025

for file in ~/Desktop/Public_mosquito_rnaseq/SRR*/fastq

do

	/home/megan/src/FastQC/fastqc *1.fastq.gz
	/home/megan/src/FastQC/fastqc *2.fastq.gz

done

