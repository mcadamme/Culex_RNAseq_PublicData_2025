#! /bin/bash

#This is my script to run STAR v. 2.7.11b
#MLF 05262025


##Generate Genome Index
/home/megan/src/STAR-2.7.11b/source/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /home/megan/Desktop/Public_mosquito_rnaseq/Genome \
--genomeFastaFiles /home/megan/Desktop/Public_mosquito_rnaseq/Genome/GCA_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna \
--sjdbGTFfile /home/megan/Desktop/Public_mosquito_rnaseq/Genome/merged_sorted_locusName_06042025.gff3 --sjdbOverhang 150 --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 13.5


# Setting the base directory for read mapping
BASE_DIR="/home/megan/Desktop/Public_mosquito_rnaseq"

for first_level in "$BASE_DIR"/SRR*/; do

    for second_level in "$first_level"*/; do
    
    	#getting SRR name
	SUBDIR_NAME=$(basename "$(dirname "${second_level}")")
               
        # Define inputput file names
    	INFILE1="${second_level}${SUBDIR_NAME}"_trimmed_R1_paired.fastq.gz
    	INFILE2="${second_level}${SUBDIR_NAME}"_trimmed_R2_paired.fastq.gz
    
	#Run Mapping Job with ENCODE standard options
	/home/megan/src/STAR-2.7.11b/source/STAR --runThreadN 7 --genomeDir /home/megan/Desktop/Public_mosquito_rnaseq/Genome \
	--readFilesCommand gunzip -c --readFilesIn "${INFILE1}" "${INFILE2}" --outFilterType BySJout --outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 \
	--alignIntronMax 1000000 --alignMatesGapMax 1000000

	mv Aligned.out.sam "${second_level}${SUBDIR_NAME}".sam
	cat Log.final.out >> Log_summary
	
	done
done

cp ~/Log_summary "$BASE_DIR"

