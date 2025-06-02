#script to count reads using R subread
#MF 06012025



###### This is a description of the files used for comparison ######

# Ryazansky female antennae (n=30): SRR23829117, SRR23829116
# Ryazansky female proboscis (n=30): SRR23829115

# Leal female antennae (n=500): SRR991016
# Leal female hindleg (n=500): SRR991017



####### Loading libraries ######

library(Rsubread); library(DESeq2); library(pheatmap); library(dplyr)



###### Setting working directory ######

#Using Rstudio for Windows, made my working directory the one my script is found in. 
setwd("~/OneDrive/Desktop/Menna_Public_RNAseq_Proj_2025")



###### Formatting genome annotation information ######

#loading genome annotation data
z <- read.table("./genome_dir/merged_sorted_locusName.gff3", sep="\t", quote="", stringsAsFactors = F)
head(z)

z <- z[z$V3=="gene",] #subsetting to pull only gene 
head(z)# we can check subset worked using the head command

ids <- gsub(" gene_ID ","",sapply(strsplit(z$V9,";"),.subset,1)) #getting out gene IDs
ids2 <- sapply(strsplit(ids,"="),.subset,2)

names <- gsub(" gene_name ","",sapply(strsplit(z$V9,";"),.subset,2)) #getting out gene names
names2 <- sapply(strsplit(names,"="),.subset,2)

#generating dataframe that matches gene ids and names
id_name_match <- data.frame(GeneID=as.character(ids2),GeneName=as.character(names2))
head(id_name_match)

#searching id_name_match for specific gene families and assigning them to different categories

# Define the mapping of search strings to factor labels
patterns <- c("Or\\d{1,3}", "Ir\\d{1,3}", "Gr\\d{1,3}", "Obp\\d{1,3}")
labels <- c("OR", "IR", "GR", "OBP")

# Initialize new column with "Other"
id_name_match$gene_fam <- "Other"


# Assign labels based on pattern matches
for (i in seq_along(patterns)) {
  match_indices <- grepl(patterns[i], id_name_match$GeneName, ignore.case = FALSE)
  id_name_match$gene_fam[match_indices] <- labels[i]
}

id_name_match$gene_fam <- factor(id_name_match$gene_fam, levels = c(labels, "Other"))
head(id_name_match)
str(id_name_match)


#generating a simple annotation format file required for feature counting - header names critical
anno <- data.frame(GeneID=as.character(ids2), Chr=z$V1, Start=z$V4, End=z$V5, Strand=z$V7)


###### Counting reads ######

#counting reads for bam files from each tissue.
SRR2389115<- featureCounts(files = "./Input_BAM_files/SRR23829115.bam", annot.ext=anno, isGTFAnnotationFile=FALSE, 
                      countMultiMappingReads=FALSE, isPairedEnd = TRUE, countReadPairs = TRUE)

SRR2389116<- featureCounts(files = "./Input_BAM_files/SRR23829116.bam", annot.ext=anno, isGTFAnnotationFile=FALSE, 
                           countMultiMappingReads=FALSE, isPairedEnd = TRUE, countReadPairs = TRUE)

SRR2389117<- featureCounts(files = "./Input_BAM_files/SRR23829117.bam", annot.ext=anno, isGTFAnnotationFile=FALSE, 
                           countMultiMappingReads=FALSE, isPairedEnd = TRUE, countReadPairs = TRUE)

SRR991016<- featureCounts(files = "./Input_BAM_files/SRR991016.bam", annot.ext=anno, isGTFAnnotationFile=FALSE, 
                           countMultiMappingReads=FALSE, isPairedEnd = TRUE, countReadPairs = TRUE)

SRR991017<- featureCounts(files = "./Input_BAM_files/SRR991017.bam", annot.ext=anno, isGTFAnnotationFile=FALSE, 
                           countMultiMappingReads=FALSE, isPairedEnd = TRUE, countReadPairs = TRUE)


#this is the raw count table that we need for DESeq2 analysis
tissue_counts <- data.frame(cbind(SRR2389115$counts,SRR2389116$counts,SRR2389117$counts,SRR991016$counts,SRR991017$counts))


#giving the count table row and column names
rownames(tissue_counts) <- ids2
colnames(tissue_counts) <- c("Prob", "Ant1", "Ant2", "Ant3", "Leg")


#writing out the raw count table file
write.table(tissue_counts, file = "./CQtissue_raw_read_counts.txt")



###### Running differential gene expression analysis ######

#Loading metadata table
metadata <- read.table("./RNAseq_MetaData.txt", row.names = colnames(tissue_counts), header = T)
metadata$Exp <- as.factor(metadata$Exp)
metadata$pool <- as.factor(metadata$pool)
metadata$tissue <- as.factor(metadata$tissue)


#Preparing the DESeq data object
Full_Count_Table <- DESeqDataSetFromMatrix(countData = tissue_counts, colData = metadata, design = ~ Exp + tissue) 

head(Full_Count_Table)

Full_Count_Table$tissue <- relevel( Full_Count_Table$tissue, "Ant")#making antennal data the reference

#Applying the DESeq function to get diff expressed genes.
Full_Count_Table.out <- DESeq(Full_Count_Table)

#writing out diff exp results
diff_exp_results <- results(Full_Count_Table.out)

write.table(diff_exp_results, file = "./tissue_foldChange_and_pValues")



###### Prepping data for visualization ######

#plotting overall diff exp in 2 ways
png(file = "./MA_plot.png", units = "px", height = 400, width = 600)
plotMA (diff_exp_results, ylim = c(-2, 2), cex = 1.5 )
dev.off()

png(file = "./heatmap_plot.png", units = "px", height = 400, width = 600)
heatmap(as.matrix(tissue_counts), Rowv = NA, Colv=NA, scale = "row" )
dev.off()


#plotting expression of our gene families of interest
#Get normalized counts
normalized_counts <- counts(Full_Count_Table.out, normalized = TRUE)
rld <- rlog(normalized_counts)
rld2 <- assay(rld)
normalized_counts_withID <- data.frame(cbind(id = as.character(rownames(rld2)), rld2))
rownames(normalized_counts_withID) <- NULL
normalized_counts_withID$id <- as.character(normalized_counts_withID$id)
str(normalized_counts_withID)


#Converting read count cols to numeric
cols_to_convert <- c("Prob", "Ant1", "Ant2", "Ant3", "Leg")

# Convert each column in-place
for (col in cols_to_convert) {
  normalized_counts_withID[[col]] <- as.numeric(normalized_counts_withID[[col]])
}

str(normalized_counts_withID)


#acquiring gene families of interest
data_for_split <- normalized_counts_withID %>%
  left_join(id_name_match, by = c("id" = "GeneID"))


data_after_split <- subset(data_for_split, gene_fam != "Other")
data_after_split$id <- as.character(data_after_split$id)
str(data_after_split)

table(data_after_split$gene_fam) #sanity check to make sure the numbers of genes look reasonable.


#looks like a lot of genes have 0 for normalized expression values - how many?
data_after_split$exp_sums <- rowSums(data_after_split[,c(2:6)])

noexp_data_after_split <- subset(data_after_split, exp_sums == 0)
table(noexp_data_after_split$gene_fam)

exp_data_after_split <- subset(data_after_split, exp_sums != 0)
table(exp_data_after_split$gene_fam)


mat <- as.matrix(exp_data_after_split[, c("Prob", "Ant1", "Ant2", "Ant3", "Leg")])
rownames(mat) <- exp_data_after_split$GeneName #can switch to gene name if we want


#calibrating - this function implements row-wise Z-score normalization. This transforms the values to Z-scores, 
#where transformed vector has mean = 0 and a standard deviation = 1.

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}


#heatmap plot of chemosensory gene families only
png(file = "./Chemosens_heatmap_plot.png", units = "px", height = 4000, width = 6000)
data_subset_norm <- t(apply(mat, 1, cal_z_score))
pheatmap(data_subset_norm,
         fontsize = 30,
         fontsize_row = 15,
         fontsize_col = 50)
dev.off()


