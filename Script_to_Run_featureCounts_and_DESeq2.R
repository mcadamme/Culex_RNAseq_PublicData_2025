#script to count reads using R subread
#MF 06012025



###### This is a description of the files used for comparison ######

# Ryazansky female antennae (n=30): SRR23829117, SRR23829116
# Ryazansky female proboscis (n=30): SRR23829115

# Leal female antennae (n=500): SRR991016
# Leal female hindleg (n=500): SRR991017



###### PREPPING WORKSPACE #######

x <- c("dplyr", "Rsubread", "DESeq2", "ggplot2", "pheatmap","reshape2", "ashr") 
lapply(x, FUN = function(X) {do.call("library", list(X))}) #loading libraries


# Using Rstudio for Windows, made my working directory the one my script is found in. 
setwd("~/OneDrive/Desktop/Menna_Public_RNAseq_Proj_2025")



###### Formatting genome annotation information ######

# loading genome annotation data
z <- read.table("./Annot_with_geneLengths.txt", sep="\t", stringsAsFactors = F)
colnames(z) <- c("GeneID", "Chr", "Annot_type", "Start", "End", "Strand","V9", "for_sorting", "CDS_length", "gene_fam")

head(z)
 

#use these commands to get Gene_ID information from V9 column - will use either this or gene names for DeSeq2 analysis
ids <- gsub(" gene_ID ","",sapply(strsplit(z$V9,";"),.subset,1)) #getting out gene IDs
ids2 <- sapply(strsplit(ids,"="),.subset,2)


#adding ids2 to z for later filtering
z <- data.frame(cbind(z, ids2))

duplicated_values <- z$ids2[duplicated(z$ids2)] 
print(duplicated_values)#shows four instances of "parker23.1095181383,tmenna.1075933603" as gene ID, 
#so need to use the gene names as the "GeneID" in the annotation file,


#generating a simple annotation format file required for feature counting - header names critical
anno <- data.frame(GeneID=as.character(z$GeneID), Chr=z$Chr, Start=z$Start, End=z$End, Strand=z$Strand)
#anno <- data.frame(GeneID=as.character(z$ids2), Chr=z$Chr, Start=z$Start, End=z$End, Strand=z$Strand)



###### Counting reads ######

# counting reads for bam files from each tissue.
SRR2389115<- featureCounts(files = "./Input_BAM_files/SRR23829115.bam", annot.ext=anno, isGTFAnnotationFile=FALSE, 
                      countMultiMappingReads=FALSE, isPairedEnd = TRUE, countReadPairs = TRUE)


# double-checking output to make featureCounts ran properly on first bam file
str(SRR2389115)
SRR2389115$stat  # Summary of assigned/unassigned reads
head(SRR2389115$counts) # Checking head of count matrix (genes x samples)


#running featureCounts on remaining files
SRR2389116<- featureCounts(files = "./Input_BAM_files/SRR23829116.bam", annot.ext=anno, isGTFAnnotationFile=FALSE, 
                           countMultiMappingReads=FALSE, isPairedEnd = TRUE, countReadPairs = TRUE)

SRR2389117<- featureCounts(files = "./Input_BAM_files/SRR23829117.bam", annot.ext=anno, isGTFAnnotationFile=FALSE, 
                           countMultiMappingReads=FALSE, isPairedEnd = TRUE, countReadPairs = TRUE)

SRR991016<- featureCounts(files = "./Input_BAM_files/SRR991016.bam", annot.ext=anno, isGTFAnnotationFile=FALSE, 
                           countMultiMappingReads=FALSE, isPairedEnd = TRUE, countReadPairs = TRUE)

SRR991017<- featureCounts(files = "./Input_BAM_files/SRR991017.bam", annot.ext=anno, isGTFAnnotationFile=FALSE, 
                           countMultiMappingReads=FALSE, isPairedEnd = TRUE, countReadPairs = TRUE)


# this is the raw count table that we need for DESeq2 analysis
tissue_counts <- data.frame(cbind(SRR2389115$counts,SRR2389116$counts,SRR2389117$counts,SRR991016$counts,SRR991017$counts))


# writing out the raw count table file
write.table(tissue_counts, file = "./CQtissue_raw_read_counts.txt")



###### Running differential gene expression analysis ######
#tissue_counts <- read.table("CQtissue_raw_read_counts.txt", header = T) #in case I want to skip featureCounts, which takes 20 min.



# giving the count table row and column names
colnames(tissue_counts) <- c("Prob", "Ant1", "Ant2", "Ant3", "Leg")


# Loading metadata table
metadata <- read.table("./RNAseq_MetaData.txt", row.names = colnames(tissue_counts), header = T)
metadata$Exp <- as.factor(metadata$Exp)
metadata$pool <- as.factor(metadata$pool)
metadata$tissue <- as.factor(metadata$tissue)



# Preparing the DESeq data object
Full_Count_Table <- DESeqDataSetFromMatrix(countData = tissue_counts, colData = metadata, design = ~ Exp + tissue) 

head(Full_Count_Table)

Full_Count_Table$tissue <- relevel(Full_Count_Table$tissue, "Ant") #making antennal data the reference


#Prefilter so only genes with at least 10 reads in at least 1 sample are considered
keep <- rowSums(counts(Full_Count_Table) >= 10) >= 1

Full_Count_Table <- Full_Count_Table[keep,]

dim(Full_Count_Table)#new filtered num genes in dataset and verifying num samples



#Applying the DESeq function to get diff expressed genes.
Full_Count_Table.out <- DESeq(Full_Count_Table)

res<-results(Full_Count_Table.out)

resultsNames(Full_Count_Table.out)

#writing out diff exp results
write.table(res, file = "./tissue_foldChange_and_pValues")

Norm_counts <- counts(Full_Count_Table.out, normalized=TRUE)
Norm_counts_rounded <- round(Norm_counts, digits = 0)
dim(Norm_counts_rounded)

#getting filtered read count data in two formats
dds <- estimateSizeFactors(Full_Count_Table.out)
dim(dds)

filt_locus_names <- data.frame(rownames(dds))
colnames(filt_locus_names) <- c("filt_locus_names")



#getting gene lengths for each locus in filtered dataset
merged_filt_gene_lengths <- merge(filt_locus_names, z, by.x = "filt_locus_names", by.y = "GeneID")
merged_filt_gene_lengths <- merged_filt_gene_lengths[order(merged_filt_gene_lengths$for_sorting), ]#re-sorting by chromosome and position


#checking that the gene names match for dds and merged_filt_gene_lengths

head(rownames(dds))
head(merged_filt_gene_lengths$filt_locus_names)

tail(rownames(dds))
tail(merged_filt_gene_lengths$filt_locus_names)



###### Getting FPKM Values for heatmap ######


counts_matrix <- counts(Full_Count_Table.out)
size_factors <- sizeFactors(Full_Count_Table.out)

gene_lengths_kb <- merged_filt_gene_lengths$CDS_length / 1000

# Normalize counts to library size using size factors
norm_counts <- sweep(counts_matrix, 2, size_factors, FUN = "/")

# Calculate FPKM
fpkm <- sweep(norm_counts, 1, gene_lengths_kb, FUN = "/")

# View result
head(fpkm)



###### Tissue Diff Exp Comp ######

#prob vs. antenna
resLFC <- lfcShrink(dds, contrast = c("tissue", "Ant", "Prob"), type="ashr")

#Look at summary values
summary(resLFC)

#how many significantly DE genes? Looked at adjusted pval of 0.05. Added a FC cutoff for shrunken FCs of 1.5x, as well.
sum(resLFC$padj < 0.05 & abs(resLFC$log2FoldChange) > 0.58, na.rm=TRUE)

#getting full diffexp gene lists with different LFCs and padj - from lfcShrink and results.  
full_genes_Prob_Ant <- data.frame(cbind(resLFC@rownames, resLFC@listData$log2FoldChange,resLFC@listData$padj))
colnames(full_genes_Prob_Ant) <- c("GENE_ID", "LFC", "PADJ")
write.table(full_genes_Prob_Ant, file = "Prob_v_Ant_AllGenes_LFCS.txt", 
            col.names = T, row.names = F, sep = "\t")


#Leg vs. antenna
resLFC2 <- lfcShrink(dds, contrast = c("tissue", "Ant", "Leg"), type="ashr")

#Look at summary values
summary(resLFC2)

#how many significantly DE genes? Looked at adjusted pval of 0.05. Added a FC cutoff for shrunken FCs of 1.5x, as well.
sum(resLFC2$padj < 0.05 & abs(resLFC2$log2FoldChange) > 0.58, na.rm=TRUE)

#getting full diffexp gene lists with different LFCs and padj - from lfcShrink and results.  
full_genes_Leg_Ant <- data.frame(cbind(resLFC2@rownames, resLFC2@listData$log2FoldChange,resLFC2@listData$padj))
colnames(full_genes_Leg_Ant) <- c("GENE_ID", "LFC", "PADJ")
write.table(full_genes_Leg_Ant, file = "Leg_v_Ant_AllGenes_LFCS.txt", 
            col.names = T, row.names = F, sep = "\t")



#Getting lists of diffexp genes in prob-ant and leg-ant comparisons
CS_genes_Prob_Ant <- merge(full_genes_Prob_Ant, z, by.x = "GENE_ID", by.y = "GeneID")
diff_gene_Prob_Ant <- subset(CS_genes_Prob_Ant, PADJ <= 0.05 & gene_fam != "Other")

CS_genes_Leg_Ant <- merge(full_genes_Leg_Ant, z, by.x = "GENE_ID", by.y = "GeneID")
diff_gene_Leg_Ant <- subset(CS_genes_Leg_Ant, PADJ <= 0.05 & gene_fam != "Other")



###### Data visualization ######

#Looking at distr of filtered count data across samples
librarySizes <- colSums(counts(Full_Count_Table))

png(file = "./figures/barplot_libSize.png", units = "px", height = 400, width = 600)
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, ylim = c(0,1.5e+08),
        main="Barplot of library sizes")

dev.off()

logcounts <- log2(counts(Full_Count_Table) + 1)
head(logcounts)


#Visualize difference between per gene counts for each of the experiments
statusCol <- as.numeric(factor(Full_Count_Table$Exp)) + 1  #make a color vector by experiment

png(file = "./figures/boxplot_pergenecounts.png", units = "px", height = 400, width = 600)
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)

abline(h=median(as.matrix(logcounts)), col="blue") #Adding median log counts
dev.off()



#Heatmap showing correlation coefficients for gene expression values between tissues.
cormat <- cor(Norm_counts, method = "spearman")#

#Function to Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)# Melt the correlation matrix

# Create a ggheatmap
png(file = "./figures/overall_heatmap_plot.png", units = "px", height = 400, width = 600)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "pink", 
                       midpoint = 0.85, limit = c(0.7,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  xlab("") + ylab("")+
  coord_fixed()

print(ggheatmap)
dev.off()


#plotting overall diff exp in 2 ways
png(file = "./figures/MA_plot.png", units = "px", height = 400, width = 600)
plotMA (res, ylim = c(-2, 2), cex = 1.5 )
dev.off()

png(file = "./figures/bygene_heatmap_plot.png", units = "px", height = 400, width = 600)
heatmap(as.matrix(tissue_counts), Rowv = NA, Colv=NA, scale = "row" )
dev.off()



#plotting expression of our gene families of interest

log_fpkm <- log2(fpkm + 1)

#subsetting by gene family
gene_fam <- data.frame(GeneID=as.character(merged_filt_gene_lengths$filt_locus_names), GeneFam=merged_filt_gene_lengths$gene_fam)
Or_fam <- subset(gene_fam, GeneFam == "OR")
Ir_fam <- subset(gene_fam, GeneFam == "IR")
Obp_fam <- subset(gene_fam, GeneFam == "OBP")
Gr_fam <- subset(gene_fam, GeneFam == "GR")


# OR heatmap
Or_fpkm <- log_fpkm[rownames(log_fpkm) %in% Or_fam$GeneID, ]

png(file = "./figures/Or_heatmap.png", units = "px", height = 1000, width = 800)
pheatmap(Or_fpkm,
         cluster_cols = TRUE,
         #scale = "row",  # normalize rows (genes)
         fontsize_col = 10,
         fontsize_row = 10,
         show_rownames = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()


# IR heatmap
Ir_fpkm <- log_fpkm[rownames(log_fpkm) %in% Ir_fam$GeneID, ]

png(file = "./figures/Ir_heatmap.png", units = "px", height = 1200, width = 800)
pheatmap(Ir_fpkm,
         cluster_cols = TRUE,
         #scale = "row",  # normalize rows (genes)
         fontsize_col = 10,
         fontsize_row = 10,
         show_rownames = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()


# OBP heatmap
Obp_fpkm <- log_fpkm[rownames(log_fpkm) %in% Obp_fam$GeneID, ]

png(file = "./figures/Obp_heatmap.png", units = "px", height = 1000, width = 800)
pheatmap(Obp_fpkm,
         cluster_cols = TRUE,
         #scale = "row",  # normalize rows (genes)
         fontsize_col = 10,
         fontsize_row = 10,
         show_rownames = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()


# GR heatmap
Gr_fpkm <- log_fpkm[rownames(log_fpkm) %in% Gr_fam$GeneID, ]

png(file = "./figures/Gr_heatmap.png", units = "px", height = 1000, width = 800)
pheatmap(Gr_fpkm,
         cluster_cols = TRUE,
         #scale = "row",  # normalize rows (genes)
         fontsize_col = 10,
         fontsize_row = 10,
         show_rownames = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()

