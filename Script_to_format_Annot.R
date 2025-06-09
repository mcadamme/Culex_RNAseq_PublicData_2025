#This is the script I used to clean and format my annotation file prior to featureCounts
#MF 06052025


###### Loading libraries and setting up workspace ######

x <- c("dplyr","reshape2", "stringr") 
lapply(x, FUN = function(X) {do.call("library", list(X))}) #loading libraries


# Working directory 
setwd("~/OneDrive/Desktop/Menna_Public_RNAseq_Proj_2025")


# Loading data

# Genome annotation data
gff_full <- read.table("./genome_dir/merged_sorted_locusName_06042025.gff3", sep="\t", quote="", stringsAsFactors = F)
head(gff_full)

# CDS length data
for_fpkm <- read.table("Cqui_CDS_lengths.txt", header = F, stringsAsFactors = F)




###### Approach 1 - Subsetting gff3 on gene ######


gff <- gff_full[gff_full$V3=="gene",] #subsetting annotation to pull only gene 
head(gff)


gff$locus_id <- str_extract(gff$V9, "(?<=locus=)[^ ;]+")#this gets the locus name

gff$for_sorting <- 1:nrow(gff)#creates a numeric vector for sorting later



# Define the mapping of search strings to factor labels
patterns <- c("Or\\d{1,3}", "Ir\\d{1,3}", "Gr\\d{1,3}", "Obp\\d{1,3}")
labels <- c("OR", "IR", "GR", "OBP")

# Initialize new column with "Other"
gff$gene_fam <- "Other"




# Assign labels based on pattern matches
for (i in seq_along(patterns)) {
  match_indices <- grepl(patterns[i], gff$locus_id, ignore.case = FALSE)
  gff$gene_fam[match_indices] <- labels[i]
}

gff$gene_fam <- factor(gff$gene_fam, levels = c(labels, "Other"))

nrow(gff) #gives 15284


###### Approach 2 - filtering CDS count file by gene name ######

for_fpkm$locus_id <- str_extract(for_fpkm$V1, "(?<=locus=)[^ ;]+")

# Initialize new column with "Other"
for_fpkm$gene_fam <- "Other"


# Assign labels based on pattern matches from above
for (i in seq_along(patterns)) {
  match_indices <- grepl(patterns[i], for_fpkm$locus_id, ignore.case = FALSE)
  for_fpkm$gene_fam[match_indices] <- labels[i]
}

for_fpkm$gene_fam <- factor(for_fpkm$gene_fam, levels = c(labels, "Other"))


for_fpkm_filtered <- for_fpkm %>%
  group_by(locus_id) %>%
  slice_max(order_by = V3, n = 1, with_ties = FALSE) %>%
  ungroup()

nrow(for_fpkm_filtered) #gives 15283

for_fpkm_filtered <- data.frame(for_fpkm_filtered)




###### Difference in "gene names" between these two approaches ######

# Find matches (shared elements)
matches <- intersect(gff$locus_id, for_fpkm_filtered$locus_id)
length(matches)

# Find differences (unique elements)
only_in_gff <- setdiff(gff$locus_id, for_fpkm_filtered$locus_id)
only_in_gff

only_in_fpkm <- setdiff(for_fpkm_filtered$locus_id, gff$locus_id)
only_in_fpkm

#neither setdiff shows a line that could be considered something that differed

#looking at structure of the two files
str(gff)
str(for_fpkm_filtered)

#checking for duplicates in gff
unique_duplicates <- unique(gff$locus_id[duplicated(gff$locus_id)])
print(unique_duplicates)

#This shows that the gff3 contains two genes labeled OBP62, but one should be OBP125.  
#We had previously asked VeuPathDB to fix this, but they have not. 
#The one at 131771762 should be obp125 and the one at 131778146 should be obp62



###### Here is the correction ######

rows_with_Obp62 <- grep("Obp62", gff$locus_id)
print(rows_with_Obp62)

print(gff[c(13164,13166),])


gff[13164,10] = "Obp125"
for_fpkm[21239,4] = "Obp125"


#re-filtering fpkm file
for_fpkm_filtered <- for_fpkm %>%
  group_by(locus_id) %>%
  slice_max(order_by = V3, n = 1, with_ties = FALSE) %>%
  ungroup()

nrow(for_fpkm_filtered) #now gives 15284

for_fpkm_filtered <- data.frame(for_fpkm_filtered)


###### Generating merged full dataset ######

sub_for_fpkm <- for_fpkm_filtered[,c(3:5)]

merged_annot <- merge(gff, sub_for_fpkm, by = "locus_id")

row_with_Obp125 <- grep("Obp125", merged_annot$locus_id) #checking that the obp renaming worked
print(row_with_Obp125)


#prepping for write out
merged_annot <- merged_annot[order(merged_annot$for_sorting), ]#re-sorting by chromosome and position
sub_merged_annot <- merged_annot[,c(-3,-7,-9,-12)]



write.table(sub_merged_annot, file = "./Annot_with_geneLengths.txt", sep="\t", col.names=FALSE, row.names=FALSE)



