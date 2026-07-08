#This script uses the median-of-ratios method implemented in the DESeq2 R library to normalize genus count data

#Load libraries
library(data.table)
library("DESeq2")
library(stringr)

#Variables to loop on 
technologies <- list(a=list("path/to/MiSeq/genus/count/table.txt", "MiSeq"), 
b=list("path/to/AVITI/genus/count/table.txt", "AVITI"),
d=list("path/to/NovaSeq/genus/count/table.txt", "NovaSeq"))

#Load metadata
metadata <- fread("/path/to/metadata.txt", header=T, data.table=T)

#Remove the "T1_" from sample names
metadata[, Sample_name := str_remove(Sample_name, "T1_")]

#Merge the 3 treatment columns
metadata[, Soil_Autoclave_Heat_root := paste(Soil, Autoclave, Heat_root, sep="_")]
metadata[, c("Soil", "Autoclave", "Heat_root") := NULL]

#Transform metadata object in a dataframe
metadata_df <- as.data.frame(metadata)
rownames(metadata_df) <- metadata_df$Sample_name
metadata_df$Sample_name <- NULL
metadata_df$Soil_Autoclave_Heat_root <- as.factor(metadata_df$Soil_Autoclave_Heat_root)

for (tech in technologies) {

#Load the table
table_tech <- fread(tech[[1]], header=T, data.table=T)

#Reorder columns
setcolorder(table_tech, c("Genus", as.character(1:36)))

#Transform the table in a dataframe
table_tech_df <- as.data.frame(table_tech)

#Set genera as rownames
rownames(table_tech_df) <- table_tech_df$Genus
table_tech_df$Genus <- NULL

#Check consistency between colnames in data table and rownames in metadata
# > all(rownames(metadata_df) == colnames(table_tech_df))
# [1] TRUE

#Build a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=table_tech_df, colData=metadata_df, design= ~ Soil_Autoclave_Heat_root)

#No pre-filtering step is performed

#Run DESeq function
dds_2 <- DESeq(dds)

#Retrieve normalized counts
tech_normalized <- counts(dds_2, normalized=TRUE)

#Refine table
tech_normalized_df <- as.data.frame(tech_normalized)
tech_normalized_df$Genus <- rownames(tech_normalized_df)
tech_normalized_dt <- as.data.table(tech_normalized_df)
setcolorder(tech_normalized_dt, c("Genus", as.character(1:36)))

#Save normalized counts
write.table(tech_normalized_dt, paste0("/path/to/the/", tech[[2]], "/result.txt"), quote=F, sep="\t", col.names=T, row.names=F) 
}