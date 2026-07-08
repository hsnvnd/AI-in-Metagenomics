#This script formats SILVA taxonomic classifications in raw ASV count tables
#"Incertae Sedis" and "Unknown Family" classifications are converted to NA
#Underscores are replaced with single spaces
# After reformatting:
# - Each taxonomic classification consists of one or more words separated by a single space
# - NA is the only value used for unclassified taxonomic levels
# - Hyphens (-) are preserved as part of taxonomic classification names

#Load libraries
library(data.table)
library(stringr)

#Useful vectors
technologies <- c("Miseq", "Aviti", "Novaseq")
level <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

for (technology in technologies){

#Load raw ASV count table
raw_table <- fread(paste0("/path/to/the/table/", technology, "/path/to/the/table.txt"), header=T, data.table=T)

#Convert taxonomic classifications containing "Incertae Sedis" or "Unknown Family" to NA
reformatted_table <- raw_table[, (level) := lapply(.SD, function(x){x[grepl("Unknown|Incertae", x)] <- NA; x}), .SDcols=level]

#Convert each underscore into a single space
reformatted_table_2 <- reformatted_table[, (level) := lapply(.SD, function(column) str_replace_all(column, "_", " ")), .SDcols=level]

#Save the table
write.table(reformatted_table_2, paste0("/path/to/the/result/", technology, "/path/to/the/result.txt"), quote=F, col.names=T, row.names=F, sep="\t")
}