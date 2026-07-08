#This script formats SILVA taxonomic classifications in ASV count tables
#"Incertae Sedis" and "Unknown Family" classifications are converted to NA
#Underscores are replaced with single spaces
# After formatting:
# - Each taxonomic classification consists of one or more words separated by a single space
# - NA is the only value used for unclassified taxonomic levels
# - Hyphens (-) are preserved as part of taxonomic classification names

#Load libraries
library(data.table)
library(stringr)

#Useful vectors
technologies <- c("MiSeq", "AVITI", "NovaSeq")
level <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

for (technology in technologies){

#Load ASV count table
starting_table <- fread(paste0("/path/to/the/", technology, "/starting_table.txt"), header=T, data.table=T)

#Convert taxonomic classifications containing "Incertae Sedis" or "Unknown Family" to NA
formatted_table <- starting_table[, (level) := lapply(.SD, function(x){x[grepl("Unknown|Incertae", x)] <- NA; x}), .SDcols=level]

#Convert each underscore into a single space
formatted_table_2 <- formatted_table[, (level) := lapply(.SD, function(column) str_replace_all(column, "_", " ")), .SDcols=level]

#Save the table
write.table(formatted_table_2, paste0("/path/to/the/", technology, "/result.txt"), quote=F, col.names=T, row.names=F, sep="\t")
}