#This script computes percentages starting from ASV count tables containing formatted taxonomic classifications

#Load library
library(data.table)

#Define percentage function
perc <- function(i) {round((i/sum(i))*100, 5)}

#Useful vector
technologies <- c("MiSeq", "AVITI", "NovaSeq")

for (technology in technologies) {

#Load the ASV count table
count_table <- fread(paste0("/path/to/the/", technology, "/starting_table.txt"), header=T, data.table=T)

#Define sample names
samples <- colnames(count_table)[2:37]

#Compute percentages
perc_table <- count_table[, (samples) := lapply(.SD, perc),.SDcols=samples]

#Save the table
write.table(perc_table, paste0("/path/to/the/", technology, "/result.txt"), quote=F, col.names=T, row.names=F, sep="\t")
}