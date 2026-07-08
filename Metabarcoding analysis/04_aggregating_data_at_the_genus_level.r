#This script aggregates ASV counts or percentages at the genus level starting from ASV tables containing formatted taxonomic classifications

#Load library
library(data.table)


#Use ASV count tables
count_tech <- list(a=list("path/to/MiSeq/ASV/count/table.txt", "MiSeq"),
b=list("path/to/AVITI/ASV/count/table.txt", "AVITI"),
d=list("path/to/NovaSeq/ASV/count/table.txt", "NovaSeq"))

for (count_technology in count_tech) {

#Load data
count_table <- fread(count_technology[[1]], header=T, data.table=T)

#Define sample names
samples_count <- colnames(count_table)[2:37]

#Discard NAs
count_table <- count_table[!is.na(Genus)]

#Sum all the counts that belong to the same genus within each sample 
count_table_2 <- count_table[,lapply(.SD,sum), by=Genus, .SDcols=samples_count] 

#Save table
write.table(count_table_2, paste0("/path/to/the/", count_technology[[2]], "/result.txt"), quote=F, sep="\t", row.names=F, col.names=T)
}


#Use ASV percentage tables
perc_tech <- list(a=list("path/to/MiSeq/ASV/percentage/table.txt", "MiSeq"),
b=list("path/to/AVITI/ASV/percentage/table.txt", "AVITI"),
d=list("path/to/NovaSeq/ASV/percentage/table.txt", "NovaSeq"))

for (perc_technology in perc_tech) {

#Load data
perc_table <- fread(perc_technology[[1]], header=T, data.table=T)

#Define sample names
samples_perc <- colnames(perc_table)[2:37]

#Discard NAs
perc_table <- perc_table[!is.na(Genus)]

#Sum all the percentages that belong to the same genus within each sample 
perc_table_2 <- perc_table[,lapply(.SD,sum), by=Genus, .SDcols=samples_perc] 

#Save table
write.table(perc_table_2, paste0("/path/to/the/", perc_technology[[2]], "/result.txt"), quote=F, sep="\t", row.names=F, col.names=T)
}