#This script is used to denoise sequencing data using the dada2 R library: ASVs are inferred and then assigned taxonomic classifications based on the SILVA database
#This script follows the dada2 R library tutorial: https://benjjneb.github.io/dada2/tutorial.html 
#Before running this script, set the following paths: the path to the reads (line 12), the path to the output directory where the results will be stored (line 14), and the path to the SILVA database (line 294)
#Use the --help flag for usage information


suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-R", "--readsdir"), type="character", default="path/to/sequencing/data",
              help="Input file", metavar="character"),
  make_option(c("-B", "--basedir"), type="character", default="path/to/the/results/directory/",
              help="Input file", metavar="character"),
  make_option(c("-P", "--previously_run"), type="logical", default=FALSE,
              help="Was read trimming and error estimate previously run?", metavar="character"),
  make_option(c("-S", "--use_silva"), type="logical", default=TRUE,
              help="Use SILVA database for classification?", metavar="character"),
  make_option(c("-F", "--filepattern"), type="character", default="",
              help="Input file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$basedir)) {
  stop("ERROR: No basedir specified with '-I' flag.")
} else {  cat ("basedir is ", opt$basedir, "\n")
  basedir <- opt$basedir  
  }

if (is.null(opt$readsdir)) {
  stop("ERROR: No readsdir specified with '-I' flag.")
} else {  cat ("readsdir is ", opt$readsdir, "\n")
  readsdir <- opt$readsdir  
  }

if (is.null(opt$previously_run)) {
  stop("ERROR: No previously_run specified with '-R' flag.")
} else {  cat ("previously_run is ", opt$previously_run, "\n")
  previously_run <- opt$previously_run  
  }

if (is.null(opt$use_silva)) {
  stop("ERROR: No use_silva specified with '-S' flag.")
} else {  cat ("use_silva is ", opt$use_silva, "\n")
  use_silva <- opt$use_silva  
  }

if (is.null(opt$filepattern)) {
  stop("ERROR: No filepattern specified with '-I' flag.")
} else {  cat ("filepattern is ", opt$filepattern, "\n")
  filepattern <- opt$filepattern  
  }

#Load libraries
library(dada2)
library(ShortRead)
library(Biostrings)
library(data.table)
library(stringr)

#Set the number of cores
ncores=8

if(!previously_run)
{

#Set read paths 
fnFs <- list.files(readsdir, pattern="R1_001.fastq.gz", full.names=TRUE) 
fnRs <- list.files(readsdir, pattern="R2_001.fastq.gz", full.names=TRUE)

#Select some samples if needed
fnFs <- fnFs[grep(filepattern, fnFs)]
fnRs <- fnRs[grep(filepattern, fnRs)]
 
#Extract sample names
find_strings <- c("Aviti", "Novaseq")
if (unique(apply(sapply(find_strings , grepl, fnFs), 1, any))){
sample.names <- str_match(fnFs, "37_\\s*(.*?)\\s*-")[,2]
}else{
sample.names <- str_match(fnFs, "46-\\s*(.*?)\\s*-")[,2]
}

#Set primers
FWD <- "CCTACGGGNBGCASCAG"
REV <- "GACTACNVGGGTATCTAATCC"

#Create all orientation sequences for the primers
allOrients <- function(primer) {
    require(Biostrings)
    dna <- DNAString(primer)
    orients <- c(Forward=dna, Complement=complement(dna), Reverse=reverse(dna), RevComp=reverseComplement(dna))
    return(sapply(orients, toString))}

#Compute all orientations for the forward primer
FWD.orients <- allOrients(FWD) 

#Compute all orientations for the reverse primer
REV.orients <- allOrients(REV)

#Remove reads containing Ns
fnFs.filtN <- file.path(paste(basedir, "01_reads_without_Ns/", basename(fnFs), sep="")) 
fnRs.filtN <- file.path(paste(basedir, "01_reads_without_Ns/", basename(fnRs), sep=""))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN=0, truncLen=0, compress=T, multithread=ncores)

#Count primers found in the reads using all the possible primer orientations
#Use only the first sample since all the samples have the same library preparation
primerHits <- function(primer, fn) {
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed=FALSE)
    return(sum(nhits > 0))}

primer_hit_before <- rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fnFs.filtN[[1]]), 
    FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fnRs.filtN[[1]]), 
    REV.ForwardReads=sapply(REV.orients, primerHits, fn=fnFs.filtN[[1]]), 
    REV.ReverseReads=sapply(REV.orients, primerHits, fn=fnRs.filtN[[1]]))

print(primer_hit_before)

#Table example from MiSeq data
                 # Forward Complement Reverse RevComp
# FWD.ForwardReads   34095          0       0       0
# FWD.ReverseReads      22          0       0       1
# REV.ForwardReads      25          0       0      24
# REV.ReverseReads   43965          0       0       0

#Remove primers
fnFs.filtN_primers <- file.path(paste(basedir, "02_trimmed_reads_without_N/", basename(fnFs), sep=""))
fnRs.filtN_primers <- file.path(paste(basedir, "02_trimmed_reads_without_N/", basename(fnRs), sep=""))
filterAndTrim(fnFs, fnFs.filtN_primers, fnRs, fnRs.filtN_primers, maxN=0, trimLeft=c(17,21), compress=T, multithread=ncores)

#Count primers again
primer_hit_after <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn=fnFs.filtN_primers[[1]]), 
    FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fnRs.filtN_primers[[1]]), 
    REV.ForwardReads=sapply(REV.orients, primerHits, fn=fnFs.filtN_primers[[1]]), 
    REV.ReverseReads=sapply(REV.orients, primerHits, fn=fnRs.filtN_primers[[1]]))

print(primer_hit_after)

#Table example from MiSeq data
                 # Forward Complement Reverse RevComp
# FWD.ForwardReads      17          0       0       0
# FWD.ReverseReads       0          0       0       1
# REV.ForwardReads       0          0       0      24
# REV.ReverseReads       0          0       0       0

#Plot base quality
dir.create(paste0(basedir, "04_plots/"))
pdf(paste0(basedir, "04_plots/", filepattern, "_quality_cut_F.pdf"))
base_quality_F <- plotQualityProfile(fnFs.filtN_primers)
print(base_quality_F)
dev.off()

pdf(paste0(basedir, "04_plots/", filepattern, "_quality_cut_R.pdf"))
base_quality_R <- plotQualityProfile(fnRs.filtN_primers)
print(base_quality_R)
dev.off()

#Trim low quality bases at the 3' end of the reads
fnFs.trimmed <- file.path(paste(basedir, "03_quality_trimmed_reads/", basename(fnFs), sep="")) 
fnRs.trimmed <- file.path(paste(basedir, "03_quality_trimmed_reads/", basename(fnRs), sep="")) 
out <- filterAndTrim(fnFs, fnFs.trimmed, fnRs, fnRs.trimmed, maxN=0, truncLen=c(250,230), trimLeft=c(17,21), compress=T, multithread=ncores)

#Plot base quality
pdf(paste0(basedir, "04_plots/", filepattern, "_quality_nolow_F.pdf"))
base_quality_F <- plotQualityProfile(fnFs.trimmed)
print(base_quality_F)
dev.off()

pdf(paste0(basedir, "04_plots/", filepattern, "_quality_nolow_R.pdf"))
base_quality_R <- plotQualityProfile(fnRs.trimmed)
print(base_quality_R)
dev.off()

#Estimate error rates
errF <- learnErrors(fnFs.trimmed, multithread=ncores)
errR <- learnErrors(fnRs.trimmed, multithread=ncores)

#Plot error rates
pdf(paste0(basedir, "04_plots/", filepattern, "_errors_nolow_F.pdf"))
error_rate_F <- plotErrors(errF, nominalQ=TRUE)
print(error_rate_F)
dev.off()

pdf(paste0(basedir, "04_plots/", filepattern, "_errors_nolow_R.pdf"))
error_rate_R <- plotErrors(errR, nominalQ=TRUE)
print(error_rate_R)
dev.off()

#Sample inference
dadaFs <- dada(fnFs.trimmed, err=errF, multithread=ncores)
dadaRs <- dada(fnRs.trimmed, err=errR, multithread=ncores)

#Merge overlapping reads
mergers <- mergePairs(dadaFs, fnFs.trimmed, dadaRs, fnRs.trimmed, verbose=TRUE, returnRejects=FALSE)

#Build ASV table
merged_seqtab <- makeSequenceTable(mergers)

#Remove chimeras
merged_seqtab.nochim <- removeBimeraDenovo(merged_seqtab, method="consensus", multithread=ncores, verbose=TRUE)

#Save the table. This table is also used to verify whether dada2 has been previously run.
merged_nochim <- t(merged_seqtab.nochim)
dir.create(paste0(basedir, "05_tables/"))
write.table(merged_nochim, paste0(basedir, "05_tables/", filepattern, "_merged_nochim.txt"), sep="\t", quote=F, row.names=T)
merged_nochim <- data.frame(cbind(row.names(merged_nochim), merged_nochim))
names(merged_nochim) <- c("seq", sample.names)

#Track how many reads passed the quality filtering steps
getN <- function(x) sum(getUniques(x))
 
if(ncol(merged_nochim)>1)
{

#Multi sample
track <- data.table(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(merged_seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
perc_track<-track[, .(perc_input=100, perc_filtered=round(100*filtered/input, 2), perc_denoisedF=round(100*denoisedF/input, 2), perc_denoisedR=round(100*denoisedR/input, 2), perc_merged=round(100*merged/input, 2), perc_nonchim=round(100*nonchim/input, 2))]
rownames(perc_track) <- sample.names

}else{

#Single sample
track <- data.table(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(merged_seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
perc_track<-track[, .(perc_input=100, perc_filtered=round(100*filtered/input, 2), perc_denoisedF=round(100*denoisedF/input, 2), perc_denoisedR=round(100*denoisedR/input, 2), perc_merged=round(100*merged/input, 2), perc_nonchim=round(100*nonchim/input, 2))]
rownames(perc_track) <- sample.names
}

#Save tables
write.table(track, paste0(basedir, "05_tables/", filepattern, "_read_numbers_merged.txt"), row.names=T)
write.table(perc_track, paste0(basedir, "05_tables/", filepattern, "_perc_read_numbers_merged.txt"), row.names=T)
}


#If dada2 has already been run, it is possible to skip the section above
if(previously_run)
{

#Load libraries
library(dada2)
library(ShortRead)
library(Biostrings)
library(data.table)
library(stringr)

#Set read paths
fnFs <- list.files(readsdir, pattern="R1_001.fastq.gz", full.names=TRUE)
fnRs <- list.files(readsdir, pattern="R2_001.fastq.gz", full.names=TRUE)

#Select some samples if needed
fnFs <- fnFs[grep(filepattern, fnFs)]
fnRs <- fnRs[grep(filepattern, fnRs)]

#Extract sample names
find_strings <- c("Aviti", "Novaseq")
if (unique(apply(sapply(find_strings , grepl, fnFs), 1, any))){
sample.names <- str_match(fnFs, "37_\\s*(.*?)\\s*-")[,2]
}else{
sample.names <- str_match(fnFs, "46-\\s*(.*?)\\s*-")[,2]
}

#Read the dataframe with ASV on rows and samples on columns
merged_nochim <- fread(paste0(basedir, "/05_tables/", filepattern, "_merged_nochim.txt"), data.table=F)
row.names(merged_nochim) <- merged_nochim$V1 
merged_nochim$V1 <- NULL
merged_seqtab.nochim <- t(merged_nochim) 
merged_nochim <- data.frame(cbind(row.names(merged_nochim), merged_nochim)) 
names(merged_nochim) <- c("seq", sample.names) 
}


#Classify reads with SILVA
if(use_silva)
{

#Perform taxonomic classification
dir.create(paste0(basedir, "/06_SILVA_tables/"))
myfasta_silva <- "path/to/SILVA/database"
merged_taxa <- assignTaxonomy(merged_seqtab.nochim, myfasta_silva, minBoot=50, multithread=ncores)

#According to dada2 the first column contains the Kingdom level. This is not true as the first column contains the Domain level. Hence the Kingdom level is replaced with the Domain level.
colnames(merged_taxa)[1] <- "Domain"
write.table(merged_taxa, paste0(basedir, "/06_SILVA_tables/", filepattern, "_Taxa_merged.txt"), quote=F, col.names=NA)

#Merge the ASV abundance table with the taxonomic classification table
merged_taxa <- data.frame(cbind(row.names(merged_taxa),merged_taxa)) 
names(merged_taxa)[1] <- "seq" 
row.names(merged_taxa) <- NULL 
myin <- merge(merged_nochim,merged_taxa,by="seq") 
write.table(myin, paste0(basedir, "/06_SILVA_tables/", filepattern, "_Taxa_merged_count.txt"), quote=F, sep="\t", row.names=F)
}