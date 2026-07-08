#This script plots the number of reads generated for each sample
#This script computes different alpha diversity metrics using DESeq2-normalized genus count data
#Before computing the alpha diversity metrics, the dataset is filtered to retain genera with a mean relative abundance of at least 0.1% across all samples 
#Relative abundances used for filtering are computed with respect to the total number of reads given to the classifier

#Load libraries
library(data.table)
library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(rstatix)
library(svglite)
library(stringr)

#Set and create the directory where all the results will be stored
basedir <- "path/to/the/result/directory"
dir.create(basedir)

#Empty list to store plots
chart_list <- list()


#################
#Number of reads#
#################

#Set paths
tech_number <- list(a=list("/path/to/the/table/containing/the/number/of/reads/per/sample/for/MiSeq.txt", "MiSeq"),
b=list("/path/to/the/table/containing/the/number/of/reads/per/sample/for/AVITI.txt", "AVITI"),
d=list("/path/to/the/table/containing/the/number/of/reads/per/sample/for/NovaSeq.txt", "NovaSeq"))

#Create a table with read numbers for all technologies
number_final <- rbindlist(lapply(tech_number, function(r){	
	read_number <- fread(r[[1]], data.table=T, header=F)
	read_number[, V1 := str_split_i(V1, "-", 2)]
	read_number <- unique(read_number)
	setnames(read_number, c("Sample", "Reads"))
	read_number <- read_number[match(as.character(seq(1, 36)), Sample)]
	read_number[, Tech := r[[2]]]}), use.names=TRUE)

#Save relevant statistics
final_range <- number_final %>% group_by(Tech) %>% summarise(min_reads=min(Reads), max_reads=max(Reads), range_reads=max(Reads) - min(Reads), mean_reads=mean(Reads), median_reads=median(Reads))
write.table(final_range, paste0(basedir, "/Range_of_reads.txt"), sep="\t", quote=F, col.names=T, row.names=F)


###############################################
#Statistics for the number of reads per sample#
###############################################

#Use a Friedman test to compare the number of reads per sample between different technologies and test if at least one technology has a different number of reads per sample compared to the others.
#Then perform a Wilcoxon signed-rank test between each pair of technologies to identify which one differs.

#Set Sample column as numeric
number_final[, Sample := as.numeric(Sample)]

#Set Tech column as factor
number_final[, Tech := factor(Tech, level=c("MiSeq", "AVITI", "NovaSeq"))]

#Order rows by Sample and then by Tech. The order of the Tech column will follow the levels previously defined.
setorder(number_final, "Sample", "Tech")

#Compute Friedman test
number_friedman <- friedman.test(y=number_final$Reads, groups=number_final$Tech, blocks=number_final$Sample)

#Save Friedman test output
capture.output(number_friedman, file=paste0(basedir, "/stats_package_Friedman_test_number_of_reads.txt"))

#Compute a Wilcoxon signed-rank test for each pair of technologies
number_wilcoxon <- pairwise.wilcox.test(number_final$Reads, number_final$Tech, paired=TRUE, p.adjust.method="BH", alternative="two.sided", exact=FALSE)

#Save Wilcoxon signed-rank test output
capture.output(number_wilcoxon, file=paste0(basedir, "/stats_package_Wilcoxon_signed-rank_test_number_of_reads.txt"))

#Perform again a Wilcoxon signed-rank test to compare the number of reads between each pair of technologies
#This time dyplr and rstatix are used to add pvalues on the charts
number_tests <- number_final %>% 
	wilcox_test(Reads ~ Tech, paired=TRUE, exact=FALSE) %>%   
	adjust_pvalue(method="BH") %>%
	add_significance("p.adj")

#Compute pvalue coordinates
number_tests <- number_tests %>% add_xy_position(x="Tech", fun="max", step.increase=0.04)
capture.output(number_tests, file=paste0(basedir, "/rstatix_package_Wilcoxon_signed-rank_test_number_of_reads.txt"))


##########################
#Plot the number of reads#
##########################

#Plot the number of reads
plot_number <- ggplot(number_final, aes(x=Tech, y=Reads, col=Tech)) +
	geom_boxplot(outlier.shape=NA, linewidth=0.8) +
	geom_point(position=position_jitter(width=0.2), alpha=0.5) +
	labs(title="Number of reads", x=NULL, y=NULL, col="Technologies ") +
	theme_classic() +
	theme(plot.title=element_text(hjust=0.5))
  ggsave(filename=paste0(basedir, "/Number_of_reads.png"), plot_number, width=7, height=7)

#Save the chart as .svg
svglite(paste0(basedir, "/Number_of_reads.svg"), width=7, height=7)
print(plot_number)
dev.off()

#Save the read plot in the list
chart_list <- c(chart_list, list(plot_number))


#########################
#Alpha diversity metrics#
#########################

#Define the abundance filter function
low.count.removal <- function(
                        data, #Genus count dataframe of size n (samples) x p (genera)
                        inputs, #Vector with the total number of reads for each sample before classification
						percent=0.1 #Cutoff chosen
						) 
  {
    if (!all(rownames(data) %in% names(inputs))) {stop("Some rownames in 'data' are not present in 'inputs'")}	
	rs <- inputs[rownames(data)]
	rel_abundances <- sweep(data, 1, rs, "/") * 100
	mean_abundances <- colMeans(rel_abundances)
	keep.genus <- mean_abundances >= percent
    data.filter <- data[,keep.genus, drop=FALSE]
    return(list(data.filter=data.filter, keep.genus=keep.genus))
}

#Alpha metrics to loop on 
alpha_metrics <- c("Observed", "Shannon", "Simpson")

#Vector to store all warnings
all_warnings <- c()

#Variables to loop on for alpha diversity
technologies <- list(a=list("/path/to/DESeq2-normalized/genus/count/table_Miseq.txt", "MiSeq", "path/to/the/table/containing/the/per-sample/number/of/sequences/given/to/the/classifier_Miseq.txt"), 
b=list("/path/to/DESeq2-normalized/genus/count/table_Aviti.txt", "AVITI", "path/to/the/table/containing/the/per-sample/number/of/sequences/given/to/the/classifier_Aviti.txt"),
d=list("/path/to/DESeq2-normalized/genus/count/table_Novaseq.txt", "NovaSeq", "path/to/the/table/containing/the/per-sample/number/of/sequences/given/to/the/classifier_Novaseq.txt"))

#Loop through the alpha metrics
for (metric in alpha_metrics){

all_technologies <- list()

for (tech in technologies){

#Load the table
starting_data <- fread(tech[[1]], data.table=T, header=T)

#Define sample names
sample_names <- as.character(seq(1:36))

#Reorder the columns from 1 to 36 samples
setcolorder(starting_data, c("Genus", sample_names))

#Retrieve the number of input reads for each sample, give them the sample name and save them
input <- fread(tech[[3]], header=T, data.table=T)
input_numbers <- setNames(input$nochim, input$Sample)

#Get a matrix with only numeric values and with samples as rows and genera as columns
starting_data_2 <- as.data.frame(starting_data)
rownames(starting_data_2) <- starting_data_2$Genus
starting_data_2$Genus <- NULL
starting_data_2_tran <- t(starting_data_2)

#Apply the genus abundance filter
all_filtered_data_set <- low.count.removal(starting_data_2_tran, input_numbers, percent=0.1)
filtered_data_set <- all_filtered_data_set$data.filter

#A matrix with genera as rows and samples as columns is required by phyloseq
filtered_data_set_2 <- t(filtered_data_set)

#Rounded numbers are needed for Observed metric
if (metric=="Observed"){filtered_data_set_2 <- round(filtered_data_set_2)}

#Load the matrix as an otu_table-class object within phyloseq
phyloseq_table <- otu_table(filtered_data_set_2, taxa_are_rows=TRUE)

#Compute the alpha diversity matrix
#Save the warning and the command in the command_warn variable and then in the all_warnings vector
alpha_metric <- withCallingHandlers({
	estimate_richness(phyloseq_table, split=TRUE, measures=metric)},
	warning=function(war) {
	  #Retrieve command and warning and save them. invokeRestart("muffleWarning") prevents the warning from being printed in the console.
      command_warn <- c(deparse(conditionCall(war)), conditionMessage(war))
	  all_warnings <<- c(all_warnings, command_warn)
      invokeRestart("muffleWarning")})	

#Reorder the data frame, retrieve sample data and rename the columns
alpha_metric$Sample <- rownames(alpha_metric)
rownames(alpha_metric) <- NULL
alpha_metric$Sample <- str_remove(alpha_metric$Sample, "X")
alpha_metric <- alpha_metric[, c("Sample", metric)]
names(alpha_metric)[names(alpha_metric)==metric] <- tech[[2]]

#Reshape data from wide to long
alpha_metric_2 <- alpha_metric %>% pivot_longer(cols=-1, names_to=tech[[2]], values_to="Value")

#Rename the Combination column
names(alpha_metric_2)[names(alpha_metric_2)==tech[[2]]] <- "Technology"

#Save the tables in a list
all_technologies <- c(all_technologies, list(alpha_metric_2))
}

#Merge tables from all technologies
alpha_final <- do.call(rbind, all_technologies)


########################################
#Statistics for alpha diversity metrics#
########################################

#Use a Friedman test to compare the relevant metric between different technologies and test if at least one technology has a different metric compared to the others.
#Then perform a Wilcoxon signed-rank test between each pair of technologies to identify which ones differ.

#Convert "alpha_final" in a data table 
alpha_final_dt <- as.data.table(alpha_final)

#Set Sample column as numeric
alpha_final_dt[, Sample := as.numeric(Sample)]

#Set Technology column as factor
alpha_final_dt[, Technology := factor(Technology, level=c("MiSeq", "AVITI", "NovaSeq"))]

#Order rows by Sample and then by Technology. The order of the Technology column will follow the levels previously defined.
setorder(alpha_final_dt, "Sample", "Technology")

#Compute Friedman test
alpha_friedman <- friedman.test(y=alpha_final_dt$Value, groups=alpha_final_dt$Technology, blocks=alpha_final_dt$Sample)

#Save Friedman test output
capture.output(alpha_friedman, file=paste0(basedir, "/stats_package_Friedman_test_", metric, ".txt"))

#Compute a Wilcoxon signed-rank test for each pair of technologies
alpha_wilcoxon <- pairwise.wilcox.test(alpha_final_dt$Value, alpha_final_dt$Technology, paired=TRUE, p.adjust.method="BH", alternative="two.sided", exact=FALSE)

#Save Wilcoxon signed-rank test output
capture.output(alpha_wilcoxon, file=paste0(basedir, "/stats_package_Wilcoxon_signed-rank_test_", metric, ".txt"))

#Perform again a Wilcoxon signed-rank test to compare alpha-diversity metric between each pair of technologies
#This time dyplr and rstatix are used to add pvalues on the charts
tests_alpha <- alpha_final_dt %>% 
	wilcox_test(Value ~ Technology, paired=TRUE, exact=FALSE) %>%   
	adjust_pvalue(method="BH") %>%
	add_significance("p.adj")

#Compute pvalue coordinates
if (metric=="Observed")
{tests_alpha <- tests_alpha %>% add_xy_position(x="Technology", fun="max", step.increase=0.09)} 
else if (metric=="Shannon")
{tests_alpha <- tests_alpha %>% add_xy_position(x="Technology", fun="max", step.increase=0.2)} 
else if (metric=="Simpson") 
{tests_alpha <- tests_alpha %>% add_xy_position(x="Technology", fun="max", step.increase=0.55)}

capture.output(tests_alpha, file=paste0(basedir, "/rstatix_package_Wilcoxon_signed-rank_test_", metric, ".txt"))


##############################
#Plot alpha diversity metrics#
##############################

#Plot the metrics
plot_alpha <- ggplot(alpha_final_dt, aes(x=Technology, y=Value, col=Technology)) + 
	geom_boxplot(outlier.shape=NA, linewidth=0.8) + 
	geom_jitter(alpha=0.5) +
	labs(title=if(metric=="Observed"){paste0(metric, " genera")}else{paste0(metric, " index")}, x=NULL, y=NULL, col="Technologies ") +
	theme_classic() +
	theme(plot.title=element_text(hjust=0.5))
  ggsave(filename=paste0(basedir, "/Alpha_diversity_Genus_", metric, ".png"), plot_alpha, width=7, height=7)
  
#Save the chart as .svg
svglite(paste0(basedir, "/Alpha_diversity_Genus_", metric, ".svg"), width=7, height=7)
print(plot_alpha)
dev.off()

#Save plot in a list
chart_list <- c(chart_list, list(plot_alpha))}


#Plot a final chart
final_chart <- ggarrange(plotlist=chart_list, labels=c("A", "B", "C", "D"), ncol=2, nrow=2, common.legend=TRUE, legend="bottom")

#Save it as .png 
ggsave(filename=paste0(basedir, "/All_plots.png"), final_chart, width=15, height=15)

#Save it as .svg
svglite(paste0(basedir, "/All_plots.svg"), width=15, height=15)
print(final_chart)
dev.off()

#Save all warnings
writeLines(all_warnings, paste0(basedir, "/All_warning_messages.txt"))