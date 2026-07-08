#This script performs PERMANOVA and PERMDISP tests to investigate the effect of soil type, soil autoclaving and root thermal treatment on genus composition for each technology dataset
#For each technology:
# 1) A genus abundance filter is applied to retain genera with a mean relative abundance of at least 0.1% across all samples
# 2) Squared root transformation is applied to genus percentages (Hellinger transformation)
# 3) PERMANOVA and PERMDISP tests are performed using Euclidean distances

#Load libraries
library(data.table)
library(stringr)
library(vegan)

#Set and create the directory where all the results will be stored
basedir <- "path/to/the/result/directory"
dir.create(basedir)
dir.create(paste0(basedir, "/MiSeq"))
dir.create(paste0(basedir, "/AVITI"))
dir.create(paste0(basedir, "/NovaSeq"))

#Define the abundance filter function
low.count.removal <- function(
                        data, #Genus percentage dataframe of size n (samples) x p (genera)
                        percent=0.1 #Cutoff chosen
                        ) 
  {
	mean_abundances <- colMeans(data)
	keep.genus <- mean_abundances >= percent
    data.filter <- data[,keep.genus, drop=FALSE]
    return(list(data.filter=data.filter, keep.genus=keep.genus))
}

#Variables to loop on
technologies <- list(a=list("/path/to/genus/percentage/table_MiSeq.txt", "MiSeq"), 
b=list("/path/to/genus/percentage/table_AVITI.txt", "AVITI"),
d=list("/path/to/genus/percentage/table_NovaSeq.txt", "NovaSeq"))


##########
#Metadata#
##########

#Load metadata
metadata <- fread("/path/to/metadata.txt", header=T, data.table=T)

#Remove the "T1_" from sample names
metadata[, Sample_name := str_remove(Sample_name, "T1_")]

#Set metadata as a data frame
metadata <- as.data.frame(metadata)


########################
#Loop through each path#
########################

for (tech in technologies) {

#Load the table
starting_data <- fread(tech[[1]], header=T, data.table=T)

#Define sample names
sample_names <- as.character(seq(1:36))

#Reorder the columns from 1 to 36 samples
setcolorder(starting_data, c("Genus", sample_names))

#Get a matrix with only numeric values and with samples as rows and Genus as columns
starting_data_2 <- as.data.frame(starting_data)
rownames(starting_data_2) <- starting_data_2$Genus
starting_data_2$Genus <- NULL
starting_data_2_tran <- t(starting_data_2)

#Apply Genus abundance filter
all_filtered_data_set <- low.count.removal(starting_data_2_tran, percent=0.1)
filtered_data_set <- all_filtered_data_set$data.filter

#Apply the squared root transformation
squared_root_filtered_data_set <- sqrt(filtered_data_set)


#####################################################
#Define the main data table for statistical analyses#
#####################################################

main <- as.data.frame(squared_root_filtered_data_set)
main$Sample_name <- as.character(rownames(main))
rownames(main) <- NULL

#Retrieve metadata info. I put "metadata" first to set its columns as the first columns in "main_final". 
main_final <- merge(metadata, main, by="Sample_name", all=TRUE, sort=FALSE)
rownames(main_final) <- main_final$Sample_name


##############################
#PERMANOVA with all variables#
##############################

#Set the dataframe
main_all_var <- main_final
main_all_var$Soil <- factor(main_all_var$Soil, levels=c("Manure", "Peat", "Sand"))
main_all_var$Autoclave <- factor(main_all_var$Autoclave, levels=c("Yes", "No"))
main_all_var$Heat_root <- factor(main_all_var$Heat_root, levels=c("Yes", "No"))

#Define the distance matrix with Euclidean distance
distance_matrix_all_var <- vegdist(main_all_var[, -(1:4)], method="euclidean")

#Compute PERMANOVA using all variables in a formula
#Samples are indipendent, therefore no permutation scheme is defined

#Use the by="margin" parameter
set.seed(123)
permanova_all_var_margin <- adonis2(distance_matrix_all_var ~ Soil + Autoclave + Heat_root + Soil:Autoclave + Soil:Heat_root + Autoclave:Heat_root, data=main_all_var, permutations=999, by="margin")

#Use the by="terms" parameter
set.seed(123)
permanova_all_var_terms <- adonis2(distance_matrix_all_var ~ Soil + Autoclave + Heat_root + Soil:Autoclave + Soil:Heat_root + Autoclave:Heat_root + Soil:Autoclave:Heat_root, data=main_all_var, permutations=999, by="terms")

#Save by="margin" PERMANOVA results
permanova_all_var_margin_df <- as.data.frame(permanova_all_var_margin)
permanova_all_var_margin_df$Element <- rownames(permanova_all_var_margin_df)
permanova_all_var_margin_df <- permanova_all_var_margin_df[, c("Element", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
write.table(permanova_all_var_margin_df, paste0(basedir, "/", tech[[2]], "/All_variables_by_margin_PERMANOVA.txt"), sep="\t", quote=F, col.names=T, row.names=F)

#Save by="terms" PERMANOVA results
permanova_all_var_terms_df <- as.data.frame(permanova_all_var_terms)
permanova_all_var_terms_df$Element <- rownames(permanova_all_var_terms_df)
permanova_all_var_terms_df <- permanova_all_var_terms_df[, c("Element", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
write.table(permanova_all_var_terms_df, paste0(basedir, "/", tech[[2]], "/All_variables_by_terms_PERMANOVA.txt"), sep="\t", quote=F, col.names=T, row.names=F)


#######################################################
#Global test of dispersion for Soil_Autoclave variable#
#######################################################

#Create a new Soil_Autoclave variable
main_all_var$Soil_Autoclave <- paste(main_all_var$Soil, main_all_var$Autoclave, sep="_")
main_all_var$Soil_Autoclave <- factor(main_all_var$Soil_Autoclave, levels=c("Manure_Yes", "Manure_No", "Peat_Yes", "Peat_No", "Sand_Yes", "Sand_No"))
main_all_var <- main_all_var[, c("Sample_name", "Soil", "Autoclave", "Heat_root", "Soil_Autoclave", setdiff(colnames(main_all_var), c("Sample_name", "Soil", "Autoclave", "Heat_root", "Soil_Autoclave")))]

#Samples are indipendent, therefore no permutation scheme is defined.
#Compute distances from centroids
set.seed(123)
betadisper_Soil_Autoclave <- betadisper(d=distance_matrix_all_var, group=main_all_var$Soil_Autoclave)

#Perform permutational ANOVA
set.seed(123)
dispersion_Soil_Autoclave <- permutest(betadisper_Soil_Autoclave, permutations=999)

#Save results
dispersion_Soil_Autoclave_df <- dispersion_Soil_Autoclave$tab
dispersion_Soil_Autoclave_df$Element <- rownames(dispersion_Soil_Autoclave_df)
dispersion_Soil_Autoclave_df <- dispersion_Soil_Autoclave_df[, c("Element", "Df", "Sum Sq", "Mean Sq", "F", "N.Perm", "Pr(>F)")]
write.table(dispersion_Soil_Autoclave_df, paste0(basedir, "/", tech[[2]], "/Soil_Autoclave_variable_global_PERMDSIP.txt"), sep="\t", quote=F, col.names=T, row.names=F)


##########################################################################
#Pairwise comparisons for Soil_Autoclave variable: PERMANOVA and PERMDISP#
##########################################################################

#List to save results
Soil_Autoclave_list <- list()

#Set a counter
counter_Soil_Autoclave <- 1

#Set all possible combinations of 2 elements for Soil_Autoclave variable
combinations_Soil_Autoclave <- combn(levels(main_all_var$Soil_Autoclave), 2, simplify=FALSE)

for (SA in combinations_Soil_Autoclave) {
	
	#Filter the main data set to retain only the two required variables
	main_all_var_Soil_Autoclave <- main_all_var[main_all_var$Soil_Autoclave==SA[1] | main_all_var$Soil_Autoclave==SA[2],]
	main_all_var_Soil_Autoclave <- droplevels(main_all_var_Soil_Autoclave)
	
	#Define the distance matrix with Euclidean distance
	distance_matrix_all_var_Soil_Autoclave <- vegdist(main_all_var_Soil_Autoclave[, -(1:5)], method="euclidean")

	#Compute pairwise PERMANOVA
	#Samples are indipendent, therefore no permutation scheme is defined
	set.seed(123)
	pair_permanova_all_var_Soil_Autoclave <- adonis2(distance_matrix_all_var_Soil_Autoclave ~ Soil_Autoclave, data=main_all_var_Soil_Autoclave, permutations=999)

	#Compute pairwise PERMDISP
	#Samples are indipendent, therefore no permutation scheme is defined
	#Compute distances from centroids
	set.seed(123)
	pair_betadisper_all_var_Soil_Autoclave <- betadisper(d=distance_matrix_all_var_Soil_Autoclave, group=main_all_var_Soil_Autoclave$Soil_Autoclave)

	#Perform permutational ANOVA
	set.seed(123)
	pair_dispersion_all_var_Soil_Autoclave <- permutest(pair_betadisper_all_var_Soil_Autoclave, permutations=999)
	
	#Save results
	Soil_Autoclave_list[[counter_Soil_Autoclave]] <- data.table(Soil_Autoclave_comparison=paste(SA, collapse=" vs "), PERMANOVA_pvalue=as.data.frame(pair_permanova_all_var_Soil_Autoclave)["Model", "Pr(>F)"], PERMDISP_pvalue=pair_dispersion_all_var_Soil_Autoclave$tab["Groups", "Pr(>F)"])
	
	#Update the counter
	counter_Soil_Autoclave	<- counter_Soil_Autoclave + 1}

#Combine all results
final_Soil_Autoclave <- rbindlist(Soil_Autoclave_list, use.names=TRUE)

#Apply BH method for FDR correction
final_Soil_Autoclave$adjusted_PERMANOVA_pvalue <- p.adjust(final_Soil_Autoclave$PERMANOVA_pvalue, method = "BH")
final_Soil_Autoclave$adjusted_PERMDISP_pvalue <- p.adjust(final_Soil_Autoclave$PERMDISP_pvalue, method = "BH")

#Save results
write.table(final_Soil_Autoclave, paste0(basedir, "/", tech[[2]], "/Soil_Autoclave_variable_pairwise_PERMANOVA_PERMDSIP_pvalues.txt"), sep="\t", quote=F, col.names=T, row.names=F)


#######################################################
#Global test of dispersion for Soil_Heat_root variable#
#######################################################

#Create a new Soil_Heat_root variable
main_all_var$Soil_Heat_root <- paste(main_all_var$Soil, main_all_var$Heat_root, sep="_")
main_all_var$Soil_Heat_root <- factor(main_all_var$Soil_Heat_root, levels=c("Manure_Yes", "Manure_No", "Peat_Yes", "Peat_No", "Sand_Yes", "Sand_No"))
main_all_var <- main_all_var[, c("Sample_name", "Soil", "Autoclave", "Heat_root", "Soil_Autoclave", "Soil_Heat_root", setdiff(colnames(main_all_var), c("Sample_name", "Soil", "Autoclave", "Heat_root", "Soil_Autoclave", "Soil_Heat_root")))]

#Samples are indipendent, therefore no permutation scheme is defined.
#Compute distances from centroids
set.seed(123)
betadisper_Soil_Heat_root <- betadisper(d=distance_matrix_all_var, group=main_all_var$Soil_Heat_root)

#Perform permutational ANOVA
set.seed(123)
dispersion_Soil_Heat_root <- permutest(betadisper_Soil_Heat_root, permutations=999)

#Save results
dispersion_Soil_Heat_root_df <- dispersion_Soil_Heat_root$tab
dispersion_Soil_Heat_root_df$Element <- rownames(dispersion_Soil_Heat_root_df)
dispersion_Soil_Heat_root_df <- dispersion_Soil_Heat_root_df[, c("Element", "Df", "Sum Sq", "Mean Sq", "F", "N.Perm", "Pr(>F)")]
write.table(dispersion_Soil_Heat_root_df, paste0(basedir, "/", tech[[2]], "/Soil_Heat_root_variable_global_PERMDSIP.txt"), sep="\t", quote=F, col.names=T, row.names=F)


##########################################################################
#Pairwise comparisons for Soil_Heat_root variable: PERMANOVA and PERMDISP#
##########################################################################

#List to save results
Soil_Heat_root_list <- list()

#Set a counter
counter_Soil_Heat_root <- 1

#Set all possible combinations of 2 elements for Soil_Heat_root variable
combinations_Soil_Heat_root <- combn(levels(main_all_var$Soil_Heat_root), 2, simplify=FALSE)

for (SH in combinations_Soil_Heat_root) {
	
	#Filter the main data set to retain only the two required variables
	main_all_var_Soil_Heat_root <- main_all_var[main_all_var$Soil_Heat_root==SH[1] | main_all_var$Soil_Heat_root==SH[2],]
	main_all_var_Soil_Heat_root <- droplevels(main_all_var_Soil_Heat_root)
	
	#Define the distance matrix with Euclidean distance
	distance_matrix_all_var_Soil_Heat_root <- vegdist(main_all_var_Soil_Heat_root[, -(1:6)], method="euclidean")

	#Compute pairwise PERMANOVA
	#Samples are indipendent, therefore no permutation scheme is defined
	set.seed(123)
	pair_permanova_all_var_Soil_Heat_root <- adonis2(distance_matrix_all_var_Soil_Heat_root ~ Soil_Heat_root, data=main_all_var_Soil_Heat_root, permutations=999)

	#Compute pairwise PERMDISP
	#Samples are indipendent, therefore no permutation scheme is defined
	#Compute distances from centroids
	set.seed(123)
	pair_betadisper_all_var_Soil_Heat_root <- betadisper(d=distance_matrix_all_var_Soil_Heat_root, group=main_all_var_Soil_Heat_root$Soil_Heat_root)

	#Perform permutational ANOVA
	set.seed(123)
	pair_dispersion_all_var_Soil_Heat_root <- permutest(pair_betadisper_all_var_Soil_Heat_root, permutations=999)
	
	#Save results
	Soil_Heat_root_list[[counter_Soil_Heat_root]] <- data.table(Soil_Heat_root_comparison=paste(SH, collapse=" vs "), PERMANOVA_pvalue=as.data.frame(pair_permanova_all_var_Soil_Heat_root)["Model", "Pr(>F)"], PERMDISP_pvalue=pair_dispersion_all_var_Soil_Heat_root$tab["Groups", "Pr(>F)"])
	
	#Update the counter
	counter_Soil_Heat_root	<- counter_Soil_Heat_root + 1}

#Combine all results
final_Soil_Heat_root <- rbindlist(Soil_Heat_root_list, use.names=TRUE)

#Apply BH method for FDR correction
final_Soil_Heat_root$adjusted_PERMANOVA_pvalue <- p.adjust(final_Soil_Heat_root$PERMANOVA_pvalue, method = "BH")
final_Soil_Heat_root$adjusted_PERMDISP_pvalue <- p.adjust(final_Soil_Heat_root$PERMDISP_pvalue, method = "BH")

#Save results
write.table(final_Soil_Heat_root, paste0(basedir, "/", tech[[2]], "/Soil_Heat_root_variable_pairwise_PERMANOVA_PERMDSIP_pvalues.txt"), sep="\t", quote=F, col.names=T, row.names=F)


############################################################
#Global test of dispersion for Autoclave_Heat_root variable#
############################################################

#Create a new Autoclave_Heat_root variable
main_all_var$Autoclave_Heat_root <- paste(main_all_var$Autoclave, main_all_var$Heat_root, sep="_")
main_all_var$Autoclave_Heat_root <- factor(main_all_var$Autoclave_Heat_root, levels=c("Yes_Yes", "Yes_No", "No_Yes", "No_No"))
main_all_var <- main_all_var[, c("Sample_name", "Soil", "Autoclave", "Heat_root", "Soil_Autoclave", "Soil_Heat_root", "Autoclave_Heat_root", setdiff(colnames(main_all_var), c("Sample_name", "Soil", "Autoclave", "Heat_root", "Soil_Autoclave", "Soil_Heat_root", "Autoclave_Heat_root")))]

#Samples are indipendent, therefore no permutation scheme is defined.
#Compute distances from centroids
set.seed(123)
betadisper_Autoclave_Heat_root <- betadisper(d=distance_matrix_all_var, group=main_all_var$Autoclave_Heat_root)

#Perform permutational ANOVA
set.seed(123)
dispersion_Autoclave_Heat_root <- permutest(betadisper_Autoclave_Heat_root, permutations=999)

#Save results
dispersion_Autoclave_Heat_root_df <- dispersion_Autoclave_Heat_root$tab
dispersion_Autoclave_Heat_root_df$Element <- rownames(dispersion_Autoclave_Heat_root_df)
dispersion_Autoclave_Heat_root_df <- dispersion_Autoclave_Heat_root_df[, c("Element", "Df", "Sum Sq", "Mean Sq", "F", "N.Perm", "Pr(>F)")]
write.table(dispersion_Autoclave_Heat_root_df, paste0(basedir, "/", tech[[2]], "/Autoclave_Heat_root_variable_global_PERMDSIP.txt"), sep="\t", quote=F, col.names=T, row.names=F)


###############################################################################
#Pairwise comparisons for Autoclave_Heat_root variable: PERMANOVA and PERMDISP#
###############################################################################

#List to save results
Autoclave_Heat_root_list <- list()

#Set a counter
counter_Autoclave_Heat_root <- 1

#Set all possible combinations of 2 elements for Autoclave_Heat_root variable
combinations_Autoclave_Heat_root <- combn(levels(main_all_var$Autoclave_Heat_root), 2, simplify=FALSE)

for (AH in combinations_Autoclave_Heat_root) {
	
	#Filter the main data set to retain only the two required variables
	main_all_var_Autoclave_Heat_root <- main_all_var[main_all_var$Autoclave_Heat_root==AH[1] | main_all_var$Autoclave_Heat_root==AH[2],]
	main_all_var_Autoclave_Heat_root <- droplevels(main_all_var_Autoclave_Heat_root)
	
	#Define the distance matrix with Euclidean distance
	distance_matrix_all_var_Autoclave_Heat_root <- vegdist(main_all_var_Autoclave_Heat_root[, -(1:7)], method="euclidean")

	#Compute pairwise PERMANOVA
	#Samples are indipendent, therefore no permutation scheme is defined
	set.seed(123)
	pair_permanova_all_var_Autoclave_Heat_root <- adonis2(distance_matrix_all_var_Autoclave_Heat_root ~ Autoclave_Heat_root, data=main_all_var_Autoclave_Heat_root, permutations=999)

	#Compute pairwise PERMDISP
	#Samples are indipendent, therefore no permutation scheme is defined
	#Compute distances from centroids
	set.seed(123)
	pair_betadisper_all_var_Autoclave_Heat_root <- betadisper(d=distance_matrix_all_var_Autoclave_Heat_root, group=main_all_var_Autoclave_Heat_root$Autoclave_Heat_root)

	#Perform permutational ANOVA
	set.seed(123)
	pair_dispersion_all_var_Autoclave_Heat_root <- permutest(pair_betadisper_all_var_Autoclave_Heat_root, permutations=999)
	
	#Save results
	Autoclave_Heat_root_list[[counter_Autoclave_Heat_root]] <- data.table(Autoclave_Heat_root_comparison=paste(AH, collapse=" vs "), PERMANOVA_pvalue=as.data.frame(pair_permanova_all_var_Autoclave_Heat_root)["Model", "Pr(>F)"], PERMDISP_pvalue=pair_dispersion_all_var_Autoclave_Heat_root$tab["Groups", "Pr(>F)"])
	
	#Update the counter
	counter_Autoclave_Heat_root	<- counter_Autoclave_Heat_root + 1}

#Combine all results
final_Autoclave_Heat_root <- rbindlist(Autoclave_Heat_root_list, use.names=TRUE)

#Apply BH method for FDR correction
final_Autoclave_Heat_root$adjusted_PERMANOVA_pvalue <- p.adjust(final_Autoclave_Heat_root$PERMANOVA_pvalue, method = "BH")
final_Autoclave_Heat_root$adjusted_PERMDISP_pvalue <- p.adjust(final_Autoclave_Heat_root$PERMDISP_pvalue, method = "BH")

#Save results
write.table(final_Autoclave_Heat_root, paste0(basedir, "/", tech[[2]], "/Autoclave_Heat_root_variable_pairwise_PERMANOVA_PERMDSIP_pvalues.txt"), sep="\t", quote=F, col.names=T, row.names=F)
}