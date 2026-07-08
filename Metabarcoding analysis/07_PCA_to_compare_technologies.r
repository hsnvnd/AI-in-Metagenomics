#This script generates a PCA biplot to compare genus composition across samples from different sequencing technologies
# 1) A genus abundance filter is applied to retain genera with a mean relative abundance of at least 0.1% across all samples after merging data from all technologies
# 2) Squared root transformation is applied to genus percentages (Hellinger transformation)
# 3) PCA is performed

#Load libraries
library(data.table)
library(stringr)
library(vegan)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(svglite)

#Set and create the directory where all the results will be stored
basedir <- "path/to/the/result/directory"
dir.create(basedir)

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


#Loop through technologies
tot_taxa <- list()

for (tech in technologies) {

#Load the table 
starting_data <- fread(tech[[1]], header=T, data.table=T)

#Define sample names
sample_names <- as.character(seq(1:36))

#Reorder the columns from 1 to 36 samples
setcolorder(starting_data, c("Genus", sample_names))

#Add technology name to sample names
colnames(starting_data)[2:37] <- paste0(colnames(starting_data)[2:37], "_", tech[[2]])

#Save tech data tables
tot_taxa <- c(tot_taxa, list(starting_data))}

#Merge tables from all technologies
all_tech <- Reduce(function (...) { merge(..., by="Genus", all=TRUE) }, tot_taxa)

#Replace NA with zeros
all_tech[is.na(all_tech)] <- 0

#Get a matrix with only numeric values and with samples as rows and genera as columns
all_tech <- as.data.frame(all_tech)
rownames(all_tech) <- all_tech$Genus
all_tech$Genus <- NULL
all_tech_tran <- t(all_tech)

#Apply genus abundance filter
all_filtered_data_set <- low.count.removal(all_tech_tran, percent=0.1)
filtered_data_set <- all_filtered_data_set$data.filter

#Apply the squared root transformation
squared_root_filtered_data_set <- sqrt(filtered_data_set)


#########################################################
#Use Singular Value Decompositon method to compute a PCA#
#########################################################

#scale option is set to FALSE
pca <- prcomp(squared_root_filtered_data_set, retx=TRUE, center=TRUE, scale.=FALSE)


#######################################
#Statistics for comparing technologies#
#######################################

#Create a dataframe to store PERMANOVA and PERMDISP pvalues
stat_results <- data.frame(matrix(NA, ncol=2, nrow=1))
colnames(stat_results) <- c("PERMDISP", "PERMANOVA")
rownames(stat_results) <- "pvalue"

#Define the main data table
main <- as.data.frame(squared_root_filtered_data_set)
main$Sample_Tech <- rownames(main)
main[c("Sample", "Technology")] <- str_split_fixed(main$Sample_Tech, "_", 2)
main <- main[, c("Sample_Tech", "Sample", "Technology", colnames(main)[1:(length(colnames(main))-3)])]
main$Sample <- as.numeric(main$Sample)
main <- main[order(main$Technology, main$Sample), ]
main$Technology <- as.factor(main$Technology)
main$Sample <- as.factor(main$Sample)

#Define the distance matrix with Euclidean distance
distance_matrix_stat <- vegdist(main[, -(1:3)], method="euclidean")


#PERMDISP
#Perform PERMDISP statistical test to verify if there are significant differences in patterns of dispersion among groups (variability among sample units).

#A global test is performed to test whether at least one technology (group) has a different dispersion.
#To account for paired (dependent) samples, a permutation scheme is defined in order to restrict permutations only within samples.
set.seed(123)
perm_scheme_dispersion <- how(within=Within(type="free"), blocks=main$Sample, nperm=999)

#Compute distances from centroids
set.seed(123)
bd_tech <- betadisper(d=distance_matrix_stat, group=main$Technology)

#Perform permutational ANOVA
set.seed(123)
global_ANOVA_dispersion <- permutest(bd_tech, permutations=perm_scheme_dispersion)

#Save PERMDISP - permutational ANOVA pvalue
stat_results[1, 1] <- global_ANOVA_dispersion$tab$`Pr(>F)`[1]

#Compare technology dispersions in a pairwise manner
pairwise_dispersion_results <- data.frame(matrix(NA, ncol=3, nrow=1))
colnames(pairwise_dispersion_results) <- c("MiSeq-AVITI", "MiSeq-NovaSeq", "AVITI-NovaSeq")
rownames(pairwise_dispersion_results) <- "pvalue"

pairwise_dispersion <- list(z1=list("MiSeq-AVITI", "NovaSeq"), z2=list("MiSeq-NovaSeq", "AVITI"), z3=list("AVITI-NovaSeq", "MiSeq"))

for (disp in pairwise_dispersion){

#Filter main table
pair_disp_main <- main[main$Technology!=disp[[2]],]
pair_disp_main$Technology <- droplevels(pair_disp_main$Technology)

#Define the distance matrix with Euclidean distance between the selected technologies
pair_disp_distance_matrix <- vegdist(pair_disp_main[, -(1:3)], method="euclidean")

#Define permutation scheme
set.seed(123)
pair_disp_perm_scheme <- how(within=Within(type="free"), blocks=pair_disp_main$Sample, nperm=999)

#Compute distances from centroids
set.seed(123)
bd_pair_tech <- betadisper(d=pair_disp_distance_matrix, group=pair_disp_main$Technology)

#Perform permutational ANOVA
set.seed(123)
pair_ANOVA_dispersion <- permutest(bd_pair_tech, permutations=pair_disp_perm_scheme)

#Save pairwise pvalues
pairwise_dispersion_results[1, disp[[1]]] <- pair_ANOVA_dispersion$tab$`Pr(>F)`[1]
}
tran_pairwise_dispersion_results <- as.data.frame(t(pairwise_dispersion_results))

#Adjust pvalues for multiple comparisons to control FDR
tran_pairwise_dispersion_results$adjusted_pvalue <- p.adjust(tran_pairwise_dispersion_results$pvalue, method="BH")


#PERMANOVA
#Permutations are restricted within samples. Technology labels are reassigned only within the same sample in order to respect the experimental design: the same 36 samples are analyzed with different technologies.
#Define permutation scheme
set.seed(123)
perm_scheme <- how(within=Within(type="free"), blocks=main$Sample, nperm=999)

#Compute PERMANOVA
set.seed(123)
permanova <- adonis2(distance_matrix_stat~Technology, data=main, permutations=perm_scheme)

#Save PERMANOVA pvalue
stat_results[1, 2] <- permanova$`Pr(>F)`[1]

#Test manually for pairwise differences between PERMANOVA groups
pairwise_permanova_results <- data.frame(matrix(NA, ncol=3, nrow=1))
colnames(pairwise_permanova_results) <- c("MiSeq-AVITI", "MiSeq-NovaSeq", "AVITI-NovaSeq")
rownames(pairwise_permanova_results) <- "pvalue"

pairwise <- list(z4=list("MiSeq-AVITI", "NovaSeq"), z5=list("MiSeq-NovaSeq", "AVITI"), z6=list("AVITI-NovaSeq", "MiSeq"))

for (p_w in pairwise){

#Filter main table
pair_main <- main[main$Technology!=p_w[[2]],]
pair_main$Technology <- droplevels(pair_main$Technology)

#Define the distance matrix with Euclidean distance between the selected technologies
pair_distance_matrix <- vegdist(pair_main[, -(1:3)], method="euclidean")

#Define permutation scheme
set.seed(123)
pair_perm_scheme <- how(within=Within(type="free"), blocks=pair_main$Sample, nperm=999)

#Compute PERMANOVA
set.seed(123)
pair_permanova <- adonis2(pair_distance_matrix~Technology, data=pair_main, permutations=pair_perm_scheme)

#Save pairwise pvalues
pairwise_permanova_results[1, p_w[[1]]] <- pair_permanova$`Pr(>F)`[1]
}
tran_pairwise_permanova_results <- as.data.frame(t(pairwise_permanova_results))

#Adjust pvalues for multiple comparisons to control FDR
tran_pairwise_permanova_results$adjusted_pvalue <- p.adjust(tran_pairwise_permanova_results$pvalue, method="BH")


#Save statistical outputs
#Save the dataframe with PERMANOVA and PERMDISP pvalues
write.table(stat_results, paste0(basedir, "/compared_tech_PERMDISP_PERMANOVA_pvalues.txt"), sep="\t", quote=F, col.names=T, row.names=T)

#Save the PERMDISP pairwise comparisons
write.table(tran_pairwise_dispersion_results, paste0(basedir, "/compared_tech_pairwise_PERMDISP_dispersion_pvalues.txt"), sep="\t", quote=F, col.names=T, row.names=T)

#Save the PERMANOVA pairwise comparisons
write.table(tran_pairwise_permanova_results, paste0(basedir, "/compared_tech_pairwise_PERMANOVA_pvalues.txt"), sep="\t", quote=F, col.names=T, row.names=T)

#Save the global PERMANOVA output 
capture.output(permanova, file=paste0(basedir, "/compared_tech_global_PERMANOVA_output.txt"))


############
#Scree plot#
############

scree_plot_df <- data.frame(Num_component=seq(1:length(pca$sdev)), Components=sort((pca$sdev^2)*100/sum(pca$sdev^2), decreasing=TRUE))

scree_plot <- ggplot(scree_plot_df, aes(x=Num_component, y=Components)) +
	geom_col() +
	labs(title=paste0("Scree plot"), x="Components", y="Variance explained (%)") +
	theme_classic() +
	theme(plot.title=element_text(hjust=0.5))
  ggsave(filename=paste0(basedir, "/Scree_plot.png"), scree_plot, width=15)

#Save scree plot as .svg
svglite(paste0(basedir, "/Scree_plot.svg"), width=15)
print(scree_plot)
dev.off()


########
#Biplot#
########

#Use autoplot to generate a biplot with sample coordinates scaled based on the genus arrows
original_plot <- autoplot(pca, scale=1, loadings=TRUE)

#Extract sample coordinates
PC_coord <- original_plot$data[, c("PC1", "PC2")]
PC_coord$Sample_Tech <- rownames(PC_coord)
rownames(PC_coord) <- NULL
PC_coord[c("Sample", "Technology")] <- str_split_fixed(PC_coord$Sample_Tech, "_", 2)
PC_coord$Sample_Tech <- NULL

#Extract genus loadings for Principal Component 1 and 2
variables_loadings <- original_plot$layers[[2]]$data[, c("PC1", "PC2")]
variables_loadings$Genus <- rownames(variables_loadings)
rownames(variables_loadings) <- NULL

#Compute the module of genus vectors on PC1 and PC2
variables_loadings$Vector_module <- sqrt((variables_loadings$PC1)^2 + (variables_loadings$PC2)^2)

#Order genera in descending order based on the vector module
variables_loadings <- variables_loadings[order(-variables_loadings$Vector_module), ]

#Save loadings
write.table(variables_loadings, paste0(basedir, "/Genus_Loadings.txt"), sep="\t", quote=F, col.names=T, row.names=F)

#Choose the top genera with the highest vector module
top_Genus <- variables_loadings[1:20,]

#Add soil information to the dataframe
PC_coord$Soil <- metadata[match(PC_coord$Sample, metadata$Sample_name), Soil]

#Plot
pca_plot <- ggplot(PC_coord, aes(x=PC1, y=PC2, color=factor(Technology, c("MiSeq", "AVITI", "NovaSeq")), label=Sample, shape=Soil)) +
	geom_point() +
	geom_text_repel(size=2.3, min.segment.length=Inf, box.padding=0.17, show.legend=FALSE) +
	labs(title=NULL, x=paste0("PC1 ", round(scree_plot_df[1, 2], 2), "%"), y=paste0("PC2 ", round(scree_plot_df[2, 2], 2), "%"), color="Technology") +
	theme_linedraw() +
	theme(panel.grid=element_blank()) +
	guides(color=guide_legend(order=1), shape=guide_legend(order=2))
	# geom_segment(data=top_Genus, aes(x=0, y=0, xend=PC1, yend=PC2), arrow=arrow(length=unit(0.2, "cm")), color="black", alpha=0.7, inherit.aes=FALSE) +
	# geom_text_repel(data=top_Genus, aes(x=PC1, y=PC2, label=Genus), inherit.aes=FALSE, size=3, box.padding=0.2, min.segment.length=1)
   ggsave(filename=paste0(basedir, "/PCA.png"), pca_plot, width=15)

#Save plot as .svg
svglite(paste0(basedir, "/PCA.svg"), width=15)
print(pca_plot)
dev.off()