args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)
options(bitmapType = 'cairo', device = 'pdf')

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}




make_gtex_posterior_vs_amish_pvalue_scatter <- function(merged_data, outlier_type) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_watershed_posterior

	corry <- cor.test(posterior, log_pval, method="spearman")
	pvalue = corry[3]
	correlation = corry[4]

	df <- data.frame(pvalue=log_pval, posterior=posterior)
	scatter <- ggplot(df, aes(x=pvalue, y=posterior)) +
 			 geom_point() +
 			 gtex_v8_figure_theme() + 
 			 labs(x=paste0("median -log10(pvalue) ", outlier_type, " in Amish cohort"),y=paste0("median GTEx Watershed ", outlier_type, " posterior"), title=paste0("Spearman rho: ", correlation, " / Spearman pvalue: ", pvalue))
 	return(scatter)

}

make_gtex_gam_posterior_vs_amish_pvalue_scatter <- function(merged_data, outlier_type) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_gam_posterior

	corry <- cor.test(posterior, log_pval, method="spearman")
	pvalue = corry[3]
	correlation = corry[4]

	df <- data.frame(pvalue=log_pval, posterior=posterior)
	scatter <- ggplot(df, aes(x=pvalue, y=posterior)) +
 			 geom_point() +
 			 gtex_v8_figure_theme() + 
 			 labs(x=paste0("median -log10(pvalue) ", outlier_type, " in Amish cohort"),y=paste0("median GTEx GAM ", outlier_type, " posterior"), title=paste0("Spearman rho: ", correlation, " / Spearman pvalue: ", pvalue))
 	return(scatter)

}

make_gtex_posterior_vs_amish_pvalue_boxplot <- function(merged_data, outlier_type) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_watershed_posterior

	high_watershed_indices = posterior > .8
	# med_watershed_indices = posterior > .4 & posterior <=.8
	low_watershed_indices = posterior < .01

	high_watershed_pvalz = log_pval[high_watershed_indices]
	# med_watershed_pvalz = log_pval[med_watershed_indices]
	low_watershed_pvalz = log_pval[low_watershed_indices]

	# get gam thresholds
	num_high = length(high_watershed_pvalz)
	num_low = length(low_watershed_pvalz)

	high_threshold = sort(merged_data$median_gam_posterior)[(length(merged_data$median_gam_posterior) - num_high)]
	low_threshold = sort(merged_data$median_gam_posterior)[(num_low + 1)]
	high_gam_indices = merged_data$median_gam_posterior > high_threshold
	low_gam_indices = merged_data$median_gam_posterior < low_threshold

	high_gam_pvalz = log_pval[high_gam_indices]
	low_gam_pvalz = log_pval[low_gam_indices]
	print(num_high)
	print(num_low)


	pval <- c(high_watershed_pvalz, low_watershed_pvalz, high_gam_pvalz, low_gam_pvalz)
	type <- c(rep("Median GTEx posterior > .8", length(high_watershed_pvalz)), rep("Median GTEx posterior < .01", length(low_watershed_pvalz)), rep("Median GTEx posterior > .8", length(high_gam_pvalz)), rep("Median GTEx posterior < .01", length(low_gam_pvalz)))
	model <- c(rep("Watershed", length(high_watershed_pvalz)), rep("Watershed", length(low_watershed_pvalz)), rep("GAM", length(high_gam_pvalz)), rep("GAM", length(low_gam_pvalz)))

	print(wilcox.test(high_watershed_pvalz, low_watershed_pvalz))
	print(wilcox.test(high_watershed_pvalz, high_gam_pvalz))
	print(wilcox.test(high_gam_pvalz, low_gam_pvalz))

	df <- data.frame(pval=pval, type=factor(type, levels=c("Median GTEx posterior < .01", "Median GTEx posterior > .8")), model=factor(model, levels=c("Watershed", "GAM")))

	p <- ggplot(df, aes(x=type, y=pval, fill=model)) + 
  		geom_boxplot() +
  		scale_fill_manual(values=c("steelblue3", "firebrick4")) + 
  		gtex_v8_figure_theme() + 
  		labs(title=outlier_type,fill="",x="", y=paste0(outlier_type, " median -log10(pvalue) in Amish cohort"))
  	return(p)

}






merged_splicing_data_set_file <- args[1]
merged_ase_data_set_file <- args[2]
output_dir <- args[3]



splicing_merged_data <- read.table(merged_splicing_data_set_file,header=TRUE)
ase_merged_data <- read.table(merged_ase_data_set_file,header=TRUE)

###########################
# Splicing
############################

# make scatterplot of watershed posteriors vs amish splicing pvalue
###########################
output_file <- paste0(output_dir, "splicing_gtex_posterior_vs_amish_pvalue_scatter.pdf")
scatter <- make_gtex_posterior_vs_amish_pvalue_scatter(splicing_merged_data, "Splicing")
ggsave(scatter, file=output_file, width=7.2, height=4,units="in")

output_file <- paste0(output_dir, "splicing_gtex_gam_posterior_vs_amish_pvalue_scatter.pdf")
scatter <- make_gtex_gam_posterior_vs_amish_pvalue_scatter(splicing_merged_data, "Splicing")
ggsave(scatter, file=output_file, width=7.2, height=4,units="in")

###########################
# make boxplot of amish splicing pvalues at different watershed thresholds
###########################
output_file <- paste0(output_dir, "splicing_gtex_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_posterior_vs_amish_pvalue_boxplot(splicing_merged_data, "Splicing")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")


###########################
# ASE
############################

###########################
# make scatterplot of watershed posteriors vs amish splicing pvalue
###########################
output_file <- paste0(output_dir, "ase_gtex_posterior_vs_amish_pvalue_scatter.pdf")
scatter <- make_gtex_posterior_vs_amish_pvalue_scatter(ase_merged_data, "ASE")
ggsave(scatter, file=output_file, width=7.2, height=4,units="in")

output_file <- paste0(output_dir, "ase_gtex_gam_posterior_vs_amish_pvalue_scatter.pdf")
scatter <- make_gtex_gam_posterior_vs_amish_pvalue_scatter(ase_merged_data, "ASE")
ggsave(scatter, file=output_file, width=7.2, height=4,units="in")
###########################
# make boxplot of amish splicing pvalues at different watershed thresholds
###########################
output_file <- paste0(output_dir, "ase_gtex_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_posterior_vs_amish_pvalue_boxplot(ase_merged_data, "ASE")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")




