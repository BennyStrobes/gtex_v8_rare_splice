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




make_gtex_posterior_vs_amish_pvalue_scatter <- function(merged_data) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_watershed_posterior

	corry <- cor.test(posterior, log_pval, method="spearman")
	pvalue = corry[3]
	correlation = corry[4]

	df <- data.frame(pvalue=log_pval, posterior=posterior)
	scatter <- ggplot(df, aes(x=pvalue, y=posterior)) +
 			 geom_point() +
 			 gtex_v8_figure_theme() + 
 			 labs(x="median -log10(pvalue) in Amish cohort",y="median GTEx Watershed splicing posterior", title=paste0("Spearman rho: ", correlation, " / Spearman pvalue: ", pvalue))
 	return(scatter)

}

make_gtex_gam_posterior_vs_amish_pvalue_scatter <- function(merged_data) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_gam_posterior

	corry <- cor.test(posterior, log_pval, method="spearman")
	pvalue = corry[3]
	correlation = corry[4]

	df <- data.frame(pvalue=log_pval, posterior=posterior)
	scatter <- ggplot(df, aes(x=pvalue, y=posterior)) +
 			 geom_point() +
 			 gtex_v8_figure_theme() + 
 			 labs(x="median -log10(pvalue) in Amish cohort",y="median GTEx GAM splicing posterior", title=paste0("Spearman rho: ", correlation, " / Spearman pvalue: ", pvalue))
 	return(scatter)

}

make_gtex_posterior_vs_amish_pvalue_boxplot <- function(merged_data) {
	log_pval = -log10(merged_data$amish_splicing_pvalue + 1e-6)
	posterior = merged_data$average_watershed_splice_posterior

	high_log_pvalue_indices = log_pval > 3
	low_log_pvalue_indices = log_pval < 1

	high_watershed_pvalz = posterior[high_log_pvalue_indices]
	low_watershed_pvalz = posterior[low_log_pvalue_indices]

	posterior <- c(high_watershed_pvalz, low_watershed_pvalz)
	type <- c(rep("outlier", length(high_watershed_pvalz)), rep("inlier", length(low_watershed_pvalz)))

	df <- data.frame(posterior=posterior, type=factor(type, levels=c("inlier", "outlier")))

	p <- ggplot(df, aes(x=type, y=log(posterior))) + 
  		geom_boxplot() +
  		gtex_v8_figure_theme() + 
  		labs(x="-log10(pvalue) of splicing outlier in Amish cohort", y="log(Watershed Posterior in GTEx)")
  	return(p)

}






merged_data_set_file <- args[1]
output_dir <- args[2]



merged_data <- read.table(merged_data_set_file,header=TRUE)

###########################
# make scatterplot of watershed posteriors vs amish splicing pvalue
###########################
output_file <- paste0(output_dir, "gtex_posterior_vs_amish_pvalue_scatter.pdf")
scatter <- make_gtex_posterior_vs_amish_pvalue_scatter(merged_data)
ggsave(scatter, file=output_file, width=7.2, height=4,units="in")

output_file <- paste0(output_dir, "gtex_gam_posterior_vs_amish_pvalue_scatter.pdf")
scatter <- make_gtex_gam_posterior_vs_amish_pvalue_scatter(merged_data)
ggsave(scatter, file=output_file, width=7.2, height=4,units="in")


###########################
# make boxplot of amish splicing pvalues at different watershed thresholds
###########################
output_file <- paste0(output_dir, "gtex_posterior_vs_amish_pvalue_boxplot.pdf")
#boxplot <- make_gtex_posterior_vs_amish_pvalue_boxplot(merged_data)
#ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")


