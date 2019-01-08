args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)


cross_tissue_branchpoint_enrichment_errorbar_plot <- function(cross_tissue_enrichment_file, output_file) {
	enrichments <- read.table(cross_tissue_enrichment_file, header=FALSE)
	num_pvalues <- dim(enrichments)[1]
	pvalue_string <- c()
	# Initialize vectors
	pvalues <- c()
	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()
	# Loop through pvalues
	for (pvalue_iter in 1:num_pvalues) {
		pvalue <- as.character(as.numeric(strsplit(as.character(enrichments[pvalue_iter,1]),"_")[[1]][3]))
		pvalue_string <- c(pvalue_string, pvalue)
		# Compute odds ratios for this line
		a <- enrichments[pvalue_iter, 2]
		b <- enrichments[pvalue_iter, 3]
		c <- enrichments[pvalue_iter, 4]
		d <- enrichments[pvalue_iter, 5]
		orat <- (a/b)/(c/d)
		# Compute error bars for this orat
		log_orat <- log(orat)
		log_bounds <- 1.96*sqrt((1.0/a) - (1.0/b) + (1.0/c) - (1.0/d))
		upper_bound <- orat*exp(log_bounds)
		lower_bound <- orat*exp(-log_bounds)


		# Add data to vectors
		pvalues <- c(pvalues, pvalue)
		odds_ratios <- c(odds_ratios, orat)
		lower_bounds <- c(lower_bounds, lower_bound)
		upper_bounds <- c(upper_bounds, upper_bound)
	}
	# Make data frame
	df <- data.frame(pvalues=factor(pvalues,levels=(as.character(pvalue_string))), odds_ratios=odds_ratios, lower_bounds=lower_bounds, upper_bounds=upper_bounds)
	dodge <- position_dodge(width=0.9)
	error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=pvalues,ymin=lower_bounds, ymax=upper_bounds), position=dodge) +
					geom_point(data=df, mapping=aes(x=pvalues, y=odds_ratios), position=dodge) +
					labs(x = "p-value", y = "Enrichment") +
					geom_hline(yintercept=1.0) + 
					theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) 

	ggsave(error_bar_plot, file=output_file,width = 19,height=10.5,units="cm")

}



branchpoint_enrichment_dir <- args[1]
visualize_branchpoint_enrichment_dir <- args[2]


distance_window <- "5"
cross_tissue_enrichment_file <- paste0(branchpoint_enrichment_dir, "cross_tissue_branchpoint_", distance_window, "_enrichment.txt")
output_file <- paste0(visualize_branchpoint_enrichment_dir, "cross_tissue_branchpoint_", distance_window, "_enrichment_errorbar.pdf")

cross_tissue_branchpoint_enrichment_errorbar_plot(cross_tissue_enrichment_file, output_file)