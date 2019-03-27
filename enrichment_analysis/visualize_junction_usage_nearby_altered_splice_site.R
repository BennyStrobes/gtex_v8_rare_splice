args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)




odds_ratio_boxplot_across_tissues <- function(aa, output_file) {
	options(bitmapType = 'cairo', device = 'pdf')
	df_to_concensus <- aa[as.character(aa$to_or_from_concensus) == "to_concensus",]
	df_from_concensus <- aa[as.character(aa$to_or_from_concensus) == "from_concensus",]
	orat <- c()
	type <- c()
	orat <- c(orat, df_to_concensus$odds_ratio, df_from_concensus$odds_ratio)
	type <- c(type, rep("Consensus Created", length(df_to_concensus$odds_ratio)), rep("Consensus Destroyed", length(df_from_concensus$odds_ratio)))
	df <- data.frame(odds_ratio=log(orat), type=as.factor(type))
	plotter <- ggplot(df, aes(x=type, y=odds_ratio, fill=type)) + geom_boxplot() + 
		theme(text = element_text(size=12),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=12), legend.title = element_text(size=11)) +
		labs(x="", y="ln(Odds Ratio)", fill="") + 
		scale_fill_manual(values=c("dodgerblue", "violetred1")) +
		theme(legend.position="none") + 
		geom_hline(yintercept=0)
	ggsave(plotter, file=output_file,width = 16,height=12,units="cm")

}



jxn_usage_nearby_altered_ss_enrichment_dir = args[1]
visualize_jxn_usage_nearby_altered_ss_enrichment_dir = args[2]


# Load in data
input_file <- paste0(jxn_usage_nearby_altered_ss_enrichment_dir, "tissue_by_tissue_outliers_with_rv_in_concensus_sites_outlier_individuals_1e-06_inlier_individuals_0.5_with_read_counts.txt")
enrichment_data <- read.table(input_file, header=TRUE)

# Visualize enrichment distribution (data aggregrated) tissues
output_file <- paste0(visualize_jxn_usage_nearby_altered_ss_enrichment_dir, "jxn_usage_enrichment_across_tissues.pdf")
odds_ratio_boxplot_across_tissues(enrichment_data, output_file)