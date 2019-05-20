args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)



make_num_variants_boxplot <- function(num_variants_file, output_file) {
	options(bitmapType = 'cairo', device = 'pdf')

	num_variants_df <- read.table(num_variants_file, header=TRUE)
	# Re-format data into reasonable setting
	novel <- num_variants_df$novel_count
	singleton <- num_variants_df$singleton_count
	doubleton <- num_variants_df$doubleton_count
	bin1 <- num_variants_df[,6] 
	bin2 <- num_variants_df[,7]

	num_variants <- c(novel, singleton, doubleton, bin1, bin2)
	labels <- c(rep("novel", length(novel)), rep("gnomad singleton", length(singleton)), rep("gnomad doubleton", length(doubleton)), rep("doubleton < MAF <= .001", length(bin1)), rep(".001 < MAF <= .01", length(bin2)))



	df <- data.frame(num_variants=num_variants, label=factor(labels, levels=c(".001 < MAF <= .01", "doubleton < MAF <= .001", "gnomad doubleton", "gnomad singleton", "novel")))

	plotter <- ggplot(df, aes(x=label, y=num_variants, fill=label)) + geom_violin() + 
		theme(text = element_text(size=12),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=12), legend.title = element_text(size=11)) +
		labs(x="", y="Number of variants per individaul", fill="") + 
		#scale_fill_manual(values=c("dodgerblue", "violetred1")) +
		theme(legend.position="none") + 
		ylim(0,50000) +
		coord_flip()
	ggsave(plotter, file=output_file,width = 28,height=8,units="cm")

}








num_variants_file <- args[1]
visualization_dir <- args[2]


output_file <- paste0(visualization_dir, "num_variants_boxplot.pdf")
make_num_variants_boxplot(num_variants_file, output_file)
