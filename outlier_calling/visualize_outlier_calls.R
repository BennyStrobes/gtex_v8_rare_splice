args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)




histogram_showing_number_of_multi_tissue_outlier_calls <- function(number_of_outliers_file, output_file) {
	number_of_outliers_data <- read.table(number_of_outliers_file, header=TRUE)

	df <- data.frame(individual_id=factor(number_of_outliers_data$individual, levels=as.character(number_of_outliers_data$individual)), number_of_outliers=number_of_outliers_data$number_of_multi_tissue_outliers)


	bar_plot <- ggplot(data=df, aes(x=individual_id,y=number_of_outliers)) + geom_bar(stat="identity",fill="mediumpurple") +
				labs(x = "Individual", y = "Number of multi-tissue outliers", title="median(p-value) <= .0456") +
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
				theme(plot.title = element_text(hjust = 0.5),text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=9)) 


	ggsave(bar_plot, file=output_file,width = 15,height=10.5,units="cm")

}


bar_plot_showing_number_of_multi_tissue_outliers_at_various_pvalue_thresholds <- function(cross_tissue_outliers, output_file) {
	pvalue_thresholds <- c(.01, .001, .0001, .00001, .000001, .0000001)
	number_outliers <- c()
	for (iter in 1:length(pvalue_thresholds)) {
		pvalue_threshold <- pvalue_thresholds[iter]
		number_outliers_at_threshold <-sum(cross_tissue_outliers < pvalue_threshold, na.rm=TRUE)
		number_outliers <- c(number_outliers, number_outliers_at_threshold)
	}

	df <- data.frame(pvalue_threshold=factor(pvalue_thresholds, levels=as.character(pvalue_thresholds)), number_of_outliers=number_outliers)


	bar_plot <- ggplot(data=df, aes(x=pvalue_threshold,y=number_of_outliers)) + geom_bar(stat="identity",fill="steelblue3") + 
				labs(x = "Median(p-value) threshold", y = "Number of multi-tissue outliers") +
				theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) 


	ggsave(bar_plot, file=output_file,width = 15,height=10.5,units="cm")

}


#######################
# Command line args
#######################

tissue_names_file = args[1]
covariate_method = args[2]
splicing_outlier_dir = args[3]
splicing_outlier_visualization_dir = args[4]

# Extract vector tissue names
tissue_names <- as.character(unlist(read.table(tissue_names_file,header=FALSE), use.names=FALSE))


cross_tissue_outlier_file <- paste0(splicing_outlier_dir, "cross_tissue_covariate_method_", covariate_method, "_no_global_outliers_ea_only_emperical_pvalue.txt")
cross_tissue_outliers <- read.table(cross_tissue_outlier_file, header=TRUE)
clusters <- cross_tissue_outliers[,1]
cross_tissue_outliers <- cross_tissue_outliers[,2:dim(cross_tissue_outliers)[2]]

# Make bar plot of number outliers at various pvalue thresholds
output_file <- paste0(splicing_outlier_visualization_dir, "number_of_multi_tissue_outliers_ea_only_at_various_pvalue_thresholds_bar_plot.pdf")
bar_plot_showing_number_of_multi_tissue_outliers_at_various_pvalue_thresholds(cross_tissue_outliers, output_file)



# Plot histogram of number of multi-tissue outliers per individual
number_of_outliers_file <- paste0(splicing_outlier_dir, "number_of_multi_tissue_outliers_ea_only.txt")  # Input file
output_file <- paste0(splicing_outlier_visualization_dir, "number_of_multi_tissue_outliers_ea_only_histogram.pdf")
histogram_showing_number_of_multi_tissue_outlier_calls(number_of_outliers_file, output_file)

# Plot histogram of number of multi-tissue outliers per individual (not just ea only)
number_of_outliers_file <- paste0(splicing_outlier_dir, "number_of_multi_tissue_outliers.txt")  # Input file
output_file <- paste0(splicing_outlier_visualization_dir, "number_of_multi_tissue_outliers_histogram.pdf")
histogram_showing_number_of_multi_tissue_outlier_calls(number_of_outliers_file, output_file)
