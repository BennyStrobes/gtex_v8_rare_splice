args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)
options(warn=1)



histogram_showing_number_of_multi_tissue_outlier_calls <- function(number_of_outliers_file, output_file) {
	number_of_outliers_data <- read.table(number_of_outliers_file, header=TRUE)

	df <- data.frame(individual_id=factor(number_of_outliers_data$individual, levels=as.character(number_of_outliers_data$individual)), number_of_outliers=number_of_outliers_data$number_of_multi_tissue_outliers)


	bar_plot <- ggplot(data=df, aes(x=individual_id,y=number_of_outliers)) + geom_bar(stat="identity",fill="mediumpurple") +
				labs(x = "Individual", y = "Number of multi-tissue outliers", title="median(p-value) <= .01") +
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

corrected_vs_uncorrected_gene_level_pvalue_scatter_plot <- function(pvalue_comparison_file, tissue_type) {
	df <- read.table(pvalue_comparison_file,header=TRUE)
	df$corrected_pvalue <- -log10(df$corrected_pvalue + .0000001)
	df$uncorrected_pvalue <- -log10(df$uncorrected_pvalue + .0000001)

	df <- df[sample(nrow(df), 100000),]

	scatter <- ggplot(df, aes(x=uncorrected_pvalue, y=corrected_pvalue, colour=number_of_clusters)) + geom_point(size=.1) + 
				labs(x = "-log10(uncorrected pvalue)", y = "-log10(corrected pvalue)", colour="Number clusters") +
				gtex_v8_figure_theme()
	return(scatter)
}

corrected_vs_uncorrected_gene_level_pvalue_histogram <- function(pvalue_comparison_file, tissue_type) {
	df <- read.table(pvalue_comparison_file,header=TRUE)
	#df <- df[sample(nrow(df), 100000),]

	pvalues <- c()
	version <- c()

	pvalues <- c(pvalues, df$corrected_pvalue)
	version <- c(version, rep("corrected", length(df$corrected_pvalue)))
	pvalues <- c(pvalues, df$uncorrected_pvalue)
	version <- c(version, rep("uncorrected", length(df$uncorrected_pvalue)))


	df2 <- data.frame(pvalue=pvalues, version=factor(version))

	histo <- ggplot(df2, aes(x=pvalue, fill=version)) + geom_histogram(alpha=.35,position="identity",breaks = seq(0,1,.01)) + 
				labs(x = "pvalue", colour="", fill="") + 
				gtex_v8_figure_theme()
	return(histo)
}


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

temp_scatter <- function(cluster_id, gold_standard_p, alt_p, output_dir) {
	output_file <- paste0(output_dir, cluster_id, "pvalue_scatter.pdf")
	df <- data.frame(standard=-log10(gold_standard_p + 1e-7), alt=-log10(alt_p + 1e-7))
	scatter <- ggplot(df, aes(x=standard, y=alt)) + geom_point(size=.1) + 
				labs(x = "-log10(standard pvalue)", y = "-log10(alt pvalue)") +
				gtex_v8_figure_theme()	
	ggsave(scatter, file=output_file, width=7.2, height=5, units="in")

}


robustness_of_outlier_calls_to_emperical_read_depth_boxplot <- function(tissue_type, covariate_method, splicing_outlier_dir) {
	gold_standard_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_num_read_hyperparam_20000_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
	gold_standard_data <- read.table(gold_standard_file, header=TRUE)
	gold_standard_cluster_ids <- as.character(gold_standard_data[,1])
	gold_standard_data <- as.matrix(gold_standard_data[,2:dim(gold_standard_data)[2]])
	
	num_read_arr <- c()
	correlation_arr <- c()
	num_read_array = c("1000")
	temp <- read.table("bad_clusters.txt", header=FALSE)
	frac <- temp[,1]

	for (num_read_iter in 1:length(num_read_array)) {
		num_reads <- num_read_array[num_read_iter]
		print(num_reads)
		temp_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_num_read_hyperparam_", num_reads, "_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
		temp_data <- read.table(temp_file, header=TRUE)
		temp_cluster_ids <- as.character(temp_data[,1])
		temp_data <- as.matrix(temp_data[,2:dim(temp_data)[2]])
		if (sum(temp_cluster_ids != gold_standard_cluster_ids) != 0.0) {
			print('assumption error')
		}
		num_clusters <- dim(gold_standard_data)[1]
		counter = 0
		for (cluster_iter in 1:num_clusters) {
			cluster_id <- gold_standard_cluster_ids[cluster_iter]
			diff <- (-log10(gold_standard_data[cluster_iter,] + 1e-6)) - (-log10(temp_data[cluster_iter,] + 1e-6))


			indices <- (diff > 2) & (-log10(gold_standard_data[cluster_iter,] + 1e-6)) < 3

			if (sum(indices) > 0) {
				print(cluster_id)
			}
			if (cluster_id == "cluster10046") {
				print(diff)
			}


			corry <- cor(gold_standard_data[cluster_iter,], temp_data[cluster_iter,], method="pearson")

			#diff <- max(abs(gold_standard_data[cluster_iter,] - temp_data[cluster_iter,]))
			correlation_arr <- c(correlation_arr, corry)
			num_read_arr <- c(num_read_arr, num_reads)
		}
		print(counter)
		print(num_clusters)
		print(counter/num_clusters)
	}

	df <- data.frame(correlation=correlation_arr, num_reads=factor(num_read_arr, levels=num_read_array))
	print(summary(df))
	p <- ggplot(df, aes(x=num_reads, y=correlation, fill=num_reads)) + 
  			geom_boxplot() +
  			gtex_v8_figure_theme() +
  			labs(x = "Number of simulated Reads", y="Spearman correlation", fill="") + 
  			theme(legend.position="none")
  	return(p)
}

robustness_of_outlier_calls_to_emperical_read_depth_scatterplot <- function(tissue_type, covariate_method, splicing_outlier_dir) {
	gold_standard_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_num_read_hyperparam_20000_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
	gold_standard_data <- read.table(gold_standard_file, header=TRUE)
	gold_standard_cluster_ids <- as.character(gold_standard_data[,1])
	gold_standard_data <- as.matrix(gold_standard_data[,2:dim(gold_standard_data)[2]])
	frac <- as.matrix(read.table(paste0(splicing_outlier_dir, tissue_type, "_max_fraction_of_reads_from_a_junction.txt"), header=FALSE))

	#indices <- rep("hi", dim(gold_standard_data)[1])
	#indices <- indices=="hi"
	#indices[4937] = FALSE
	alt_pvalue <- c()
	gs_pvalue <- c()
	frac_arr <- c()
	num_read_array = c("100000")

	for (num_read_iter in 1:length(num_read_array)) {
		num_reads <- num_read_array[num_read_iter]
		temp_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_num_read_hyperparam_", num_reads, "_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
		temp_data <- read.table(temp_file, header=TRUE)
		temp_cluster_ids <- as.character(temp_data[,1])
		temp_data <- as.matrix(temp_data[,2:dim(temp_data)[2]])
		if (sum(temp_cluster_ids != gold_standard_cluster_ids) != 0.0) {
			print('assumption error')
		}
		print("START")
		alt_pvalue = as.vector(temp_data)
		gs_pvalue = as.vector(gold_standard_data)
		frac_arr = as.vector(frac)

		num_clusters <- dim(gold_standard_data)[1]

	}
	frac_arr[frac_arr < .8] = .8
	random_indices <- sample (c(1:length(gs_pvalue)), size=length(gs_pvalue)/20.0, replace=F)
	print(cor(gs_pvalue, alt_pvalue))

	diff = abs(-log10(gs_pvalue+1e-7) - -log10(alt_pvalue+1e-7))
	print(sum(diff > 1))
	print(length(diff))
	print(sum(diff > 1)/length(diff))

	gs_pvalue = gs_pvalue[random_indices]
	alt_pvalue = alt_pvalue[random_indices]
	frac_arr = frac_arr[random_indices]

	df <- data.frame(stadard_pvalue=-log10(gs_pvalue+1e-7), alt_pvalue=-log10(alt_pvalue+1e-7), fraction=frac_arr)
	p <- ggplot(df, aes(x=stadard_pvalue, y=alt_pvalue,colour=frac_arr)) + 
  			geom_point(size=.0001) +
  			gtex_v8_figure_theme() +
  			labs(x = "-log10(pvalue gs)", y="-log10(pvalue alt)",colour="frac_arr") 
  	return(p)
}

compare_outlier_distributions <- function(gold_standard_file, comparison_file, fraction_of_reads_file, x_axis_label, y_axis_label) {
	# Load in Gold Standard outlier data
	gold_standard_data <- read.table(gold_standard_file, header=TRUE)
	gold_standard_cluster_ids <- as.character(gold_standard_data[,1])
	gold_standard_data <- as.matrix(gold_standard_data[,2:dim(gold_standard_data)[2]])

	# Load in comparison outlier data
	comparison_data <- read.table(comparison_file, header=TRUE)
	comparison_cluster_ids <- as.character(comparison_data[,1])
	comparison_data <- as.matrix(comparison_data[,2:dim(comparison_data)[2]])

	frac <- as.matrix(read.table(fraction_of_reads_file, header=FALSE))

	# Subset to clusters found in b
	indices = gold_standard_cluster_ids %in% comparison_cluster_ids
	num_clusters_failed = length(indices) - sum(indices)
	print(paste0(num_clusters_failed, " clusters did not converge when no prior was used"))
	gold_standard_data <- gold_standard_data[indices,]
	frac <- frac[indices,]
	gold_standard_cluster_ids <- gold_standard_cluster_ids[indices]

	# Double check that data is matching
	if (sum(comparison_cluster_ids != gold_standard_cluster_ids) != 0.0) {
		print('assumption error')
	}

		
	alt_pvalue = as.vector(comparison_data)
	gs_pvalue = as.vector(gold_standard_data)
	frac_arr = as.vector(frac)
	frac_arr[frac_arr < .8] = .8

	# OUTPUT SOME STATISTICS
	print(cor(gs_pvalue, alt_pvalue, method="spearman"))
	diff = abs(-log10(gs_pvalue+1e-6) - -log10(alt_pvalue+1e-6))
	print(sum(diff > 1))
	print(length(diff))
	print(sum(diff > 1)/length(diff))


	# randomly sample indices for viz
	random_indices <- sample (c(1:length(gs_pvalue)), size=length(gs_pvalue)/100.0, replace=F)
	gs_pvalue = gs_pvalue[random_indices]
	alt_pvalue = alt_pvalue[random_indices]
	frac_arr = frac_arr[random_indices]

	# Put into compact data frame
	df <- data.frame(stadard_pvalue=-log10(gs_pvalue+1e-6), alt_pvalue=-log10(alt_pvalue+1e-6), fraction=frac_arr)
	# Plot
	p <- ggplot(df, aes(x=stadard_pvalue, y=alt_pvalue,colour=frac_arr)) + 
  			geom_point(size=.0001) +
  			gtex_v8_figure_theme() +
  			theme(legend.position="bottom") + 
  			labs(x = paste0("-log10(p-value) [ ", x_axis_label, " ]"), y=paste0("-log10(p-value) [ ", y_axis_label, " ]"),colour="Fraction of reads\nfrom one junction") 
  	return(p)
}

robustness_of_outlier_calls_to_prior_scatterplot <- function(tissue_type, covariate_method, splicing_outlier_dir) {
	
	prior_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_model_version_hyperparam_standard_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
	prior_data <- read.table(prior_file, header=TRUE)
	prior_cluster_ids <- as.character(prior_data[,1])
	prior_data <- as.matrix(prior_data[,2:dim(prior_data)[2]])

	no_prior_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_model_version_hyperparam_no_prior_multiple_initializations_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
	no_prior_data <- read.table(no_prior_file, header=TRUE)
	no_prior_cluster_ids <- as.character(no_prior_data[,1])
	no_prior_data <- as.matrix(no_prior_data[,2:dim(no_prior_data)[2]])

	print(dim(prior_data))
	print(dim(no_prior_data))
	# Many clusters were not completed when no_prior was used
	# Subset prior to only those in no_prior
	indices = prior_cluster_ids %in% no_prior_cluster_ids
	
	num_clusters_failed = length(indices) - sum(indices)
	print(paste0(num_clusters_failed, " clusters did not converge when no prior was used"))

	prior_data <- prior_data[indices,]


	num_clusters <- dim(prior_data)[1]
	for (cluster_num in 1:num_clusters) {
		cluster_id <- no_prior_cluster_ids[cluster_num]
		prior <- -log10(prior_data[cluster_num,] +1e-6)
		no_prior <- -log10(no_prior_data[cluster_num,]+1e-6)
		diff = no_prior - prior
		if(sum(diff > 5) > 1) {
			print(cluster_id)
		}
	}

	# Put data in vector format

	prior_pvalue = as.vector(prior_data)
	no_prior_pvalue = as.vector(no_prior_data)
	print(cor(no_prior_pvalue, prior_pvalue))



	random_indices <- sample(c(1:length(prior_pvalue)), size=length(prior_pvalue)/20.0, replace=F)
	prior_pvalue = prior_pvalue[random_indices]
	no_prior_pvalue = no_prior_pvalue[random_indices]


	df <- data.frame(prior_pvalue=-log10(prior_pvalue+1e-6), no_prior_pvalue=-log10(no_prior_pvalue+1e-6))
	p <- ggplot(df, aes(x=prior_pvalue, y=no_prior_pvalue)) + 
  			geom_point(size=.0001) +
  			gtex_v8_figure_theme() +
  			labs(x = "-log10(pvalue) (prior used)", y="-log10(pvalue alt) (no prior used)") 
  	return(p)



}

#######################
# Command line args
#######################

tissue_names_file = args[1]
covariate_method = args[2]
splicing_outlier_dir = args[3]
splicing_outlier_visualization_dir = args[4]

options(bitmapType = 'cairo', device = 'pdf')

# Extract vector tissue names
tissue_names <- as.character(unlist(read.table(tissue_names_file,header=FALSE), use.names=FALSE))


#cross_tissue_outlier_file <- paste0(splicing_outlier_dir, "cross_tissue_covariate_method_", covariate_method, "_no_global_outliers_ea_only_emperical_pvalue.txt")
#cross_tissue_outliers <- read.table(cross_tissue_outlier_file, header=TRUE)
#clusters <- cross_tissue_outliers[,1]
#cross_tissue_outliers <- cross_tissue_outliers[,2:dim(cross_tissue_outliers)[2]]

# Make bar plot of number outliers at various pvalue thresholds
output_file <- paste0(splicing_outlier_visualization_dir, "number_of_multi_tissue_outliers_ea_only_at_various_pvalue_thresholds_bar_plot.pdf")
#bar_plot_showing_number_of_multi_tissue_outliers_at_various_pvalue_thresholds(cross_tissue_outliers, output_file)



# Plot histogram of number of multi-tissue outliers per individual
number_of_outliers_file <- paste0(splicing_outlier_dir, "number_of_multi_tissue_outliers_ea_only.txt")  # Input file
output_file <- paste0(splicing_outlier_visualization_dir, "number_of_multi_tissue_outliers_ea_only_histogram.pdf")
#histogram_showing_number_of_multi_tissue_outlier_calls(number_of_outliers_file, output_file)

# Plot histogram of number of multi-tissue outliers per individual (not just ea only)
number_of_outliers_file <- paste0(splicing_outlier_dir, "number_of_multi_tissue_outliers.txt")  # Input file
output_file <- paste0(splicing_outlier_visualization_dir, "number_of_multi_tissue_outliers_histogram.pdf")
#histogram_showing_number_of_multi_tissue_outlier_calls(number_of_outliers_file, output_file)

###############################
# Cowplot supplementary figure describing gene level correction in Muscle
###############################

# Scatter plot showing corrected and uncorrected -log10(pvalues) at gene level colored by number of clusters mapped to the gene
tissue_type="Muscle_Skeletal"
pvalue_comparison_file <- paste0(splicing_outlier_dir, tissue_type, "_covariate_method_", covariate_method,"_gene_level_method_comparison.txt")
# cluster_correction_scatter <- corrected_vs_uncorrected_gene_level_pvalue_scatter_plot(pvalue_comparison_file, tissue_type)

# histogram showing corrected and uncorrected pvalues at gene level
tissue_type="Muscle_Skeletal"
pvalue_comparison_file <- paste0(splicing_outlier_dir, tissue_type, "_covariate_method_", covariate_method,"_gene_level_method_comparison.txt")
# clustter_correction_histogram <- corrected_vs_uncorrected_gene_level_pvalue_histogram(pvalue_comparison_file, tissue_type)

output_file <- paste0(splicing_outlier_visualization_dir, "corrected_vs_uncorrected_gene_level_pvalues_", tissue_type, "_joint.pdf")
#joint_correction_plot <- plot_grid(cluster_correction_scatter, clustter_correction_histogram, ncol=1, labels=c("A","B"))
# ggsave(joint_correction_plot, file=output_file, width=7.2, height=5, units="in")


###############################
# Robustness of splicing outlier calls in Muscle-Skeletal to hyperparameters
###############################

# Compare outlier distrbutions as a function of prior choice
tissue_type <- "Muscle_Skeletal"
output_file <- paste0(splicing_outlier_visualization_dir, "compare_outlier_distributions_prior_vs_no_prior_scatter.pdf")
gold_standard_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_model_version_hyperparam_standard_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
comparison_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_model_version_hyperparam_no_prior_multiple_initializations_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
fraction_of_reads_file <- paste0(splicing_outlier_dir, tissue_type, "_max_fraction_of_reads_from_a_junction.txt")
outlier_comparison_scatter_no_prior <- compare_outlier_distributions(gold_standard_file, comparison_file, fraction_of_reads_file, "Standard prior", paste0("No prior"))
ggsave(outlier_comparison_scatter_no_prior, file=output_file, width=7.2, height=5, units="in")


# Compare outlier distrbutions as a function of emperical read depth
tissue_type <- "Muscle_Skeletal"
num_reads <- "1000"
output_file <- paste0(splicing_outlier_visualization_dir, "compare_outlier_distributions_20000_vs_", num_reads, "_scatter.pdf")
gold_standard_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_num_read_hyperparam_20000_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
comparison_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_num_read_hyperparam_", num_reads, "_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
fraction_of_reads_file <- paste0(splicing_outlier_dir, tissue_type, "_max_fraction_of_reads_from_a_junction.txt")
#outlier_comparison_scatter_1000_reads <- compare_outlier_distributions(gold_standard_file, comparison_file, fraction_of_reads_file, "20000 reads", paste0(num_reads, " reads"))
#ggsave(outlier_comparison_scatter_1000_reads, file=output_file, width=7.2, height=5, units="in")

tissue_type <- "Muscle_Skeletal"
num_reads <- "10000"
output_file <- paste0(splicing_outlier_visualization_dir, "compare_outlier_distributions_20000_vs_", num_reads, "_scatter.pdf")
gold_standard_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_num_read_hyperparam_20000_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
comparison_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_num_read_hyperparam_", num_reads, "_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
fraction_of_reads_file <- paste0(splicing_outlier_dir, tissue_type, "_max_fraction_of_reads_from_a_junction.txt")
outlier_comparison_scatter_10000_reads <- compare_outlier_distributions(gold_standard_file, comparison_file, fraction_of_reads_file, "20000 reads", paste0(num_reads, " reads"))
ggsave(outlier_comparison_scatter_10000_reads, file=output_file, width=7.2, height=5, units="in")

tissue_type <- "Muscle_Skeletal"
num_reads <- "100000"
output_file <- paste0(splicing_outlier_visualization_dir, "compare_outlier_distributions_20000_vs_", num_reads, "_scatter.pdf")
gold_standard_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_num_read_hyperparam_20000_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
comparison_file <- paste0(splicing_outlier_dir, tissue_type, "_compare_num_read_hyperparam_", num_reads, "_covariate_method_", covariate_method, "_merged_emperical_pvalue.txt")
fraction_of_reads_file <- paste0(splicing_outlier_dir, tissue_type, "_max_fraction_of_reads_from_a_junction.txt")
outlier_comparison_scatter_100000_reads <- compare_outlier_distributions(gold_standard_file, comparison_file, fraction_of_reads_file, "20000 reads", paste0(num_reads, " reads"))
ggsave(outlier_comparison_scatter_100000_reads, file=output_file, width=7.2, height=5, units="in")

# Make combined plot
output_file <- paste0(splicing_outlier_visualization_dir, "compare_outlier_distributions_combined.pdf")
joint_outlier_robustness_plot <- plot_grid(outlier_comparison_scatter_10000_reads, outlier_comparison_scatter_100000_reads, outlier_comparison_scatter_no_prior, ncol=2, labels=c("A","B", "C"))
ggsave(joint_outlier_robustness_plot, file=output_file, width=7.2, height=7, units="in")














#######################
# OLD (NO LONGER USED)
########################


tissue_type <- "Muscle_Skeletal"
output_file <- paste0(splicing_outlier_visualization_dir, "pseudocount_emperical_read_depth_robustness_boxplot.pdf")
# boxplot <- robustness_of_outlier_calls_to_emperical_read_depth_boxplot(tissue_type, covariate_method, splicing_outlier_dir)
#ggsave(boxplot, file=output_file, width=7.2, height=5, units="in")


tissue_type <- "Muscle_Skeletal"
output_file <- paste0(splicing_outlier_visualization_dir, "emperical_read_depth_100000_filter_robustness_scatter.pdf")
#scatter <- robustness_of_outlier_calls_to_emperical_read_depth_scatterplot(tissue_type, covariate_method, splicing_outlier_dir)
#ggsave(scatter, file=output_file, width=7.2, height=5, units="in")


###############################
# Robustness of splicing outlier calls in prior hyperparameters
###############################
tissue_type <- "Muscle_Skeletal"
output_file <- paste0(splicing_outlier_visualization_dir, tissue_type, "_outlier_prior_robustness_scatter.pdf")
#scatter <- robustness_of_outlier_calls_to_prior_scatterplot(tissue_type, covariate_method, splicing_outlier_dir)
#ggsave(scatter, file=output_file, width=7.2, height=5, units="in")

