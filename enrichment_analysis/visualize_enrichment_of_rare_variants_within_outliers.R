args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)






tbt_variant_outlier_enrichment_errorbar_plot <- function(tissue_by_tissue_enrichment_file, output_file,title, color_vector) {
	enrichments <- read.table(tissue_by_tissue_enrichment_file, header=FALSE)
	num_tissues <- dim(enrichments)[1]

	# Initialize vectors
	tissue_names <- c()
	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()

	# Loop through tissues
	for (tissue_number in 1:num_tissues) {
		# Extract tissue name for this line
		tissue_name <- as.character(enrichments[tissue_number, 1])
		# Compute odds ratios for this line
		a <- enrichments[tissue_number, 2]
		b <- enrichments[tissue_number, 3]
		c <- enrichments[tissue_number, 4]
		d <- enrichments[tissue_number, 5]
		orat <- (a/b)/(c/d)
		# Compute error bars for this orat
		log_orat <- log(orat)
		log_bounds <- 1.96*sqrt((1.0/a) - (1.0/b) + (1.0/c) - (1.0/d))
		#upper_bound <- exp(log_orat + log_bounds)
		#lower_bound <- exp(log_orat - log_bounds) 
		upper_bound <- orat*exp(log_bounds)
		lower_bound <- orat*exp(-log_bounds)

		# Add information to vectors
		tissue_names <- c(tissue_names, tissue_name)
		odds_ratios <- c(odds_ratios, orat)
		lower_bounds <- c(lower_bounds, lower_bound)
		upper_bounds <- c(upper_bounds, upper_bound)
	}

	# Add information to data frame
	df <- data.frame(tissue_names=factor(tissue_names), odds_ratios=odds_ratios, lower_bounds=lower_bounds, upper_bounds=upper_bounds)
	error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=tissue_names,ymin=lower_bounds, ymax=upper_bounds), colour=color_vector) +
					geom_point(data=df, mapping=aes(x=tissue_names, y=odds_ratios), colour=color_vector) +
					labs(x = "Tissue", y = "Enrichment", title=title) +
					theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5)) 

	ggsave(error_bar_plot, file=output_file,width = 20,height=12.5,units="cm")

}




get_color_vector <- function(tissue_colors, tissue_names) {
	colors <- c()
	for (iter in 1:length(tissue_names)) {
		tissue_name <- as.character(tissue_names[iter])
		if (tissue_name == "Brain_Spinal_cord_cervical_c-1") {
			tissue_name = "Brain_Spinal_cord_cervical_c1"
		}
		if (tissue_name == "Cells_EBV-transformed_lymphocytes") {
			tissue_name = "Cells_EBVtransformed_lymphocytes"
		}

		indices = tissue_name==tissue_colors$tissue_id

		color <- tissue_colors$tissue_color_hex[indices]
		
		colors <- c(colors, paste0('#',color))
	}
	return(colors)
}

cross_tissue_variant_outlier_enrichment_errorbar_plot <- function(variant_enrichment_dir, output_file, version) {
	distances <- c("4", "6", "8", "10", "100", "1000")

	# Initialize vectors
	pvalues <- c()
	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()
	distance_vec <- c()

	# Loop through distances (one enrichment file for each distance)
	for (distance_iter in 1:length(distances)) {
		distance <- distances[distance_iter]
		file_name <- paste0(variant_enrichment_dir, "cross_tissue_variant_outlier_enrichment_distance_", distance, "_version_", version, ".txt")
		enrichments <- read.table(file_name, header=FALSE)
		num_pvalues <- dim(enrichments)[1]
		pvalue_string <- c()
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
			distance_vec <- c(distance_vec, distance)
		}
	}

	# Make data frame
	df <- data.frame(pvalues=factor(pvalues,levels=(as.character(pvalue_string))), odds_ratios=odds_ratios, lower_bounds=lower_bounds, upper_bounds=upper_bounds,distance=factor(distance_vec,levels=distances))
	dodge <- position_dodge(width=0.9)
	error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=distance,ymin=lower_bounds, ymax=upper_bounds, colour=pvalues), position=dodge) +
					geom_point(data=df, mapping=aes(x=distance, y=odds_ratios, colour=pvalues), position=dodge) +
					labs(x = "Distance", y = "Enrichment", colour="p-value") +
					theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) 

	ggsave(error_bar_plot, file=output_file,width = 19,height=10.5,units="cm")

}





###############################
# Command Line Args
###############################
visualize_variant_enrichment_dir = args[1]  # Output Directory
variant_enrichment_dir = args[2]  # Input directory
tissue_names_file = args[3]  # File containing names of gtex tissues
tissue_colors_file = args[4] # File containing mapping from gtex tissue names to tissue colors

# Extract vector tissue names
tissue_names <- as.character(unlist(read.table(tissue_names_file,header=FALSE), use.names=FALSE))

tissue_colors = read.table(tissue_colors_file, header = T, stringsAsFactors = F, sep = "\t")

# Get vector of hex colors in correct order
color_vector <- get_color_vector(tissue_colors, tissue_names)


distance <- '8'
pvalue_threshold <- '.000001'
version <- "all"

tissue_by_tissue_enrichment_file <- paste0(variant_enrichment_dir, "tbt_variant_outlier_enrichment_pvalue_", pvalue_threshold, "_distance_", distance, "_version_",version,".txt")
output_file <- paste0(visualize_variant_enrichment_dir, "tbt_variant_outlier_enrichment_pvalue_", pvalue_threshold, "_distance_", distance, "_version_",version, "_errorbar.pdf")
title <- paste0("distance=",distance, " / pvalue=", pvalue_threshold)
tbt_variant_outlier_enrichment_errorbar_plot(tissue_by_tissue_enrichment_file, output_file, title, color_vector)


output_file <- paste0(visualize_variant_enrichment_dir, "cross_tissue_variant_outlier_enrichment_version_",version,"_errorbar.pdf")
cross_tissue_variant_outlier_enrichment_errorbar_plot(variant_enrichment_dir, output_file, version)

