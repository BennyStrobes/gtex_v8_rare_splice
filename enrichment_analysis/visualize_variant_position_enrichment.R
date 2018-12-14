args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)





# Make density plot of distance between rare variants and splice sites 
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
make_distance_density_plot <- function(inlier_distance_file, outlier_distance_file, title, output_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)
	distances <- c()
	versions <- c()
	# Extract outlier distances
	num_outliers <- dim(outlier_distances)[1]
	#for (iter in 1:num_outliers) {
	#	distance <- outlier_distances[iter,1]
	#	distances <- c(distances, distance)
	#	versions <- c(versions, "outlier")
	#}
	# Extract inlier distances
	num_inliers <- dim(inlier_distances)[1]
	#for (iter in 1:num_inliers) {
#		distance <- inlier_distances[iter, 1]
	#	distances <- c(distances, distance)
	#	versions <- c(versions, "background")
	#}
	distances <- c(as.numeric(outlier_distances[,1]), as.numeric(inlier_distances[,1]))
	versions <- c(rep("outlier", num_outliers), rep("background", num_inliers))

	# Put into compact data frame
	df <- data.frame(distance=distances, version=factor(versions, levels=c("outlier", "background")))

	density_plot <- ggplot(df, aes(distance, colour = version, fill=version)) + geom_density(alpha=.1,adjust=3/4) +
		labs(x = "Distance from splice site (BP)", y = "Density", title=title,colour=" ", fill=" ") +
		#theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) +
		theme_bw() + theme(text = element_text(size=10),axis.text=element_text(size=9),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) +
		geom_hline(yintercept=0, colour="white", size=1)

	ggsave(density_plot, file=output_file,width = 16,height=10.5,units="cm")

}











#######################
# Command Line args
#######################
variant_position_enrichment_dir <- args[1]  # Input dir
visualize_variant_position_enrichment_dir <- args[2]  # Output dir



# Make density plot of distance between rare variants and splice sites 
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "20"
version <- "exon_annotation"
pvalue_thresholds <- c("1e-05", "0.0001", "0.001")
for (iter in 1:length(pvalue_thresholds)) {
	pvalue_threshold <- pvalue_thresholds[iter]
	outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
	inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
	output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "density_plot.pdf")
	title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
	make_distance_density_plot(inlier_distance_file, outlier_distance_file, title, output_file)
}

distance <- "100"
version <- "exon_annotation"
pvalue_thresholds <- c("1e-05", "0.0001", "0.001")
for (iter in 1:length(pvalue_thresholds)) {
	pvalue_threshold <- pvalue_thresholds[iter]
	outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
	inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
	output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "density_plot.pdf")
	title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
	make_distance_density_plot(inlier_distance_file, outlier_distance_file, title, output_file)
}

distance <- "20"
version <- "observed_splice_site"
pvalue_thresholds <- c("1e-05", "0.0001", "0.001")
for (iter in 1:length(pvalue_thresholds)) {
	pvalue_threshold <- pvalue_thresholds[iter]
	outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
	inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
	output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "density_plot.pdf")
	title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
	make_distance_density_plot(inlier_distance_file, outlier_distance_file, title, output_file)
}

distance <- "100"
version <- "observed_splice_site"
pvalue_thresholds <- c("1e-05", "0.0001", "0.001")
for (iter in 1:length(pvalue_thresholds)) {
	pvalue_threshold <- pvalue_thresholds[iter]
	outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
	inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
	output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "density_plot.pdf")
	title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
	make_distance_density_plot(inlier_distance_file, outlier_distance_file, title, output_file)
}
