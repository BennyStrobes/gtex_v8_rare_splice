args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)

extract_data_no_splice_site_type <- function(distance_object, distance_window) {
	dist_to_ss <- c()
	density <- c()
	distances <- distance_object[,1]
	for (distance_iter in -distance_window:distance_window) {
		temp_density = sum(distances == distance_iter)/length(distances)
		density <- c(density, temp_density)
		dist_to_ss <- c(dist_to_ss, distance_iter)
	}
	temp_density <- sum(distances == distance_window)/length(distances)
	density <- c(density, temp_density)
	dist_to_ss <- c(dist_to_ss, distance_window + 1)
	df <- data.frame(density=density, dist_to_ss=dist_to_ss-.5)
	return(df)
}


# Make density plot of distance between rare variants and splice sites 
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
make_distance_density_plot <- function(distance_window, inlier_distance_file, outlier_distance_file, title, output_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	outlier_df <- extract_data_no_splice_site_type(outlier_distances, distance_window)
	inlier_df <- extract_data_no_splice_site_type(inlier_distances, distance_window)

	num_rows <- distance_window*2 + 2

	# Combine outliers and non-outliers into a compact data frame
	density <- c(outlier_df$density, inlier_df$density)
	dist_to_ss <- c(outlier_df$dist_to_ss, inlier_df$dist_to_ss)
	version <- c(rep("outlier", num_rows), rep("background", num_rows))

	df <- data.frame(density=density, dist_to_ss=dist_to_ss, version=factor(version, levels=c("outlier","background")))
	density_plot <- ggplot(df, aes(x=dist_to_ss,y=density,colour=version)) + geom_step(size=1.5) +
		scale_color_manual(values=c("darkorchid", "grey42")) +
		labs(x = "Distance from splice site (BP)", y = "Density", title=title,colour=" ", fill=" ") +
		theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) +
		geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
		geom_vline(xintercept = -2.5, size=.00001,linetype="dashed")
	
	ggsave(density_plot, file=output_file,width = 16,height=10.5,units="cm")


}




# Make density plot of distance between rare variants and splice sites 
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
make_distance_density_plot_seperated_by_ss_type <- function(distance_window, inlier_distance_file, outlier_distance_file, title, output_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)


	outlier_donor_distances <- outlier_distances[outlier_distances$splice_site_type=="donor",]
	outlier_acceptor_distances <- outlier_distances[outlier_distances$splice_site_type=="acceptor",]
	inlier_donor_distances <- inlier_distances[inlier_distances$splice_site_type=="donor",]
	inlier_acceptor_distances <- inlier_distances[inlier_distances$splice_site_type=="acceptor",]

	outlier_donor_df <- extract_data_no_splice_site_type(outlier_donor_distances, distance_window)
	outlier_acceptor_df <- extract_data_no_splice_site_type(outlier_acceptor_distances, distance_window)
	inlier_donor_df <- extract_data_no_splice_site_type(inlier_donor_distances, distance_window)
	inlier_acceptor_df <- extract_data_no_splice_site_type(inlier_acceptor_distances, distance_window)


	num_rows <- distance_window*2 + 2

	#Combine outliers and non-outliers into a compact data frame
	density <- c(outlier_donor_df$density, inlier_donor_df$density, outlier_acceptor_df$density, inlier_acceptor_df$density)
	dist_to_ss <- c(outlier_donor_df$dist_to_ss, inlier_donor_df$dist_to_ss, outlier_acceptor_df$dist_to_ss, inlier_acceptor_df$dist_to_ss)
	version <- c(rep("outlier", num_rows), rep("background", num_rows), rep("outlier", num_rows), rep("background", num_rows))
	ss_type <- c(rep("donor",num_rows*2), rep("acceptor", num_rows*2))

	df <- data.frame(density=density, dist_to_ss=dist_to_ss, ss_type=factor(ss_type, levels=c("donor","acceptor")), version=factor(version, levels=c("outlier","background")))

    #p <- ggplot(melted, aes(x=X, y=value, color=CellLineCluster)) + geom_line() + facet_wrap(~ L, ncol=10) 
    #if (ci){
    #    p <- p + geom_ribbon(data=data, aes(ymin=lower,ymax=upper, fill=CellLineCluster),alpha=0.3, colour=NA)
    #}
    #p <- p + labs(x="Day",y = "Expression",fill="",color="") + theme(legend.position = "bottom") + scale_fill_manual(values=c("dodgerblue3","chartreuse4"))+ scale_color_manual(values=c("dodgerblue3","chartreuse4"))
    #p <- p + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    #return(p)

	density_plot <- ggplot(df, aes(x=dist_to_ss,y=density,colour=version)) + geom_step(size=1.5) + 
		facet_wrap( ~ ss_type) +
		scale_color_manual(values=c("darkorchid", "grey42")) +
		labs(x = "Distance from splice site (BP)", y = "Density", title=title,colour=" ", fill=" ") +
		theme(text = element_text(size=11),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=11)) +
		theme(legend.position="bottom") +
		geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
		geom_vline(xintercept = -2.5, size=.00001,linetype="dashed")
	
	ggsave(density_plot, file=output_file,width = 17,height=11,units="cm")


}











#######################
# Command Line args
#######################
variant_position_enrichment_dir <- args[1]  # Input dir
visualize_variant_position_enrichment_dir <- args[2]  # Output dir



# Make density plot of distance between rare variants and splice sites 
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "10"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
make_distance_density_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)

# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites)
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "10"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
make_distance_density_plot_seperated_by_ss_type(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)



