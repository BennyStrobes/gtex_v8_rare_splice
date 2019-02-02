args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)

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

extract_odds_ratio_data <- function(outlier_distances, inlier_distances, distance_window) {
	dist_to_ss <- c()
	odds_ratios <- c()
	upper_bounds <- c()
	lower_bounds <- c()
	for (distance_iter in -distance_window:(distance_window-1)) {
		rv_outlier <- sum(outlier_distances[,1] == distance_iter)
		rv_inlier <- sum(inlier_distances[,1] == distance_iter)
		no_rv_outlier <- sum(outlier_distances[,1] != distance_iter)
		no_rv_inlier <- sum(inlier_distances[,1] != distance_iter)
		# odds_ratio_old <- (rv_outlier/rv_inlier)/(no_rv_outlier/no_rv_inlier)
		a <- rv_outlier  + 1
		b <- rv_outlier + no_rv_outlier +1
		c <- rv_inlier +1
		d <- rv_inlier + no_rv_inlier +1
		orat <- (a/b)/(c/d)

		log_bounds <- 1.96*sqrt((1.0/a) - (1.0/b) + (1.0/c) - (1.0/d))
		#upper_bound <- exp(log_orat + log_bounds)
		#lower_bound <- exp(log_orat - log_bounds) 
		upper_bound <- orat*exp(log_bounds)
		lower_bound <- orat*exp(-log_bounds)
		if (orat == 0) {
			upper_bound = 0
			lower_bound = 0
		}

		# Add data to array
		odds_ratios <- c(odds_ratios, orat)
		dist_to_ss <- c(dist_to_ss, distance_iter)
		upper_bounds <- c(upper_bounds, upper_bound)
		lower_bounds <- c(lower_bounds, lower_bound)
	}
	df <- data.frame(odds_ratio=odds_ratios, dist_to_ss=dist_to_ss, lower_bounds=lower_bounds, upper_bounds=upper_bounds)

	return(df)

}



make_variant_allele_odds_ratio_plot <- function(distance, inlier_distance_file, outlier_distance_file, output_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)


	ss_type="acceptor"
	position=-2

	position_outlier_distances <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position,]
	position_inlier_distances <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance == position,]

	print(position_inlier_distances)
	alleles <- c("A", "C", "T", "G")

	variant_types <- c()
	odds_ratios <- c()
	upper_bounds <- c()
	lower_bounds <- c()

	# Loop through all possible variants
	for (iter1 in 1:length(alleles)) {
		for (iter2 in 1:length(alleles)) {
			major_allele = alleles[iter1]
			variant_allele = alleles[iter2]
			if (major_allele != variant_allele) {
				variant_type <- paste0(major_allele, "->",variant_allele)
				a <- sum(position_outlier_distances$major_allele == major_allele & position_outlier_distances$variant_allele == variant_allele) 
				b <- dim(position_outlier_distances)[1] 
				c <- sum(position_inlier_distances$major_allele == major_allele & position_inlier_distances$variant_allele == variant_allele) 
				d <- dim(position_inlier_distances)[1]
				orat <- (a/b)/(c/d)
				print(variant_type)
				print(a)
				print(b)
				print(c)
				print(d)

				log_bounds <- 1.96*sqrt((1.0/a) - (1.0/b) + (1.0/c) - (1.0/d))
				upper_bound <- orat*exp(log_bounds)
				lower_bound <- orat*exp(-log_bounds)

				variant_types <- c(variant_types, variant_type)
				odds_ratios <- c(odds_ratios, orat)
				lower_bounds <- c(lower_bounds, lower_bound)
				upper_bounds <- c(upper_bounds, upper_bound)
			}

		}
	}

	df <- data.frame(variant_types=factor(variant_types), odds_ratios=odds_ratios, lower_bounds=lower_bounds, upper_bounds=upper_bounds)
	print(df)
}




# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites) using odds ratios of real vs background
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
make_distance_odds_ratio_density_plot_seperated_by_ss_type <- function(distance_window, inlier_distance_file, outlier_distance_file, title, output_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	outlier_donor_distances <- outlier_distances[as.character(outlier_distances$splice_site_type)=="donor",]
	outlier_acceptor_distances <- outlier_distances[as.character(outlier_distances$splice_site_type)=="acceptor",]
	inlier_donor_distances <- inlier_distances[as.character(inlier_distances$splice_site_type)=="donor",]
	inlier_acceptor_distances <- inlier_distances[as.character(inlier_distances$splice_site_type)=="acceptor",]



	donor_df <- extract_odds_ratio_data(outlier_donor_distances, inlier_donor_distances, 10)
	acceptor_df <- extract_odds_ratio_data(outlier_acceptor_distances, inlier_acceptor_distances, 10)



	#Combine outliers and non-outliers into a compact data frame
	odds_ratio <- c(donor_df$odds_ratio, acceptor_df$odds_ratio)
	dist_to_ss <- c(donor_df$dist_to_ss, acceptor_df$dist_to_ss)
	ss_type <- c(rep("donor",length(donor_df$dist_to_ss)), rep("acceptor", length(acceptor_df$dist_to_ss)))
	lower_bounds <- c(donor_df$lower_bounds, acceptor_df$lower_bounds)
	upper_bounds <- c(donor_df$upper_bounds, acceptor_df$upper_bounds)

	df <- data.frame(odds_ratio=odds_ratio, dist_to_ss=dist_to_ss, lower_bound=lower_bounds, upper_bound=upper_bounds, ss_type=factor(ss_type, levels=c("donor","acceptor")))

	df_acceptor <- df[as.character(df$ss_type) == "acceptor",]
	df_donor <- df[as.character(df$ss_type) == "donor", ]

	options(bitmapType = 'cairo', device = 'pdf')

	acceptor_plot <-  ggplot() + geom_errorbar(data=df_acceptor, mapping=aes(x=dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="darkorchid") +
					geom_point(data=df_acceptor, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
					labs(x = "Distance from splice site (BP)", y = "Enrichment", title="Acceptor splice site") +
					geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
					geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					theme(axis.text.x=element_text(angle=45,hjust=1), text = element_text(size=11),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=11)) +
					scale_x_continuous(breaks=-10:9, labels=c("A-10","A-9","A-8","A-7","A-6","A-5","A-4","A-3","A-2","A-1","A+1","A+2","A+3","A+4","A+5","A+6", "A+7","A+8","A+9","A+10"))


	donor_plot <-  ggplot() + geom_errorbar(data=df_donor, mapping=aes(x=dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="darkorchid") +
					geom_point(data=df_donor, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
					labs(x = "Distance from splice site (BP)", y = "Enrichment", title="Donor splice site") +
					geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
					geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					theme(axis.text.x=element_text(angle=45,hjust=1),text = element_text(size=11),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=11)) +
					scale_x_continuous(breaks=-10:9, labels=c("D+10","D+9","D+8","D+7","D+6","D+5","D+4","D+3","D+2","D+1","D-1","D-2","D-3","D-4","D-5","D-6", "D-7","D-8","D-9","D-10"))


	error_bar_plot <- plot_grid(donor_plot,acceptor_plot, ncol=1)
    #error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="darkorchid") +
	#				geom_point(data=df, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
	#				facet_wrap( ~ ss_type,nrow=2,scales = "free") +
	#				labs(x = "Distance from splice site (BP)", y = "Enrichment") +
	#				geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
	#				geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
	#				geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
	#				theme(text = element_text(size=11),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=11))




	options(bitmapType = 'cairo', device = 'pdf')

	
	ggsave(error_bar_plot, file=output_file,width = 17,height=14,units="cm")
}










#######################
# Command Line args
#######################
variant_position_enrichment_dir <- args[1]  # Input dir
visualize_variant_position_enrichment_dir <- args[2]  # Output dir


#######################
# Specific allele plots
#########################
distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_variant_alleles_odds_ratio_plot.pdf")
# make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file)


#######################
# Positional distribution plots for RV
#########################
# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites) using odds ratios of real vs background
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_odds_ratio_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
make_distance_odds_ratio_density_plot_seperated_by_ss_type(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)

distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-06"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_odds_ratio_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
make_distance_odds_ratio_density_plot_seperated_by_ss_type(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)




