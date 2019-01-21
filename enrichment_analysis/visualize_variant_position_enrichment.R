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

extract_odds_ratio_data <- function(outlier_distances, inlier_distances, distance_window) {
	dist_to_ss <- c()
	odds_ratios <- c()
	upper_bounds <- c()
	lower_bounds <- c()
	for (distance_iter in -distance_window:distance_window) {
		rv_outlier <- sum(outlier_distances[,1] == distance_iter)
		rv_inlier <- sum(inlier_distances[,1] == distance_iter)
		no_rv_outlier <- sum(outlier_distances[,1] != distance_iter)
		no_rv_inlier <- sum(inlier_distances[,1] != distance_iter)
		# odds_ratio_old <- (rv_outlier/rv_inlier)/(no_rv_outlier/no_rv_inlier)
		a <- rv_outlier
		b <- rv_outlier + no_rv_outlier
		c <- rv_inlier
		d <- rv_inlier + no_rv_inlier
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
	#rv_outlier <- sum(outlier_distances[,1] == distance_window)
	#rv_inlier <- sum(inlier_distances[,1] == distance_window)
	#no_rv_outlier <- sum(outlier_distances[,1] != distance_window)
	#no_rv_inlier <- sum(inlier_distances[,1] != distance_window)
	#odds_ratio <- (rv_outlier/rv_inlier)/(no_rv_outlier/no_rv_inlier)
	#odds_ratios <- c(odds_ratios, odds_ratio)
	#dist_to_ss <- c(dist_to_ss, distance_window + 1)
	#df <- data.frame(odds_ratio=odds_ratios, dist_to_ss=dist_to_ss-.5)
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

# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites) using odds ratios of real vs background
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
make_distance_odds_ratio_density_plot_seperated_by_ss_type <- function(distance_window, inlier_distance_file, outlier_distance_file, title, output_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)


	outlier_donor_distances <- outlier_distances[outlier_distances$splice_site_type=="donor",]
	outlier_acceptor_distances <- outlier_distances[outlier_distances$splice_site_type=="acceptor",]
	inlier_donor_distances <- inlier_distances[inlier_distances$splice_site_type=="donor",]
	inlier_acceptor_distances <- inlier_distances[inlier_distances$splice_site_type=="acceptor",]



	donor_df <- extract_odds_ratio_data(outlier_donor_distances, inlier_donor_distances, distance_window)
	acceptor_df <- extract_odds_ratio_data(outlier_acceptor_distances, inlier_acceptor_distances, distance_window)



	#Combine outliers and non-outliers into a compact data frame
	odds_ratio <- c(donor_df$odds_ratio, acceptor_df$odds_ratio)
	dist_to_ss <- c(donor_df$dist_to_ss, acceptor_df$dist_to_ss)
	ss_type <- c(rep("donor",length(donor_df$dist_to_ss)), rep("acceptor", length(acceptor_df$dist_to_ss)))
	lower_bounds <- c(donor_df$lower_bounds, acceptor_df$lower_bounds)
	upper_bounds <- c(donor_df$upper_bounds, acceptor_df$upper_bounds)

	df <- data.frame(odds_ratio=odds_ratio, dist_to_ss=dist_to_ss, lower_bound=lower_bounds, upper_bound=upper_bounds, ss_type=factor(ss_type, levels=c("donor","acceptor")))

    error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="darkorchid") +
					geom_point(data=df, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
					facet_wrap( ~ ss_type,nrow=2,scales = "free") +
					labs(x = "Distance from splice site (BP)", y = "Enrichment") +
					geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
					geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					theme(text = element_text(size=11),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=11))

	#density_plot <- ggplot(df, aes(x=dist_to_ss,y=odds_ratio)) + geom_step(size=1.5,color="darkorchid") + 
	#	facet_wrap( ~ ss_type) +
	#	labs(x = "Distance from splice site (BP)", y = "Odds ratio", title=title) +
	#	theme(text = element_text(size=11),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=11)) +
	#	theme(legend.position="bottom") +
	#	geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
	#	geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
	#	geom_hline(yintercept = 1, size=.00001,linetype="dashed")
	
	ggsave(error_bar_plot, file=output_file,width = 17,height=14,units="cm")
}



concensus_splice_site_plot <- function(inlier_distances, outlier_distances, consensus_allele, position, ss_type, output_file) {
	# Filter distance objects to only contain the given ss type and the given position
	inlier_distances = inlier_distances[inlier_distances$splice_site_type==ss_type,]
	inlier_distances = inlier_distances[inlier_distances$distance == position,]
	outlier_distances = outlier_distances[outlier_distances$splice_site_type==ss_type,]
	outlier_distances = outlier_distances[outlier_distances$distance == position,]

	inlier_total = dim(inlier_distances)[1]

	outlier_total = dim(outlier_distances)[1]

	fraction <- c()  
	version <- c()  # outlier or inlier
	class <- c() # "to_concensus", "from_concensus", "neither"

	to_concensus_outlier_counts <- sum(outlier_distances$alt == consensus_allele)/outlier_total
	from_concensus_outlier_counts <- sum(outlier_distances$ref == consensus_allele)/outlier_total
	neither_outlier_counts <- sum(outlier_distances$ref != consensus_allele & outlier_distances$alt != consensus_allele)/outlier_total

	to_concensus_inlier_counts <- sum(inlier_distances$alt == consensus_allele)/inlier_total
	from_concensus_inlier_counts <- sum(inlier_distances$ref == consensus_allele)/inlier_total
	neither_inlier_counts <- sum(inlier_distances$ref != consensus_allele & inlier_distances$alt != consensus_allele)/inlier_total

	# Add one at a time
	fraction <- c(fraction, to_concensus_outlier_counts)
	version <- c(version, "outlier")
	class <- c(class, "to_concensus")

	fraction <- c(fraction, from_concensus_outlier_counts)
	version <- c(version, "outlier")
	class <- c(class, "from_concensus")


	fraction <- c(fraction, neither_outlier_counts)
	version <- c(version, "outlier")
	class <- c(class, "neither")

	fraction <- c(fraction, to_concensus_inlier_counts)
	version <- c(version, "background")
	class <- c(class, "to_concensus")

	fraction <- c(fraction, from_concensus_inlier_counts)
	version <- c(version, "background")
	class <- c(class, "from_concensus")


	fraction <- c(fraction, neither_inlier_counts)
	version <- c(version, "background")
	class <- c(class, "neither")

	title <- paste0(ss_type, ' ', position)
	df <- data.frame(fraction=fraction, version= factor(version,levels=c("outlier","background")), class=factor(class, levels=c("to_concensus","from_concensus","neither")))

	bar_plot <- ggplot(df, aes(x=version,y=fraction,fill=class,colour=class)) + geom_bar(stat="identity") +
		labs(x = "", y = "Fraction", title=title,colour=" ", fill=" ") +
		theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10))

	
	ggsave(bar_plot, file=output_file,width = 16,height=10.5,units="cm")

}


make_consensus_splice_site_plots <- function(inlier_distance_file, outlier_distance_file, output_root) {
	# Load in data
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	# Make for Donor -1 (concensus='G')
	consensus_allele <- "G"
	position <- -1
	ss_type <- "donor"
	output_file <- paste0(output_root, ss_type, "_", position, "_concensus_rate_plot.pdf")
	concensus_splice_site_plot(inlier_distances, outlier_distances, consensus_allele, position, ss_type, output_file)

	# Make for Donor -2 (concensus='T')
	consensus_allele <- "T"
	position <- -2
	ss_type <- "donor"
	output_file <- paste0(output_root, ss_type, "_", position, "_concensus_rate_plot.pdf")
	concensus_splice_site_plot(inlier_distances, outlier_distances, consensus_allele, position, ss_type, output_file)

	# Make for Acceptor -1 (concensus='G')
	consensus_allele <- "G"
	position <- -1
	ss_type <- "acceptor"
	output_file <- paste0(output_root, ss_type, "_", position, "_concensus_rate_plot.pdf")
	concensus_splice_site_plot(inlier_distances, outlier_distances, consensus_allele, position, ss_type, output_file)

	# Make for Acceptor -2 (concensus='A')
	consensus_allele <- "A"
	position <- -2
	ss_type <- "acceptor"
	output_file <- paste0(output_root, ss_type, "_", position, "_concensus_rate_plot.pdf")
	concensus_splice_site_plot(inlier_distances, outlier_distances, consensus_allele, position, ss_type, output_file)


}










#######################
# Command Line args
#######################
variant_position_enrichment_dir <- args[1]  # Input dir
visualize_variant_position_enrichment_dir <- args[2]  # Output dir


#######################
# Positional distribution plots for RV
#########################
# Make density plot of distance between rare variants and splice sites 
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "10"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
#make_distance_density_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)

# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites)
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "10"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
#make_distance_density_plot_seperated_by_ss_type(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)

# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites) using odds ratios of real vs background
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "10"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_odds_ratio_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
make_distance_odds_ratio_density_plot_seperated_by_ss_type(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)

distance <- "10"
version <- "observed_splice_site"
pvalue_threshold <- "1e-06"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_odds_ratio_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
make_distance_odds_ratio_density_plot_seperated_by_ss_type(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)















































######################################
# Old (currrently retired) scripts
######################################


#######################
# Concensus site changes for RV
#########################
# Make density plot of distance between rare variants and splice sites 
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "10"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_root <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_")
# make_consensus_splice_site_plots(inlier_distance_file, outlier_distance_file, output_root)





#######################
# Positional distribution plots for sQTLs
#########################
# Make density plot of distance between rare variants and splice sites 
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "10"
version <- "observed_splice_site"
pvalue_threshold <- "1e-11"
tissue_name <- "Muscle_Skeletal"
outlier_distance_file <- paste0(variant_position_enrichment_dir, tissue_name, "_sqtl_significant_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, tissue_name, "_sqtl_background_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, tissue_name, "_sqtl_distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
#make_distance_density_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)

# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites)
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "10"
version <- "observed_splice_site"
pvalue_threshold <- "1e-11"
tissue_name <- "Muscle_Skeletal"
outlier_distance_file <- paste0(variant_position_enrichment_dir, tissue_name, "_sqtl_significant_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, tissue_name, "_sqtl_background_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, tissue_name, "_sqtl_distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
#make_distance_density_plot_seperated_by_ss_type(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)



