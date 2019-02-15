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

make_novel_annotated_seperated_variant_allele_odds_ratio_plot <- function(distance, inlier_distance_file, outlier_distance_file, output_file, ss_type, position,title) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	print(paste0(ss_type, " ", position))

	position_outlier_distances_novel <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position & as.character(outlier_distances$annotated_splice_site)=="novel",]
	position_inlier_distances_novel <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance == position & as.character(inlier_distances$annotated_splice_site)=="novel",]

	alleles <- c("A", "C", "T", "G")

	variant_types <- c()
	outlier_status <- c()
	counts <- c()
	# Loop through all possible variants
	for (iter1 in 1:length(alleles)) {
		for (iter2 in 1:length(alleles)) {
			major_allele = alleles[iter1]
			variant_allele = alleles[iter2]
			if (major_allele != variant_allele) {
				variant_type <- paste0(major_allele, "->",variant_allele)
				a <- sum(position_outlier_distances_novel$major_allele == major_allele & position_outlier_distances_novel$variant_allele == variant_allele) 
				c <- sum(position_inlier_distances_novel$major_allele == major_allele & position_inlier_distances_novel$variant_allele == variant_allele) 

				# Add outliers
				variant_types <- c(variant_types, variant_type)
				outlier_status <- c(outlier_status, "outlier")
				counts <- c(counts, a)
				# Add inliers
				variant_types <- c(variant_types, variant_type)
				outlier_status <- c(outlier_status, "inlier")
				counts <- c(counts, c)

			}

		}
	}

	df_novel <- data.frame(variant_types=factor(variant_types), outlier_status = factor(outlier_status), counts=counts)
	
	position_outlier_distances_annotated <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position & as.character(outlier_distances$annotated_splice_site)=="annotated",]
	position_inlier_distances_annotated <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance == position & as.character(inlier_distances$annotated_splice_site)=="annotated",]

	alleles <- c("A", "C", "T", "G")

	variant_types <- c()
	outlier_status <- c()
	counts <- c()
	# Loop through all possible variants
	for (iter1 in 1:length(alleles)) {
		for (iter2 in 1:length(alleles)) {
			major_allele = alleles[iter1]
			variant_allele = alleles[iter2]
			if (major_allele != variant_allele) {
				variant_type <- paste0(major_allele, "->",variant_allele)
				a <- sum(position_outlier_distances_annotated$major_allele == major_allele & position_outlier_distances_annotated$variant_allele == variant_allele) 
				c <- sum(position_inlier_distances_annotated$major_allele == major_allele & position_inlier_distances_annotated$variant_allele == variant_allele) 

				# Add outliers
				variant_types <- c(variant_types, variant_type)
				outlier_status <- c(outlier_status, "outlier")
				counts <- c(counts, a)
				# Add inliers
				variant_types <- c(variant_types, variant_type)
				outlier_status <- c(outlier_status, "inlier")
				counts <- c(counts, c)

			}

		}
	}

	df_annotated <- data.frame(variant_types=factor(variant_types), outlier_status = factor(outlier_status), counts=counts)



	options(bitmapType = 'cairo', device = 'pdf')

	outlier_novel_counts <- dim(position_outlier_distances_novel)[1] 
	inlier_novel_counts <- dim(position_inlier_distances_novel)[1]

	outlier_annotated_counts <- dim(position_outlier_distances_annotated)[1] 
	inlier_annotated_counts <- dim(position_inlier_distances_annotated)[1]

	plotter_novel <- ggplot(df_novel,aes(x=outlier_status, y=counts, fill=variant_types)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="",y="proportion",fill="", title=paste0("Novel ", title,"\n")) + 
    	scale_fill_manual(values=c("dodgerblue4", "dodgerblue3", "deepskyblue2", "palegreen4", "palegreen3", "palegreen1","violetred4", "violetred3","violetred1", "goldenrod4", "goldenrod3", "goldenrod1")) + 
    	theme(text = element_text(size=12),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=12), legend.title = element_text(size=11)) 

   legend <- get_legend(plotter_novel + theme(legend.position="bottom"))


    plotter_novel <- ggdraw(plotter_novel + theme(legend.position="none")) + draw_label(as.character(inlier_novel_counts), x=.37,y=.81) + draw_label(as.character(outlier_novel_counts), x=.76,y=.81)

    	#y,x
    #plotter <- ggdraw(plotter) + draw_label(as.character(inlier_counts), x=.31,y=.85) + draw_label(as.character(outlier_counts), x=.61,y=.85)

    plotter_annotated <- ggplot(df_annotated,aes(x=outlier_status, y=counts, fill=variant_types)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="",y="proportion",fill="", title=paste0("Annotated ", title,"\n")) + 
    	scale_fill_manual(values=c("dodgerblue4", "dodgerblue3", "deepskyblue2", "palegreen4", "palegreen3", "palegreen1","violetred4", "violetred3","violetred1", "goldenrod4", "goldenrod3", "goldenrod1")) + 
    	theme(text = element_text(size=12),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=12), legend.title = element_text(size=11)) 


    plotter_annotated <- ggdraw(plotter_annotated + theme(legend.position="none")) + draw_label(as.character(inlier_annotated_counts), x=.37,y=.81) + draw_label(as.character(outlier_annotated_counts), x=.76,y=.81)


	plotter <- plot_grid(plotter_novel, plotter_annotated, legend, ncol=1, rel_heights=c(1,1,.25))
	ggsave(plotter, file=output_file,width = 16,height=16,units="cm")







}

make_variant_allele_odds_ratio_plot <- function(distance, inlier_distance_file, outlier_distance_file, output_file, ss_type, position,title) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	print(paste0(ss_type, " ", position))

	position_outlier_distances <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position,]
	position_inlier_distances <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance == position,]

	alleles <- c("A", "C", "T", "G")

	variant_types <- c()
	outlier_status <- c()
	counts <- c()
	# Loop through all possible variants
	for (iter1 in 1:length(alleles)) {
		for (iter2 in 1:length(alleles)) {
			major_allele = alleles[iter1]
			variant_allele = alleles[iter2]
			if (major_allele != variant_allele) {
				variant_type <- paste0(major_allele, "->",variant_allele)
				a <- sum(position_outlier_distances$major_allele == major_allele & position_outlier_distances$variant_allele == variant_allele) 
				c <- sum(position_inlier_distances$major_allele == major_allele & position_inlier_distances$variant_allele == variant_allele) 

				# Add outliers
				variant_types <- c(variant_types, variant_type)
				outlier_status <- c(outlier_status, "outlier")
				counts <- c(counts, a)
				# Add inliers
				variant_types <- c(variant_types, variant_type)
				outlier_status <- c(outlier_status, "inlier")
				counts <- c(counts, c)

			}

		}
	}

	df <- data.frame(variant_types=factor(variant_types), outlier_status = factor(outlier_status), counts=counts)
	options(bitmapType = 'cairo', device = 'pdf')

	outlier_counts <- dim(position_outlier_distances)[1] 
	inlier_counts <- dim(position_inlier_distances)[1]
	

	plotter <- ggplot(df,aes(x=outlier_status, y=counts, fill=variant_types)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="",y="proportion",fill="", title=paste0(title,"\n")) + 
    	scale_fill_manual(values=c("dodgerblue4", "dodgerblue3", "deepskyblue2", "palegreen4", "palegreen3", "palegreen1","violetred4", "violetred3","violetred1", "goldenrod4", "goldenrod3", "goldenrod1")) + 
    	theme(text = element_text(size=12),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=12), legend.title = element_text(size=11)) 

    	#y,x
    plotter <- ggdraw(plotter) + draw_label(as.character(inlier_counts), x=.31,y=.85) + draw_label(as.character(outlier_counts), x=.61,y=.85)

	ggsave(plotter, file=output_file,width = 15,height=10,units="cm")

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
					theme(axis.text.x=element_text(angle=45,hjust=1), text = element_text(size=14),axis.text=element_text(size=13), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=13), legend.title = element_text(size=14)) +
					scale_x_continuous(breaks=-10:9, labels=c("A-10","A-9","A-8","A-7","A-6","A-5","A-4","A-3","A-2","A-1","A+1","A+2","A+3","A+4","A+5","A+6", "A+7","A+8","A+9","A+10"))


	donor_plot <-  ggplot() + geom_errorbar(data=df_donor, mapping=aes(x=dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="darkorchid") +
					geom_point(data=df_donor, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
					labs(x = "Distance from splice site (BP)", y = "Enrichment", title="Donor splice site") +
					geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
					geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					theme(axis.text.x=element_text(angle=45,hjust=1),text = element_text(size=14),axis.text=element_text(size=13), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=13), legend.title = element_text(size=14)) +
					scale_x_continuous(breaks=-10:9, labels=c("D+10","D+9","D+8","D+7","D+6","D+5","D+4","D+3","D+2","D+1","D-1","D-2","D-3","D-4","D-5","D-6", "D-7","D-8","D-9","D-10"))


	error_bar_plot <- plot_grid(donor_plot,acceptor_plot, ncol=1)



	options(bitmapType = 'cairo', device = 'pdf')

	
	ggsave(error_bar_plot, file=output_file,width = 15,height=14,units="cm")
}

# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites) using odds ratios of real vs background
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
make_distance_odds_ratio_density_plot_seperated_by_ss_type_and_annotation_vs_novel <- function(distance_window, inlier_distance_file, outlier_distance_file, title, output_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	outlier_donor_distances_novel <- outlier_distances[as.character(outlier_distances$splice_site_type)=="donor" & as.character(outlier_distances$annotated_splice_site)=="novel",]
	outlier_acceptor_distances_novel <- outlier_distances[as.character(outlier_distances$splice_site_type)=="acceptor" & as.character(outlier_distances$annotated_splice_site)=="novel",]
	inlier_donor_distances_novel <- inlier_distances[as.character(inlier_distances$splice_site_type)=="donor" & as.character(inlier_distances$annotated_splice_site)=="novel",]
	inlier_acceptor_distances_novel <- inlier_distances[as.character(inlier_distances$splice_site_type)=="acceptor" & as.character(inlier_distances$annotated_splice_site)=="novel",]

	outlier_donor_distances_annotated <- outlier_distances[as.character(outlier_distances$splice_site_type)=="donor" & as.character(outlier_distances$annotated_splice_site)=="annotated",]
	outlier_acceptor_distances_annotated <- outlier_distances[as.character(outlier_distances$splice_site_type)=="acceptor" & as.character(outlier_distances$annotated_splice_site)=="annotated",]
	inlier_donor_distances_annotated <- inlier_distances[as.character(inlier_distances$splice_site_type)=="donor" & as.character(inlier_distances$annotated_splice_site)=="annotated",]
	inlier_acceptor_distances_annotated <- inlier_distances[as.character(inlier_distances$splice_site_type)=="acceptor" & as.character(inlier_distances$annotated_splice_site)=="annotated",]

	donor_df_novel <- extract_odds_ratio_data(outlier_donor_distances_novel, inlier_donor_distances_novel, 10)
	acceptor_df_novel <- extract_odds_ratio_data(outlier_acceptor_distances_novel, inlier_acceptor_distances_novel, 10)

	donor_df_annotated <- extract_odds_ratio_data(outlier_donor_distances_annotated, inlier_donor_distances_annotated, 10)
	acceptor_df_annotated <- extract_odds_ratio_data(outlier_acceptor_distances_annotated, inlier_acceptor_distances_annotated, 10)

	
	options(bitmapType = 'cairo', device = 'pdf')

	acceptor_novel_plot <-  ggplot() + geom_errorbar(data=acceptor_df_novel, mapping=aes(x=dist_to_ss,ymin=lower_bounds, ymax=upper_bounds),color="darkorchid") +
					geom_point(data=acceptor_df_novel, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
					labs(x = "Distance from splice site (BP)", y = "Enrichment", title="Novel acceptor splice site") +
					geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
					geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					theme(axis.text.x=element_text(angle=45,hjust=1), text = element_text(size=14),axis.text=element_text(size=13), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=13), legend.title = element_text(size=14)) +
					scale_x_continuous(breaks=-10:9, labels=c("A-10","A-9","A-8","A-7","A-6","A-5","A-4","A-3","A-2","A-1","A+1","A+2","A+3","A+4","A+5","A+6", "A+7","A+8","A+9","A+10"))


	donor_novel_plot <-  ggplot() + geom_errorbar(data=donor_df_novel, mapping=aes(x=dist_to_ss,ymin=lower_bounds, ymax=upper_bounds),color="darkorchid") +
					geom_point(data=donor_df_novel, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
					labs(x = "Distance from splice site (BP)", y = "Enrichment", title="Novel donor splice site") +
					geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
					geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					theme(axis.text.x=element_text(angle=45,hjust=1),text = element_text(size=14),axis.text=element_text(size=13), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=13), legend.title = element_text(size=14)) +
					scale_x_continuous(breaks=-10:9, labels=c("D+10","D+9","D+8","D+7","D+6","D+5","D+4","D+3","D+2","D+1","D-1","D-2","D-3","D-4","D-5","D-6", "D-7","D-8","D-9","D-10"))



	acceptor_annotated_plot <-  ggplot() + geom_errorbar(data=acceptor_df_annotated, mapping=aes(x=dist_to_ss,ymin=lower_bounds, ymax=upper_bounds),color="darkorchid") +
					geom_point(data=acceptor_df_annotated, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
					labs(x = "Distance from splice site (BP)", y = "Enrichment", title="Annotated acceptor splice site") +
					geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
					geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					theme(axis.text.x=element_text(angle=45,hjust=1), text = element_text(size=14),axis.text=element_text(size=13), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=13), legend.title = element_text(size=14)) +
					scale_x_continuous(breaks=-10:9, labels=c("A-10","A-9","A-8","A-7","A-6","A-5","A-4","A-3","A-2","A-1","A+1","A+2","A+3","A+4","A+5","A+6", "A+7","A+8","A+9","A+10"))


	donor_annotated_plot <-  ggplot() + geom_errorbar(data=donor_df_annotated, mapping=aes(x=dist_to_ss,ymin=lower_bounds, ymax=upper_bounds),color="darkorchid") +
					geom_point(data=donor_df_annotated, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
					labs(x = "Distance from splice site (BP)", y = "Enrichment", title="Annotated donor splice site") +
					geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
					geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					theme(axis.text.x=element_text(angle=45,hjust=1),text = element_text(size=14),axis.text=element_text(size=13), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=13), legend.title = element_text(size=14)) +
					scale_x_continuous(breaks=-10:9, labels=c("D+10","D+9","D+8","D+7","D+6","D+5","D+4","D+3","D+2","D+1","D-1","D-2","D-3","D-4","D-5","D-6", "D-7","D-8","D-9","D-10"))



	error_bar_plot <- plot_grid(donor_novel_plot, acceptor_novel_plot, donor_annotated_plot, acceptor_annotated_plot, ncol=2)




	
	ggsave(error_bar_plot, file=output_file,width = 26,height=16,units="cm")
}












#######################
# Command Line args
#######################
variant_position_enrichment_dir <- args[1]  # Input dir
visualize_variant_position_enrichment_dir <- args[2]  # Output dir






distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_novel_annotated_seperated_d+1_plot.pdf")
ss_type <- "donor"
position <- -1
make_novel_annotated_seperated_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+1 (Consensus=G)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_novel_annotated_seperated_d+2_plot.pdf")
ss_type <- "donor"
position <- -2
make_novel_annotated_seperated_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+2 (Consensus=T)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_novel_annotated_seperated_a-1_plot.pdf")
ss_type <- "acceptor"
position <- -1
make_novel_annotated_seperated_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "A-1 (Consensus=G)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_novel_annotated_seperated_a-2_plot.pdf")
ss_type <- "acceptor"
position <- -2
make_novel_annotated_seperated_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "A-2 (Consensus=A)")


output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_novel_annotated_seperated_a+1_plot.pdf")
ss_type <- "acceptor"
position <- 0
make_novel_annotated_seperated_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "A+1 (Consensus=G)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_novel_annotated_seperated_d-1_plot.pdf")
ss_type <- "donor"
position <- 0
make_novel_annotated_seperated_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D-1 (Consensus=G)")


output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_novel_annotated_seperated_d+3_plot.pdf")
ss_type <- "donor"
position <- -3
make_novel_annotated_seperated_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+3 (Consensus=A)")


output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_novel_annotated_seperated_d+4_plot.pdf")
ss_type <- "donor"
position <- -4
make_novel_annotated_seperated_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+4 (Consensus=A)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_novel_annotated_seperated_d+5_plot.pdf")
ss_type <- "donor"
position <- -5
make_novel_annotated_seperated_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+5 (Consensus=G)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_novel_annotated_seperated_d+6_plot.pdf")
ss_type <- "donor"
position <- -6
make_novel_annotated_seperated_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+6 (Consensus=T)")



#######################
# Specific allele plots
#########################

distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_d+1_plot.pdf")
ss_type <- "donor"
position <- -1
make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+1 (Consensus=G)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_d+2_plot.pdf")
ss_type <- "donor"
position <- -2
make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+2 (Consensus=T)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_a-1_plot.pdf")
ss_type <- "acceptor"
position <- -1
make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "A-1 (Consensus=G)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_a-2_plot.pdf")
ss_type <- "acceptor"
position <- -2
make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "A-2 (Consensus=A)")


output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_a+1_plot.pdf")
ss_type <- "acceptor"
position <- 0
make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "A+1 (Consensus=G)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_d-1_plot.pdf")
ss_type <- "donor"
position <- 0
make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D-1 (Consensus=G)")


output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_d+3_plot.pdf")
ss_type <- "donor"
position <- -3
make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+3 (Consensus=A)")


output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_d+4_plot.pdf")
ss_type <- "donor"
position <- -4
make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+4 (Consensus=A)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_d+5_plot.pdf")
ss_type <- "donor"
position <- -5
make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+5 (Consensus=G)")

output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_alleles_proportion_d+6_plot.pdf")
ss_type <- "donor"
position <- -6
make_variant_allele_odds_ratio_plot(as.numeric(distance), inlier_distance_file, outlier_distance_file, output_file, ss_type, position, "D+6 (Consensus=T)")

#######################
# Positional distribution plots for RV seperated by donor vs acceptor splice sites
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

#######################
# Positional distribution plots for RV sepertated by donor vs acceptor splice sites and annotated vs unannotated splice sites
#########################
# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites) using odds ratios of real vs background
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_annotated_vs_novel_odds_ratio_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
make_distance_odds_ratio_density_plot_seperated_by_ss_type_and_annotation_vs_novel(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)

distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-06"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_annotated_vs_novel_odds_ratio_density_plot.pdf")
title <- paste0("pvalue=", pvalue_threshold, " / version=", version)
make_distance_odds_ratio_density_plot_seperated_by_ss_type_and_annotation_vs_novel(as.numeric(distance), inlier_distance_file, outlier_distance_file, title, output_file)





