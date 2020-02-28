args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)
library(png)
library(ggseqlogo)

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

mutation_type_bar_plot <- function(inlier_distances, outlier_distances, ss_type, position, title, y_axis) {
	position_outlier_distances <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position,]
	position_inlier_distances <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance == position,]
	alleles <- c("A", "C", "T", "G")
	# Initialize output arrays
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
				outlier_status <- c(outlier_status, "Outlier")
				counts <- c(counts, a)
				# Add inliers
				variant_types <- c(variant_types, variant_type)
				outlier_status <- c(outlier_status, "Inlier")
				counts <- c(counts, c)
			}
		}
	}

	# Organize into compact data frame
	df <- data.frame(variant_types=factor(variant_types), outlier_status = factor(outlier_status), counts=counts)
	# Count number of outliers
	outlier_counts <- dim(position_outlier_distances)[1] 
	# Count number of inliers
	inlier_counts <- dim(position_inlier_distances)[1]

	plotter <- ggplot(df,aes(x=outlier_status, y=counts, fill=variant_types)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="",y="Proportion",fill="", title="") + 
    	scale_fill_manual(values=c("dodgerblue4", "dodgerblue3", "deepskyblue2", "palegreen4", "palegreen3", "palegreen1","violetred4", "violetred3","violetred1", "goldenrod4", "goldenrod3", "goldenrod1")) + 
    	gtex_v8_figure_theme()
    	#y,x
    x_pos_1 = .61
    x_pos_2 = .84
    x_pos_3 = .74
    if (y_axis == FALSE) {
    	plotter <- plotter + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank())
        x_pos_1 = .32
        x_pos_2 = .76
        x_pos_3 = .55
    }
    plotter <- ggdraw(plotter + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(legend.position="none")) +
     		   draw_label(as.character(inlier_counts), x=x_pos_1,y=.87,size=6) + 
     		   draw_label(as.character(outlier_counts), x=x_pos_2,y=.87,size=6) + 
     		   draw_label(title, x=x_pos_3, y=.12, size=9)
    return(plotter)
}

get_mutation_bar_plot_legend <- function(inlier_distances, outlier_distances, ss_type, position, title, y_axis) {
	position_outlier_distances <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position,]
	position_inlier_distances <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance == position,]
	alleles <- c("A", "C", "T", "G")
	# Initialize output arrays
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
				outlier_status <- c(outlier_status, "Outlier")
				counts <- c(counts, a)
				# Add inliers
				variant_types <- c(variant_types, variant_type)
				outlier_status <- c(outlier_status, "Inlier")
				counts <- c(counts, c)
			}
		}
	}

	# Organize into compact data frame
	df <- data.frame(variant_types=factor(variant_types), outlier_status = factor(outlier_status), counts=counts)
	# Count number of outliers
	outlier_counts <- dim(position_outlier_distances)[1] 
	# Count number of inliers
	inlier_counts <- dim(position_inlier_distances)[1]

	plotter <- ggplot(df,aes(x=outlier_status, y=counts, fill=variant_types)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="",y="Proportion",fill="", title=paste0(title,"\n")) + 
    	scale_fill_manual(values=c("dodgerblue4", "dodgerblue3", "deepskyblue2", "palegreen4", "palegreen3", "palegreen1","violetred4", "violetred3","violetred1", "goldenrod4", "goldenrod3", "goldenrod1")) + 
    	gtex_v8_figure_theme()
    	#y,x
 	legend <- get_legend(plotter + theme(plot.title=element_text(margin=margin(0,0,0,0)), legend.text = element_text(size=6),legend.key.size = unit(0.1, "cm")))
    return(legend)
}



gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


mutation_type_bar_plot_across_concensus_sites <- function(inlier_distance_file, outlier_distance_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	d_6_plot <- mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	d_5_plot <- mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	d_4_plot <- mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	d_3_plot <- mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -2
	ss_name <- "D+2"
	d_2_plot <- mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -1
	ss_name <- "D+1"
	d_1_plot <- mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	d_minus_1_plot <- mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, TRUE)

	ss_type <- "acceptor"
	position <- -2
	ss_name <- "A-2"
	a_minus_2_plot <- mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, FALSE)

	ss_type <- "acceptor"
	position <- -1
	ss_name <- "A-1"
	a_minus_1_plot <- mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, FALSE)

	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	a_1_plot <- mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, FALSE)

	mutation_bar_plot_legend <- get_mutation_bar_plot_legend(inlier_distances, outlier_distances, "acceptor", 0, "A+1", TRUE)

	gg <- plot_grid(d_minus_1_plot, d_1_plot, d_2_plot, d_3_plot, d_4_plot, d_5_plot, d_6_plot, a_minus_2_plot, a_minus_1_plot, a_1_plot,mutation_bar_plot_legend, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1, 1, 1, 1, 1, .58))

	return(gg)

}

mutation_type_pwm_plot_across_annotated_concensus_sites_two_row <- function(input_file) {
	outlier_distances <- read.table(input_file, header=TRUE)

	mutations <- as.character(outlier_distances$variant_allele)

	refs <- c()
	alts <- c()
	for (iter in 1:length(mutations)) {
		mutation <- mutations[iter]
		ref <- strsplit(mutation, "->")[[1]][1]
		alt <- strsplit(mutation, "->")[[1]][2]
		refs <- c(refs, ref)
		alts <- c(alts, alt)
	}

	alleles <- c("A", "C", "T", "G")

	ss_names <- c("D-1","D+1", "D+2", "D+3", "D+4", "D+5", "D+6")

	#####################################################
	# Outlier variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		variant_alleles <- alts[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(variant_alleles==allele)
		}
	}
	
	outlier_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier variants at donor sites") + 
	        labs(x="",y="Frequency",title="") + 
	        gtex_v8_figure_theme() +
	       	theme(plot.title = element_text(hjust = 0.5)) +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names)


	#####################################################
	# Background variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		outlier_background_alleles <- refs[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(outlier_background_alleles==allele)
		}
	}

	background_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier background at donor sites") + 
	        labs(x="",y="Frequency",title="") + 
	        theme(plot.title = element_text(hjust = 0.5)) +
	        gtex_v8_figure_theme() +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names)
	

	combined_donor_plotter <- plot_grid(background_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),nrow=2)


	ss_names <- c("A-2","A-1", "A+1")


	#####################################################
	# Outlier variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		variant_alleles <- alts[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(variant_alleles==allele)
		}
	}
	
	outlier_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier variants at acceptor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        gtex_v8_figure_theme() +
	       	theme(plot.title = element_text(hjust = 0.5)) +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names) +
	        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank())


	#####################################################
	# Background variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		outlier_background_alleles <- refs[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]

		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(outlier_background_alleles==allele)
		}
	}

	background_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier background at acceptor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        theme(plot.title = element_text(hjust = 0.5)) +
	        gtex_v8_figure_theme() +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names) +
	        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank())

	

	combined_acceptor_plotter <- plot_grid(background_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),nrow=2)

	combined <- plot_grid(combined_donor_plotter, combined_acceptor_plotter, ncol=2, rel_widths=c(7,3)) +
				draw_text("Background for splice sites with low junction usage", x=.55, y=.97, size=8)+
				draw_text("sOutlier variants for splice sites with low junction usage", x=.55, y=.47, size=8)

	return(combined)

}


mutation_type_pwm_plot_across_annotated_concensus_sites <- function(input_file) {
	outlier_distances <- read.table(input_file, header=TRUE)

	mutations <- as.character(outlier_distances$variant_allele)

	refs <- c()
	alts <- c()
	for (iter in 1:length(mutations)) {
		mutation <- mutations[iter]
		ref <- strsplit(mutation, "->")[[1]][1]
		alt <- strsplit(mutation, "->")[[1]][2]
		refs <- c(refs, ref)
		alts <- c(alts, alt)
	}

	alleles <- c("A", "C", "T", "G")

	ss_names <- c("D-1","D+1", "D+2", "D+3", "D+4", "D+5", "D+6")

	#####################################################
	# Outlier variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		variant_alleles <- alts[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(variant_alleles==allele)
		}
	}
	
	outlier_plotter_donor <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier variants at donor sites") + 
	        labs(x="",y="Frequency",title="") + 
	        gtex_v8_figure_theme() +
	       	theme(plot.title = element_text(hjust = 0.5)) +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names)


	#####################################################
	# Background variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		outlier_background_alleles <- refs[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(outlier_background_alleles==allele)
		}
	}

	background_plotter_donor <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier background at donor sites") + 
	        labs(x="",y="Frequency",title="") + 
	        theme(plot.title = element_text(hjust = 0.5)) +
	        gtex_v8_figure_theme() +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names)
	

	#combined_donor_plotter <- plot_grid(background_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),nrow=2)


	ss_names <- c("A-2","A-1", "A+1")


	#####################################################
	# Outlier variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		variant_alleles <- alts[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(variant_alleles==allele)
		}
	}
	
	outlier_plotter_acceptor <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier variants at acceptor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        gtex_v8_figure_theme() +
	       	theme(plot.title = element_text(hjust = 0.5)) +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names) +
	        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank())


	#####################################################
	# Background variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		outlier_background_alleles <- refs[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]

		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(outlier_background_alleles==allele)
		}
	}

	background_plotter_acceptor <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier background at acceptor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        theme(plot.title = element_text(hjust = 0.5)) +
	        gtex_v8_figure_theme() +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names) +
	        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank())

	

	#combined_acceptor_plotter <- plot_grid(background_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),nrow=2)

	#combined <- plot_grid(combined_donor_plotter, combined_acceptor_plotter, ncol=2, rel_widths=c(7,3)) +
	#			draw_text("Background for splice sites with low junction usage", x=.55, y=.97, size=8)+
	#			draw_text("Outlier variants for splice sites with low junction usage", x=.55, y=.47, size=8)
	combined <- plot_grid(background_plotter_donor + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), background_plotter_acceptor + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_plotter_donor + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_plotter_acceptor + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol=4, rel_widths=c(8.1,3,8.1,3)) +
				draw_text("Background for splice sites with low junction usage", x=.28, y=.95, size=8) +
				draw_text("sOutlier variants for splice sites with low junction usage", x=.78, y=.95, size=8)

	return(combined)

}

mutation_type_pwm_plot_across_novel_concensus_sites <- function(input_file) {
	outlier_distances <- read.table(input_file, header=TRUE)

	mutations <- as.character(outlier_distances$variant_allele)

	refs <- c()
	alts <- c()
	for (iter in 1:length(mutations)) {
		mutation <- mutations[iter]
		ref <- strsplit(mutation, "->")[[1]][1]
		alt <- strsplit(mutation, "->")[[1]][2]
		refs <- c(refs, ref)
		alts <- c(alts, alt)
	}

	alleles <- c("A", "C", "T", "G")

	ss_names <- c("D-1","D+1", "D+2", "D+3", "D+4", "D+5", "D+6")

	#####################################################
	# Outlier variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		variant_alleles <- alts[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio > 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(variant_alleles==allele)
		}
	}
	
	outlier_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier variants at donor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        gtex_v8_figure_theme() +
	       	theme(plot.title = element_text(hjust = 0.5)) +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names)


	#####################################################
	# Background variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		outlier_background_alleles <- refs[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio > 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(outlier_background_alleles==allele)
		}
	}

	background_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Background at donor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        theme(plot.title = element_text(hjust = 0.5)) +
	        gtex_v8_figure_theme() +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names)
	

	combined_donor_plotter <- plot_grid(background_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),nrow=2)


	ss_names <- c("A-2","A-1", "A+1")


	#####################################################
	# Outlier variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		variant_alleles <- alts[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio > 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(variant_alleles==allele)
		}
	}
	
	outlier_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier variants at acceptor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        gtex_v8_figure_theme() +
	       	theme(plot.title = element_text(hjust = 0.5)) +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names) +
	        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank())


	#####################################################
	# Background variant PWM plot
	#####################################################
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		outlier_background_alleles <- refs[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio > 1.0]

		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(outlier_background_alleles==allele)
		}
	}

	background_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        labs(x="",y="Frequency", title="") + 
	        #labs(x="",y="Frequency", title="Background at acceptor sites") + 
	        theme(plot.title = element_text(hjust = 0.5)) +
	        gtex_v8_figure_theme() +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names) +
	        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank())

	

	combined_acceptor_plotter <- plot_grid(background_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),nrow=2)

	combined <- plot_grid(combined_donor_plotter, combined_acceptor_plotter, ncol=2, rel_widths=c(7,3)) +
				draw_text("Background for splice sites with high junction usage", x=.55, y=.97, size=8)+
				draw_text("Outlier variants for splice sites with high junction usage", x=.55, y=.47, size=8)

	return(combined)
}


mutation_type_bar_plot_across_concensus_sites_seperated_by_novel_and_rare <- function(inlier_distance_file, outlier_distance_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	inlier_distances_novel <- inlier_distances[as.character(inlier_distances$annotated_splice_site) == "novel", ]
	outlier_distances_novel <- outlier_distances[as.character(outlier_distances$annotated_splice_site) == "novel", ]

	inlier_distances_annotated <- inlier_distances[as.character(inlier_distances$annotated_splice_site) == "annotated", ]
	outlier_distances_annotated <- outlier_distances[as.character(outlier_distances$annotated_splice_site) == "annotated", ]


	#################
	# Annotated
	#################
	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	d_6_plot_ano <- mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	d_5_plot_ano <- mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	d_4_plot_ano <- mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	d_3_plot_ano <- mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -2
	ss_name <- "D+2"
	d_2_plot_ano <- mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -1
	ss_name <- "D+1"
	d_1_plot_ano <- mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	d_minus_1_plot_ano <- mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, TRUE)

	ss_type <- "acceptor"
	position <- -2
	ss_name <- "A-2"
	a_minus_2_plot_ano <- mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, FALSE)

	ss_type <- "acceptor"
	position <- -1
	ss_name <- "A-1"
	a_minus_1_plot_ano <- mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, FALSE)

	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	a_1_plot_ano <- mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, FALSE)

	#################
	# NOVEL
	#################
	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	d_6_plot <- mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	d_5_plot <- mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	d_4_plot <- mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	d_3_plot <- mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -2
	ss_name <- "D+2"
	d_2_plot <- mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- -1
	ss_name <- "D+1"
	d_1_plot <- mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, FALSE)

	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	d_minus_1_plot <- mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, TRUE)

	ss_type <- "acceptor"
	position <- -2
	ss_name <- "A-2"
	a_minus_2_plot <- mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, FALSE)

	ss_type <- "acceptor"
	position <- -1
	ss_name <- "A-1"
	a_minus_1_plot <- mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, FALSE)

	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	a_1_plot <- mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, FALSE)

	mutation_bar_plot_legend <- get_mutation_bar_plot_legend(inlier_distances_novel, outlier_distances_novel, "acceptor", 0, "A+1", TRUE)


	gg_ano <- plot_grid(d_minus_1_plot_ano, d_1_plot_ano, d_2_plot_ano, d_3_plot_ano, d_4_plot_ano, d_5_plot_ano, d_6_plot_ano, a_minus_2_plot_ano, a_minus_1_plot_ano, a_1_plot_ano, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1, 1, 1, 1, 1))

	gg_novel <- plot_grid(d_minus_1_plot, d_1_plot, d_2_plot, d_3_plot, d_4_plot, d_5_plot, d_6_plot, a_minus_2_plot, a_minus_1_plot, a_1_plot, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1, 1, 1, 1, 1))

	combined_plot <- plot_grid(plot_grid(gg_ano, gg_novel, ncol=1), mutation_bar_plot_legend,ncol=2,rel_widths=c(1,.05)) +
					draw_label("Annotated Splice Site", x=.5,y=.98,size=9) + 
					draw_label("Novel Splice Site", x=.5,y=.49,size=9) 

	return(combined_plot)

}


broad_mutation_type_bar_plot <- function(inlier_distances, outlier_distances, ss_type, position, title,concensus_allele, y_axis) {
	position_outlier_distances <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position,]
	position_inlier_distances <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance == position,]
	# Initialize output arrays
	variant_types <- c()
	outlier_status <- c()
	counts <- c()

	###############
	# Variant changes to concensus allele
	###############
	variant_type <- "Concencus Created"
	a <- sum(position_outlier_distances$variant_allele == concensus_allele) 
	c <- sum(position_inlier_distances$variant_allele == concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)

	###############
	# Variant changes from concensus allele
	###############
	variant_type <- "Concencus Destroyed"
	a <- sum(position_outlier_distances$major_allele == concensus_allele) 
	c <- sum(position_inlier_distances$major_allele == concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)

	###############
	# Variant changes neither to or from concensus allele
	###############
	variant_type <- "Neither"
	a <- sum(position_outlier_distances$major_allele != concensus_allele & position_outlier_distances$variant_allele != concensus_allele) 
	c <- sum(position_inlier_distances$major_allele != concensus_allele & position_inlier_distances$variant_allele != concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)



	
	# Organize into compact data frame
	df <- data.frame(variant_types=factor(variant_types,levels=c("Neither", "Concencus Destroyed", "Concencus Created")), outlier_status = factor(outlier_status), counts=counts)
	# Count number of outliers
	outlier_counts <- dim(position_outlier_distances)[1] 
	# Count number of inliers
	inlier_counts <- dim(position_inlier_distances)[1]

	plotter <- ggplot(df,aes(x=outlier_status, y=counts, fill=variant_types)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="",y="Proportion",fill="", title="") + 
    	scale_fill_manual(values=c("grey", "palevioletred3","skyblue3")) + 
    	gtex_v8_figure_theme()
    	#y,x
    x_pos_1 = .61
    x_pos_2 = .84
    x_pos_3 = .75
    if (y_axis == FALSE) {
    	plotter <- plotter + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank())
        x_pos_1 = .32
        x_pos_2 = .76
        x_pos_3 = .55
    }
    plotter <- ggdraw(plotter + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(legend.position="none")) +
     		   draw_label(title, x=x_pos_3, y=.12, size=8)
    return(plotter)
}

broad_mutation_type_bar_plot_outlier_only <- function(outlier_distances, ss_type, position, title,concensus_allele, y_axis) {
	position_outlier_distances <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position,]
	# Initialize output arrays
	variant_types <- c()
	outlier_status <- c()
	counts <- c()

	###############
	# Variant changes to concensus allele
	###############
	variant_type <- "to_concensus"
	a <- sum(position_outlier_distances$variant_allele == concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, title)
	counts <- c(counts, a)

	###############
	# Variant changes from concensus allele
	###############
	variant_type <- "from_concensus"
	a <- sum(position_outlier_distances$major_allele == concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, title)
	counts <- c(counts, a)

	###############
	# Variant changes neither to or from concensus allele
	###############
	variant_type <- "neither"
	a <- sum(position_outlier_distances$major_allele != concensus_allele & position_outlier_distances$variant_allele != concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, title)
	counts <- c(counts, a)

	
	# Organize into compact data frame
	df <- data.frame(variant_types=factor(variant_types,levels=c("from_concensus", "to_concensus","neither")), outlier_status = factor(outlier_status), counts=counts)
	# Count number of outliers
	#outlier_counts <- dim(position_outlier_distances)[1] 
	# Count number of inliers
	#inlier_counts <- dim(position_inlier_distances)[1]

	plotter <- ggplot(df,aes(x=outlier_status, y=counts, fill=variant_types)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="",y="proportion",fill="", title="") + 
    	scale_fill_manual(values=c("palevioletred3","skyblue3", "grey")) + 
    	gtex_v8_figure_theme() +
    	theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    	#y,x
    if (y_axis == FALSE) {
    	plotter <- plotter + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank())
        x_pos_1 = .32
        x_pos_2 = .76
        x_pos_3 = .55
    }
    return(plotter)
}


get_broad_mutation_type_bar_plot_legend <- function(inlier_distances, outlier_distances, ss_type, position, title,concensus_allele, y_axis) {
	position_outlier_distances <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position,]
	position_inlier_distances <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance == position,]
	# Initialize output arrays
	variant_types <- c()
	outlier_status <- c()
	counts <- c()

	###############
	# Variant changes to concensus allele
	###############
	variant_type <- "Consensus Created"
	a <- sum(position_outlier_distances$variant_allele == concensus_allele) 
	c <- sum(position_inlier_distances$variant_allele == concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)

	###############
	# Variant changes from concensus allele
	###############
	variant_type <- "Consensus Destroyed"
	a <- sum(position_outlier_distances$major_allele == concensus_allele) 
	c <- sum(position_inlier_distances$major_allele == concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)

	###############
	# Variant changes neither to or from concensus allele
	###############
	variant_type <- "Neither"
	a <- sum(position_outlier_distances$major_allele != concensus_allele & position_outlier_distances$variant_allele != concensus_allele) 
	c <- sum(position_inlier_distances$major_allele != concensus_allele & position_inlier_distances$variant_allele != concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)



	
	# Organize into compact data frame
	df <- data.frame(variant_types=factor(variant_types,levels=c("Consensus Created", "Consensus Destroyed", "Neither")), outlier_status = factor(outlier_status), counts=counts)
	# Count number of outliers
	outlier_counts <- dim(position_outlier_distances)[1] 
	# Count number of inliers
	inlier_counts <- dim(position_inlier_distances)[1]
	plotter <- ggplot(df,aes(x=outlier_status, y=counts, fill=variant_types)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="",y="proportion",fill="", title="") + 
    	scale_fill_manual(values=c("skyblue3", "palevioletred3","grey")) + 
    	gtex_v8_figure_theme()
    	#y,x

    return(get_legend(plotter + theme(legend.text = element_text(size=7),legend.key.size = unit(0.1, "cm"), legend.position="bottom")))
}

broad_mutation_type_bar_plot_across_concensus_sites_seperated_by_novel_and_rare_v2 <- function(inlier_distance_file, outlier_distance_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	inlier_distances_novel <- inlier_distances[as.character(inlier_distances$annotated_splice_site) == "novel", ]
	outlier_distances_novel <- outlier_distances[as.character(outlier_distances$annotated_splice_site) == "novel", ]

	inlier_distances_annotated <- inlier_distances[as.character(inlier_distances$annotated_splice_site) == "annotated", ]
	outlier_distances_annotated <- outlier_distances[as.character(outlier_distances$annotated_splice_site) == "annotated", ]


	################
	# Annotated
	################
	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	concensus <- "T"
	d_6_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	concensus <- "G"
	d_5_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	concensus <- "A"
	d_4_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	concensus <- "A"
	d_3_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -2
	ss_name <- "D+2"
	concensus <- "T"
	d_2_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -1
	ss_name <- "D+1"
	concensus <- "G"
	d_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	concensus <- "G"
	d_minus_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, TRUE)

	ss_type <- "acceptor"
	position <- -2
	ss_name <- "A-2"
	concensus <- "A"
	a_minus_2_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- -1
	ss_name <- "A-1"
	concensus <- "G"
	a_minus_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	concensus <- "G"
	a_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	################
	# NOVEL
	################
	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	concensus <- "T"
	d_6_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	concensus <- "G"
	d_5_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	concensus <- "A"
	d_4_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	concensus <- "A"
	d_3_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -2
	ss_name <- "D+2"
	concensus <- "T"
	d_2_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -1
	ss_name <- "D+1"
	concensus <- "G"
	d_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	concensus <- "G"
	d_minus_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, TRUE)

	ss_type <- "acceptor"
	position <- -2
	ss_name <- "A-2"
	concensus <- "A"
	a_minus_2_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- -1
	ss_name <- "A-1"
	concensus <- "G"
	a_minus_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	concensus <- "G"
	a_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)




	mutation_bar_plot_legend <- get_broad_mutation_type_bar_plot_legend(inlier_distances_novel, outlier_distances_novel, "acceptor", 0, "A+1", "G", TRUE)

	gg_novel <- plot_grid(d_minus_1_plot, d_3_plot, d_4_plot, d_5_plot, d_6_plot, a_1_plot, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1))
	gg_ano <- plot_grid(d_minus_1_plot_ano, d_3_plot_ano, d_4_plot_ano, d_5_plot_ano, d_6_plot_ano, a_1_plot_ano, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1))

	combined_gg <- plot_grid(gg_ano, gg_novel, ncol=1)+
					draw_label("Annotated Splice Site", x=.5,y=.93,size=8) + 
					draw_label("Novel Splice Site", x=.5,y=.43,size=8) 

	#combined_gg <- ggdraw() + 
	#			draw_plot(gg_ano,0,.43,1,.55) +
	#			draw_plot(gg_novel,0,.02,1,.55) + 
	#			draw_plot(mutation_bar_plot_legend,.33,-.43,1,1) +
	#			draw_label("Annotated splice site", .53, .915, size=8) +
	#			draw_label("Novel splice site", .53, .51, size=8) 

	#return(combined_gg)
	return(plot_grid(combined_gg, mutation_bar_plot_legend, nrow=2, rel_heights=c(1,.05)))

}

make_panel_2d <- function(outlier_distance_file) {
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)
	outlier_distances_novel <- outlier_distances[as.character(outlier_distances$annotated_splice_site) == "novel", ]
	outlier_distances_annotated <- outlier_distances[as.character(outlier_distances$annotated_splice_site) == "annotated", ]

	# Make first row (annotated splice sites)

	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	concensus <- "T"
	d_6_anno_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_annotated, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	concensus <- "G"
	d_5_anno_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_annotated, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	concensus <- "A"
	d_4_anno_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_annotated, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	concensus <- "A"
	d_3_anno_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_annotated, ss_type, position, ss_name, concensus,  FALSE)


	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	concensus <- "G"
	d_minus_1_anno_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_annotated, ss_type, position, ss_name, concensus,  TRUE)


	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	concensus <- "G"
	a_1_anno_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_annotated, ss_type, position, ss_name, concensus,  FALSE)

	gg_anno <- plot_grid(d_minus_1_anno_plot + theme(legend.position="none"), d_3_anno_plot+ theme(legend.position="none"), d_4_anno_plot+ theme(legend.position="none"), d_5_anno_plot+ theme(legend.position="none"), d_6_anno_plot+ theme(legend.position="none"), a_1_anno_plot+ theme(legend.position="none"), nrow=1, rel_widths=c(1.38,1, 1, 1, 1, 1))


	# Make second row (novel splice sites)

	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	concensus <- "T"
	d_6_novel_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_novel, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	concensus <- "G"
	d_5_novel_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_novel, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	concensus <- "A"
	d_4_novel_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_novel, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	concensus <- "A"
	d_3_novel_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_novel, ss_type, position, ss_name, concensus,  FALSE)


	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	concensus <- "G"
	d_minus_1_novel_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_novel, ss_type, position, ss_name, concensus,  TRUE)


	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	concensus <- "G"
	a_1_novel_plot <- broad_mutation_type_bar_plot_outlier_only(outlier_distances_novel, ss_type, position, ss_name, concensus,  FALSE)

	legend <- get_legend(a_1_anno_plot + theme(legend.position="bottom"))

	gg_novel <- plot_grid(d_minus_1_novel_plot + theme(legend.position="none"), d_3_novel_plot+ theme(legend.position="none"), d_4_novel_plot+ theme(legend.position="none"), d_5_novel_plot+ theme(legend.position="none"), d_6_novel_plot+ theme(legend.position="none"), a_1_novel_plot+ theme(legend.position="none"), nrow=1, rel_widths=c(1.38,1, 1, 1, 1, 1))

	panel_2d <- ggdraw() + 
				draw_plot(gg_anno,0,.43,1,.55) +
				draw_plot(gg_novel,0,.02,1,.55) + 
				draw_plot(legend,.33,-.43,1,1) +
				draw_label("Annotated splice site", .53, .915, size=8) +
				draw_label("Novel splice site", .53, .51, size=8) 

	return(panel_2d)


	#return(plot_grid(gg, mutation_bar_plot_legend, nrow=2, rel_heights=c(1,.1)))
	#return(gg)

}

broad_mutation_type_bar_plot_across_concensus_sites_seperated_by_novel_and_rare_v1 <- function(inlier_distance_file, outlier_distance_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	inlier_distances_novel <- inlier_distances[as.character(inlier_distances$annotated_splice_site) == "novel", ]
	outlier_distances_novel <- outlier_distances[as.character(outlier_distances$annotated_splice_site) == "novel", ]

	inlier_distances_annotated <- inlier_distances[as.character(inlier_distances$annotated_splice_site) == "annotated", ]
	outlier_distances_annotated <- outlier_distances[as.character(outlier_distances$annotated_splice_site) == "annotated", ]


	################
	# Annotated
	################
	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	concensus <- "T"
	d_6_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	concensus <- "G"
	d_5_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	concensus <- "A"
	d_4_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	concensus <- "A"
	d_3_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -2
	ss_name <- "D+2"
	concensus <- "T"
	d_2_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -1
	ss_name <- "D+1"
	concensus <- "G"
	d_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	concensus <- "G"
	d_minus_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, TRUE)

	ss_type <- "acceptor"
	position <- -2
	ss_name <- "A-2"
	concensus <- "A"
	a_minus_2_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- -1
	ss_name <- "A-1"
	concensus <- "G"
	a_minus_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	concensus <- "G"
	a_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	################
	# NOVEL
	################
	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	concensus <- "T"
	d_6_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	concensus <- "G"
	d_5_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	concensus <- "A"
	d_4_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	concensus <- "A"
	d_3_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -2
	ss_name <- "D+2"
	concensus <- "T"
	d_2_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -1
	ss_name <- "D+1"
	concensus <- "G"
	d_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	concensus <- "G"
	d_minus_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, TRUE)

	ss_type <- "acceptor"
	position <- -2
	ss_name <- "A-2"
	concensus <- "A"
	a_minus_2_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- -1
	ss_name <- "A-1"
	concensus <- "G"
	a_minus_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	concensus <- "G"
	a_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)




	mutation_bar_plot_legend <- get_broad_mutation_type_bar_plot_legend(inlier_distances_novel, outlier_distances_novel, "acceptor", 0, "A+1", "G", TRUE)

	gg_novel <- plot_grid(d_minus_1_plot, d_1_plot, d_2_plot, d_3_plot, d_4_plot, d_5_plot, d_6_plot, a_minus_2_plot, a_minus_1_plot, a_1_plot, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1, 1, 1, 1, 1))
	gg_ano <- plot_grid(d_minus_1_plot_ano, d_1_plot_ano, d_2_plot_ano, d_3_plot_ano, d_4_plot_ano, d_5_plot_ano, d_6_plot_ano, a_minus_2_plot_ano, a_minus_1_plot_ano, a_1_plot_ano, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1, 1, 1, 1, 1))

	combined_gg <- plot_grid(gg_ano, gg_novel, ncol=1)+
					draw_label("Annotated Splice Site", x=.5,y=.98,size=9) + 
					draw_label("Novel Splice Site", x=.5,y=.49,size=9) 

	return(plot_grid(combined_gg, mutation_bar_plot_legend, nrow=2, rel_heights=c(1,.05)))

}
broad_mutation_type_bar_plot_across_concensus_sites_v2 <- function(inlier_distance_file, outlier_distance_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	concensus <- "T"
	d_6_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	concensus <- "G"
	d_5_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	concensus <- "A"
	d_4_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	concensus <- "A"
	d_3_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)


	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	concensus <- "G"
	d_minus_1_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, TRUE)


	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	concensus <- "G"
	a_1_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	mutation_bar_plot_legend <- get_broad_mutation_type_bar_plot_legend(inlier_distances, outlier_distances, "acceptor", 0, "A+1", "G", TRUE)

	gg <- plot_grid(d_minus_1_plot, d_3_plot, d_4_plot, d_5_plot, d_6_plot, a_1_plot, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1))

	return(plot_grid(gg, mutation_bar_plot_legend, nrow=2, rel_heights=c(1,.1)))

}


broad_mutation_type_bar_plot_across_concensus_sites <- function(inlier_distance_file, outlier_distance_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	concensus <- "T"
	d_6_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	concensus <- "G"
	d_5_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	concensus <- "A"
	d_4_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	concensus <- "A"
	d_3_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -2
	ss_name <- "D+2"
	concensus <- "T"
	d_2_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -1
	ss_name <- "D+1"
	concensus <- "G"
	d_1_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	concensus <- "G"
	d_minus_1_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, TRUE)

	ss_type <- "acceptor"
	position <- -2
	ss_name <- "A-2"
	concensus <- "A"
	a_minus_2_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- -1
	ss_name <- "A-1"
	concensus <- "G"
	a_minus_1_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	concensus <- "G"
	a_1_plot <- broad_mutation_type_bar_plot(inlier_distances, outlier_distances, ss_type, position, ss_name,concensus, FALSE)

	mutation_bar_plot_legend <- get_broad_mutation_type_bar_plot_legend(inlier_distances, outlier_distances, "acceptor", 0, "A+1", "G", TRUE)

	gg <- plot_grid(d_minus_1_plot, d_1_plot, d_2_plot, d_3_plot, d_4_plot, d_5_plot, d_6_plot, a_minus_2_plot, a_minus_1_plot, a_1_plot, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1, 1, 1, 1, 1))

	return(plot_grid(gg, mutation_bar_plot_legend, nrow=2, rel_heights=c(1,.1)))

}


extract_positional_odds_ratio_data_for_all_positions_in_window <- function(outlier_distances, inlier_distances, intron_start_pos, exon_start_pos) {
	dist_to_ss <- c()
	odds_ratios <- c()
	upper_bounds <- c()
	lower_bounds <- c()
	for (distance_iter in -intron_start_pos:(exon_start_pos-1)) {
		rv_outlier <- sum(outlier_distances[,1] == distance_iter)
		rv_inlier <- sum(inlier_distances[,1] == distance_iter)
		no_rv_outlier <- sum(outlier_distances[,1] != distance_iter)
		no_rv_inlier <- sum(inlier_distances[,1] != distance_iter)
		# odds_ratio_old <- (rv_outlier/rv_inlier)/(no_rv_outlier/no_rv_inlier)
		a <- rv_outlier  + 1
		b <- rv_outlier + no_rv_outlier + 2
		c <- rv_inlier +1
		d <- rv_inlier + no_rv_inlier + 2
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

# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites) using odds ratios of real vs background
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
make_positional_odds_ratio_plot_seperated_by_ss_type <- function(inlier_distance_file, outlier_distance_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	outlier_donor_distances <- outlier_distances[as.character(outlier_distances$splice_site_type)=="donor",]
	outlier_acceptor_distances <- outlier_distances[as.character(outlier_distances$splice_site_type)=="acceptor",]
	inlier_donor_distances <- inlier_distances[as.character(inlier_distances$splice_site_type)=="donor",]
	inlier_acceptor_distances <- inlier_distances[as.character(inlier_distances$splice_site_type)=="acceptor",]



	donor_df <- extract_positional_odds_ratio_data_for_all_positions_in_window(outlier_donor_distances, inlier_donor_distances, 7, 3)
	acceptor_df <- extract_positional_odds_ratio_data_for_all_positions_in_window(outlier_acceptor_distances, inlier_acceptor_distances, 3,2)



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

	acceptor_plot <-  ggplot() + geom_point(data=df_acceptor, mapping=aes(x=dist_to_ss, y=odds_ratio), color="#0D324D") +
					geom_errorbar(data=df_acceptor, mapping=aes(x=dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="#0D324D",width=0.0) +
					labs(x = "", y = "") +
					#geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
					#geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() +
					theme(axis.text.x=element_text(angle=45,hjust=1)) +
					scale_x_continuous(breaks=-3:1, labels=c("A-3","A-2","A-1","A+1", "A+2")) +
					theme(axis.title.y=element_blank())


	donor_plot <-  ggplot() + geom_point(data=df_donor, mapping=aes(x=-dist_to_ss, y=odds_ratio), color="#0D324D") +
					geom_errorbar(data=df_donor, mapping=aes(x=-dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="#0D324D",width=0.0) +
					labs(x = "", y = "Relative risk") +
					#geom_vline(xintercept = 2.5, size=.00001,linetype="dashed") +
					#geom_vline(xintercept = .5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() + 
					theme(axis.text.x=element_text(angle=45,hjust=1)) +
					scale_x_continuous(breaks=-2:7, labels=c("D-3", "D-2", "D-1", "D+1", "D+2", "D+3","D+4", "D+5", "D+6", "D+7"))


	error_bar_plot <- plot_grid(donor_plot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),acceptor_plot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol=2, rel_widths=c(1.9,1))

	return(error_bar_plot)


}


# Make density plot of distance between rare variants and splice sites (with seperate plots for 5' and 3' splice sites) using odds ratios of real vs background
# Positive distance corresponds to variant being on exon while negative distance corresponds to variant being in intron
make_positional_odds_ratio_plot_seperated_by_ss_type_no_y_axis <- function(inlier_distance_file, outlier_distance_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	outlier_donor_distances <- outlier_distances[as.character(outlier_distances$splice_site_type)=="donor",]
	outlier_acceptor_distances <- outlier_distances[as.character(outlier_distances$splice_site_type)=="acceptor",]
	inlier_donor_distances <- inlier_distances[as.character(inlier_distances$splice_site_type)=="donor",]
	inlier_acceptor_distances <- inlier_distances[as.character(inlier_distances$splice_site_type)=="acceptor",]



	donor_df <- extract_positional_odds_ratio_data_for_all_positions_in_window(outlier_donor_distances, inlier_donor_distances, 7, 3)
	acceptor_df <- extract_positional_odds_ratio_data_for_all_positions_in_window(outlier_acceptor_distances, inlier_acceptor_distances, 3,2)



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

	acceptor_plot <-  ggplot() + geom_point(data=df_acceptor, mapping=aes(x=dist_to_ss, y=odds_ratio), color="#0D324D") +
					geom_errorbar(data=df_acceptor, mapping=aes(x=dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="#0D324D",width=0.0) +
					labs(x = "", y = "") +
					#geom_vline(xintercept = -.5, size=.00001,linetype="dashed") +
					#geom_vline(xintercept = -2.5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() +
					theme(axis.text.x=element_text(angle=45,hjust=1)) +
					scale_x_continuous(breaks=-3:1, labels=c("A-3","A-2","A-1","A+1", "A+2")) +
					theme(axis.title.y=element_blank())


	donor_plot <-  ggplot() + geom_point(data=df_donor, mapping=aes(x=-dist_to_ss, y=odds_ratio), color="#0D324D") +
					geom_errorbar(data=df_donor, mapping=aes(x=-dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="#0D324D",width=0.0) +
					labs(x = "", y = "") +
					#geom_vline(xintercept = 2.5, size=.00001,linetype="dashed") +
					#geom_vline(xintercept = .5, size=.00001,linetype="dashed") + 
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() + 
					theme(axis.text.x=element_text(angle=45,hjust=1)) +
					scale_x_continuous(breaks=-2:7, labels=c("D-3", "D-2", "D-1", "D+1", "D+2", "D+3","D+4", "D+5", "D+6", "D+7"))


	error_bar_plot <- plot_grid(donor_plot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),acceptor_plot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol=2, rel_widths=c(1.9,1))

	return(error_bar_plot)


}


odds_ratio_concensus_boxplot_across_tissues <- function(aa) {
	df_to_concensus <- aa[as.character(aa$to_or_from_concensus) == "to_concensus",]
	df_from_concensus <- aa[as.character(aa$to_or_from_concensus) == "from_concensus",]
	orat <- c()
	type <- c()
	orat <- c(orat, df_to_concensus$odds_ratio, df_from_concensus$odds_ratio)
	type <- c(type, rep("Consensus Created", length(df_to_concensus$odds_ratio)), rep("Consensus Destroyed", length(df_from_concensus$odds_ratio)))
	df <- data.frame(odds_ratio=log(orat), type=as.factor(type))
	plotter <- ggplot(df, aes(x=type, y=odds_ratio, fill=type)) + geom_boxplot() + 
		labs(x="", y="Junction usage", fill="") + 
		scale_fill_manual(values=c("dodgerblue", "violetred1")) +
		theme(legend.position="none") + 
		geom_hline(yintercept=0) + 
		gtex_v8_figure_theme()
	return(plotter)
}

odds_ratio_ppt_boxplot_across_tissues <- function(aa) {
	df_to_concensus <- aa[as.character(aa$variant_type) == "purine->pyrimidine",]
	df_from_concensus <- aa[as.character(aa$variant_type) == "pyrimidine->purine",]
	orat <- c()
	type <- c()
	orat <- c(orat, df_to_concensus$odds_ratio, df_from_concensus$odds_ratio)
	type <- c(type, rep("Purine to Pyrimidine", length(df_to_concensus$odds_ratio)), rep("Pyrimidine to Purine", length(df_from_concensus$odds_ratio)))
	df <- data.frame(odds_ratio=log(orat), type=as.factor(type))
	plotter <- ggplot(df, aes(x=type, y=odds_ratio, fill=type)) + geom_boxplot() + 
		theme(text = element_text(size=12),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=12), legend.title = element_text(size=11)) +
		labs(x="", y="Junction usage", fill="") + 
		scale_fill_manual(values=c("dodgerblue", "violetred1")) +
		theme(legend.position="none") + 
		geom_hline(yintercept=0) +
		gtex_v8_figure_theme()
	return(plotter)
}

tbt_variant_outlier_enrichment_errorbar_plot <- function(tissue_by_tissue_enrichment_file, color_vector) {
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

	tissue_names <- gsub("_", " ", tissue_names)
	# Add information to data frame
	df <- data.frame(tissue_names=factor(tissue_names), odds_ratios=odds_ratios, lower_bounds=lower_bounds, upper_bounds=upper_bounds)
	error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=tissue_names,ymin=lower_bounds, ymax=upper_bounds), width=0, colour=color_vector) +
					geom_point(data=df, mapping=aes(x=tissue_names, y=odds_ratios), colour=color_vector) +
					labs(x = "Tissue", y = "Relative risk") +
					geom_hline(yintercept=1) + 
					gtex_v8_figure_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))

	return(error_bar_plot)
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

cross_tissue_variant_outlier_enrichment_errorbar_plot <- function(variant_enrichment_dir, stem) {
	distances <- c("2", "4", "6", "8", "10", "100")

	# Initialize vectors
	pvalues <- c()
	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()
	distance_vec <- c()

	# Loop through distances (one enrichment file for each distance)
	for (distance_iter in 1:length(distances)) {
		distance <- distances[distance_iter]
		file_name <- paste0(variant_enrichment_dir, stem, "distance_", distance, "_version_all.txt")
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
	dodge <- position_dodge(width=0.5)
	error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=distance,ymin=lower_bounds, ymax=upper_bounds, colour=pvalues), position=dodge,width=0.0) +
					geom_point(data=df, mapping=aes(x=distance, y=odds_ratios, colour=pvalues), position=dodge) +
					labs(x = "Base pair window around splice site", y = "Relative risk", colour="p-value") +
					geom_hline(yintercept=1) + 
					gtex_v8_figure_theme() + theme(legend.position = c(0.8, 0.8))

	return(error_bar_plot)
}

cross_tissue_variant_outlier_enrichment_mutually_exclusive_errorbar_plot <- function(variant_enrichment_dir, stem) {
	distances <- c("2", "4", "6", "8", "10", "100")
	me_distances <- c("[0,2]", "[3,4]", "[5,6]", "[7,8]", "[9,10]", "[11,100]")

	# Initialize vectors
	pvalues <- c()
	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()
	distance_vec <- c()

	# Loop through distances (one enrichment file for each distance)
	for (distance_iter in 1:length(distances)) {
		distance <- distances[distance_iter]
		me_distance <- me_distances[distance_iter]
		file_name <- paste0(variant_enrichment_dir, stem, "distance_", distance, "_version_all.txt")
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
			distance_vec <- c(distance_vec, me_distance)
		}
	}

	# Make data frame
	df <- data.frame(pvalues=factor(pvalues,levels=(as.character(pvalue_string))), odds_ratios=odds_ratios, lower_bounds=lower_bounds, upper_bounds=upper_bounds,distance=factor(distance_vec,levels=me_distances))
	dodge <- position_dodge(width=0.5)
	error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=distance,ymin=lower_bounds, ymax=upper_bounds, colour=pvalues), position=dodge,width=0.0) +
					geom_point(data=df, mapping=aes(x=distance, y=odds_ratios, colour=pvalues), position=dodge) +
					labs(x = "Base pair window around splice site", y = "Relative risk", colour="p-value") +
					geom_hline(yintercept=1) + 
					gtex_v8_figure_theme() + theme(legend.position = c(0.8, 0.8)) 

	return(error_bar_plot)
}



purine_pyrimidine_variant_changes_in_specific_region <- function(outlier_distance_file, inlier_distance_file, ss_type, ppt_start, ppt_end) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)


	inlier_distances_region <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance > ppt_start & inlier_distances$distance <= ppt_end & as.character(inlier_distances$annotated_splice_site)=="annotated",]
	outlier_distances_region <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance > ppt_start & outlier_distances$distance <= ppt_end & as.character(outlier_distances$annotated_splice_site)=="annotated",]

	num_outliers <- dim(outlier_distances_region)[1]
	num_inliers <- dim(inlier_distances_region)[1]

	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()
	variant_class <- c()

	# Purine -> pyrimidine
	pu_py_outlier <- sum((outlier_distances_region$major_allele == "A" | outlier_distances_region$major_allele == "G") & (outlier_distances_region$variant_allele == "C" | outlier_distances_region$variant_allele == "T"))
	pu_py_inlier <- sum((inlier_distances_region$major_allele == "A" | inlier_distances_region$major_allele == "G") & (inlier_distances_region$variant_allele == "C" | inlier_distances_region$variant_allele == "T"))
	orat <- (pu_py_outlier/num_outliers)/(pu_py_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/pu_py_outlier) - (1.0/num_outliers) + (1.0/pu_py_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Purine to Pyrimidine")

	# Purine -> purine
	pu_pu_outlier <- sum((outlier_distances_region$major_allele == "A" | outlier_distances_region$major_allele == "G") & (outlier_distances_region$variant_allele == "A" | outlier_distances_region$variant_allele == "G"))
	pu_pu_inlier <- sum((inlier_distances_region$major_allele == "A" | inlier_distances_region$major_allele == "G") & (inlier_distances_region$variant_allele == "A" | inlier_distances_region$variant_allele == "G"))
	orat <- (pu_pu_outlier/num_outliers)/(pu_pu_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/pu_pu_outlier) - (1.0/num_outliers) + (1.0/pu_pu_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Purine to Purine")

	# pyrimidine -> purine
	py_pu_outlier <- sum((outlier_distances_region$major_allele == "C" | outlier_distances_region$major_allele == "T") & (outlier_distances_region$variant_allele == "A" | outlier_distances_region$variant_allele == "G"))
	py_pu_inlier <- sum((inlier_distances_region$major_allele == "C" | inlier_distances_region$major_allele == "T") & (inlier_distances_region$variant_allele == "A" | inlier_distances_region$variant_allele == "G"))
	orat <- (py_pu_outlier/num_outliers)/(py_pu_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/py_pu_outlier) - (1.0/num_outliers) + (1.0/py_pu_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Pyrimidine to Purine")


	# pyrimidine -> pyrimidine
	py_py_outlier <- sum((outlier_distances_region$major_allele == "C" | outlier_distances_region$major_allele == "T") & (outlier_distances_region$variant_allele == "C" | outlier_distances_region$variant_allele == "T")) 
	py_py_inlier <- sum((inlier_distances_region$major_allele == "C" | inlier_distances_region$major_allele == "T") & (inlier_distances_region$variant_allele == "C" | inlier_distances_region$variant_allele == "T")) 
	orat <- (py_py_outlier/num_outliers)/(py_py_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/py_py_outlier) - (1.0/num_outliers) + (1.0/py_py_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Pyrimidine to Pyrimidine")

	df <- data.frame(odds_ratio=odds_ratios, dist=1:length(variant_class), lower_bound=lower_bounds, upper_bound=upper_bounds, variant_class=factor(variant_class))



	error_bar_plot_annotated <-  ggplot() + geom_errorbar(data=df, mapping=aes(x=dist,ymin=lower_bound, ymax=upper_bound),color="darkorchid", width=0) +
					geom_point(data=df, mapping=aes(x=dist, y=odds_ratio), color="darkorchid") +
					labs(x = "", y = "Relative risk", title="Annotated splice site") +
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
					scale_x_continuous(breaks=1:length(df$dist), labels=df$variant_class)





	inlier_distances_region <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance > ppt_start & inlier_distances$distance <= ppt_end & as.character(inlier_distances$annotated_splice_site)=="novel",]
	outlier_distances_region <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance > ppt_start & outlier_distances$distance <= ppt_end & as.character(outlier_distances$annotated_splice_site)=="novel",]

	num_outliers <- dim(outlier_distances_region)[1]
	num_inliers <- dim(inlier_distances_region)[1]

	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()
	variant_class <- c()

	# Purine -> pyrimidine
	pu_py_outlier <- sum((outlier_distances_region$major_allele == "A" | outlier_distances_region$major_allele == "G") & (outlier_distances_region$variant_allele == "C" | outlier_distances_region$variant_allele == "T"))
	pu_py_inlier <- sum((inlier_distances_region$major_allele == "A" | inlier_distances_region$major_allele == "G") & (inlier_distances_region$variant_allele == "C" | inlier_distances_region$variant_allele == "T"))
	orat <- (pu_py_outlier/num_outliers)/(pu_py_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/pu_py_outlier) - (1.0/num_outliers) + (1.0/pu_py_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Purine to Pyrimidine")

	# Purine -> purine
	pu_pu_outlier <- sum((outlier_distances_region$major_allele == "A" | outlier_distances_region$major_allele == "G") & (outlier_distances_region$variant_allele == "A" | outlier_distances_region$variant_allele == "G"))
	pu_pu_inlier <- sum((inlier_distances_region$major_allele == "A" | inlier_distances_region$major_allele == "G") & (inlier_distances_region$variant_allele == "A" | inlier_distances_region$variant_allele == "G"))
	orat <- (pu_pu_outlier/num_outliers)/(pu_pu_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/pu_pu_outlier) - (1.0/num_outliers) + (1.0/pu_pu_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Purine to Purine")

	# pyrimidine -> purine
	py_pu_outlier <- sum((outlier_distances_region$major_allele == "C" | outlier_distances_region$major_allele == "T") & (outlier_distances_region$variant_allele == "A" | outlier_distances_region$variant_allele == "G"))
	py_pu_inlier <- sum((inlier_distances_region$major_allele == "C" | inlier_distances_region$major_allele == "T") & (inlier_distances_region$variant_allele == "A" | inlier_distances_region$variant_allele == "G"))
	orat <- (py_pu_outlier/num_outliers)/(py_pu_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/py_pu_outlier) - (1.0/num_outliers) + (1.0/py_pu_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Pyrimidine to Purine")


	# pyrimidine -> pyrimidine
	py_py_outlier <- sum((outlier_distances_region$major_allele == "C" | outlier_distances_region$major_allele == "T") & (outlier_distances_region$variant_allele == "C" | outlier_distances_region$variant_allele == "T")) 
	py_py_inlier <- sum((inlier_distances_region$major_allele == "C" | inlier_distances_region$major_allele == "T") & (inlier_distances_region$variant_allele == "C" | inlier_distances_region$variant_allele == "T")) 
	orat <- (py_py_outlier/num_outliers)/(py_py_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/py_py_outlier) - (1.0/num_outliers) + (1.0/py_py_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Pyrimidine to Pyrimidine")

	df <- data.frame(odds_ratio=odds_ratios, dist=1:length(variant_class), lower_bound=lower_bounds, upper_bound=upper_bounds, variant_class=factor(variant_class))



	error_bar_plot_novel <-  ggplot() + geom_errorbar(data=df, mapping=aes(x=dist,ymin=lower_bound, ymax=upper_bound),color="darkorchid", width=0.0) +
					geom_point(data=df, mapping=aes(x=dist, y=odds_ratio), color="darkorchid") +
					labs(x = "", y = "Relative risk", title="Novel splice site") +
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
					scale_x_continuous(breaks=1:length(df$dist), labels=df$variant_class)

	error_bar_plot <- plot_grid(error_bar_plot_annotated, error_bar_plot_novel, ncol=2, labels=c("A","B"))

	return(error_bar_plot)
}

ppt_enrichment_acceptor <- function(outlier_distance_file, inlier_distance_file, ss_type) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	position_labels <- c()
	odds_ratios <- c()
	upper_bounds <- c()
	lower_bounds <- c()

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -35, -5, ss_type)	
	print(orat_data)
	if (FALSE) {
	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -10, -5, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[A-5, A-9]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -15, -10, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[A-10, A-14]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -20, -15, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[A-15, A-19]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -25, -20, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[A-20, A-24]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -30, -25, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[A-25, A-29]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -35, -30, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[A-30, A-34]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -40, -35, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[A-35, A-39]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -45, -40, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[A-40, A-44]")


	df <- data.frame(odds_ratio=odds_ratios,upper_bound=upper_bounds, lower_bound=lower_bounds, dist_to_ss=1:length(position_labels), position=factor(position_labels, levels=position_labels))

	error_bar_plot <-  ggplot() + geom_errorbar(data=df, mapping=aes(x=dist_to_ss,ymin=lower_bounds, ymax=upper_bounds),color="darkorchid") +
					geom_point(data=df, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
					labs(x = "Distance from splice site (BP)", y = "Relative risk", title="Acceptor Splice Site") +
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() + theme(axis.text.x=element_text(angle=45,hjust=1)) +
					scale_x_continuous(breaks=1:length(position_labels), labels=position_labels)

	return(error_bar_plot)
	}
}

ppt_enrichment_donor <- function(outlier_distance_file, inlier_distance_file, ss_type, output_file) {
	inlier_distances <- read.table(inlier_distance_file, header=TRUE)
	outlier_distances <- read.table(outlier_distance_file, header=TRUE)

	position_labels <- c()
	odds_ratios <- c()
	upper_bounds <- c()
	lower_bounds <- c()


	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -10, -7, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[D+7, D+9]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -15, -10, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[D+10, D+14]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -20, -15, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[D+15, D+19]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -25, -20, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[D+20, D+24]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -30, -25, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[D+25, D+29]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -35, -30, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[D+30, D+34]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -40, -35, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[D+35, D+39]")

	orat_data <- extract_ppt_odds_ratio_data(outlier_distances, inlier_distances, -45, -40, ss_type)
	odds_ratios <- c(odds_ratios, orat_data$odds_ratio)
	upper_bounds <- c(upper_bounds, orat_data$upper_bound)
	lower_bounds <- c(lower_bounds, orat_data$lower_bound)
	position_labels <- c(position_labels, "[D+40, D+44]")


	df <- data.frame(odds_ratio=odds_ratios,upper_bound=upper_bounds, lower_bound=lower_bounds, dist_to_ss=1:length(position_labels), position=factor(position_labels, levels=position_labels))

	error_bar_plot <-  ggplot() + geom_errorbar(data=df, mapping=aes(x=dist_to_ss,ymin=lower_bounds, ymax=upper_bounds),color="darkorchid") +
					geom_point(data=df, mapping=aes(x=dist_to_ss, y=odds_ratio), color="darkorchid") +
					labs(x = "Distance from splice site (BP)", y = "Relative risk", title="Donor Splice Site") +
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() + theme(axis.text.x=element_text(angle=45,hjust=1)) +
					scale_x_continuous(breaks=1:length(position_labels), labels=position_labels)

	return(error_bar_plot)
}


extract_ppt_odds_ratio_data <- function(outlier_distances, inlier_distances, ppt_start, ppt_end, ss_type) {
	aa <- sum(as.character(outlier_distances$splice_site_type)==ss_type & outlier_distances$distance <= ppt_end & outlier_distances$distance > ppt_start) 

	bb <- sum(as.character(outlier_distances$splice_site_type)==ss_type) 

	cc <- sum(as.character(inlier_distances$splice_site_type)==ss_type & inlier_distances$distance <= ppt_end & inlier_distances$distance > ppt_start) 

	dd <- sum(as.character(inlier_distances$splice_site_type)==ss_type) 


	orat <- (aa/bb)/(cc/dd)
	log_bounds <- 1.96*sqrt((1.0/aa) - (1.0/bb) + (1.0/cc) - (1.0/dd))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)

	df <- data.frame(odds_ratio=orat, upper_bound=upper_bound, lower_bound=lower_bound)
	return(df)

}

make_panel_2a <- function(figure_2_ab_data) {
	load(figure_2_ab_data)
	dcols = c('#BFCDE0', '#7F5A83', '#0D324D')
	#dist.data$Type = sapply(dist.data$Type, function(x) ifelse(x == 'SNPs', 'SNVs', as.character(x)))
	dist.data$Type = factor(dist.data$Type, levels=c('SNVs', 'indels', 'SVs'))
	dist.data$Method = sapply(dist.data$Method, function(x) ifelse(x == 'Splicing', 'sOutlier', as.character(x)))
	dist.data$Method = sapply(dist.data$Method, function(x) ifelse(x == 'ASE', 'aseOutlier', as.character(x)))
	dist.data$Method = sapply(dist.data$Method, function(x) ifelse(x == 'MEDZ', 'eOutlier', as.character(x)))
	dist.data$Method = factor(dist.data$Method, levels=c('eOutlier', 'aseOutlier', 'sOutlier'))

	dist.data.filtered1 = dist.data[dist.data$Riskratio > 1 & dist.data$Type=='SNVs',]


	fig_2A1 = ggplot(dist.data.filtered1, aes(x=Window, y=Riskratio,Group=Method)) +
		geom_point(size=1.7,aes(color=Method), position=position_dodge(width=0.5)) +
		geom_errorbar(aes(ymin = Lower, ymax = Upper),size=.2, width=0,position=position_dodge(width=0.5)) +
		theme_bw() + ylab('Relative risk') + xlab('') +
		geom_line(aes(group=Method),size=0.2) +
		geom_hline(yintercept=1,color='grey',size=.2) +
		scale_color_manual(values=dcols) + guides(color=F) +
		theme(axis.text.x=element_text(angle=45,hjust=1),strip.text.x=element_text(size=8)) +
		gtex_v8_figure_theme() +
		theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
			panel.border = element_blank()) +
		facet_wrap(Type~.,scales='free', ncol=3) + theme(strip.background = element_blank())

	dist.data.filtered2 = dist.data[dist.data$Riskratio > 1 & dist.data$Type!='SNVs',]
	fig_2A2 = ggplot(dist.data.filtered2, aes(x=Window, y=Riskratio,Group=Method)) +
		geom_point(size=1.7,aes(color=Method), position=position_dodge(width=0.5)) +
		geom_errorbar(aes(ymin = Lower, ymax = Upper),size=.2, width=0,position=position_dodge(width=0.5)) +
		theme_bw() + ylab('') + geom_line(aes(group=Method),size=0.2) +
		#xlab('Distance upstream from gene (kb)') + ylab('') +
		xlab('') + ylab('') +
		scale_y_log10() + geom_hline(yintercept=1,size=.2,color='grey') +
		scale_color_manual(values=dcols) +
		theme(axis.text.x=element_text(angle=45,hjust=1),
        	axis.title.x=element_text(hjust=0),
        	strip.text.x=element_text(size=8),
        	legend.position=c(0.83,0.89),
        	panel.border = element_blank(),
        	legend.key.size=unit(0.05,'in')) +
        gtex_v8_figure_theme() +
        facet_wrap(Type~.,scales='free', ncol=3) + 
        theme(plot.margin = unit(c(0, .5, 0, 0), "cm")) +
        theme(legend.title=element_blank(), strip.background = element_blank())

    fig_2A <- plot_grid(fig_2A1, fig_2A2, ncol=2, rel_widths=c(.96,2)) + draw_text('Distance upstream from gene (kb)', x=.53, y=.04,size=8)

   return(fig_2A)
}

make_panel_2b <- function(figure_2_ab_data) {
  load(figure_2_ab_data)

plot.tss.data$promoter_motif = sapply(plot.tss.data$promoter_motif, function(x) ifelse(x == 'other_motif', 'Other motif', as.character(x)))
### panel B ###
vcols = c(brewer.pal(11,'Spectral')[1:7],'grey',brewer.pal(11,'Spectral')[8:11])
names(vcols) = unique(plot.tss.data$promoter_motif)
plot.tss.data$ExpBin = factor(plot.tss.data$ExpBin, levels=c('Under', 'Control', 'Over'))
plot.tss.data$promoter_motif = factor(plot.tss.data$promoter_motif, 
                                      levels=rev(c('Other motif', 'Cmyc', 'CTCF', 'E2F4', 'E2F6', 'ELF1', 'Gabp', 'Nrf1', 'Nrsf', 'PU1', 'SP1', 'Srf', 'Tr4', 'USF1', 'Yy1')))
fig_2B = ggplot(plot.tss.data, aes(x=ExpBin,y=NumBin)) + 
  geom_bar(aes(fill=promoter_motif),stat='identity',color='black') + theme_bw() +
  scale_fill_manual(values=vcols) + xlab('eOutlier bin') + 
  ylab('Proportion with RV in promoter') +
  theme(panel.border = element_blank(),
        legend.key.height = unit(0.05, "in")) +
  gtex_v8_figure_theme() + theme(legend.title=element_blank(),
                                 legend.position=c(0.5,0.6))
  return(fig_2B)
}

make_panel_2b_presentation_version <- function(figure_2_ab_data) {
  load(figure_2_ab_data)


### panel B ###
vcols = c(brewer.pal(11,'Spectral')[1:7],'grey',brewer.pal(11,'Spectral')[8:11])
names(vcols) = unique(plot.tss.data$promoter_motif)
plot.tss.data$ExpBin = factor(plot.tss.data$ExpBin, levels=c('Under', 'Control', 'Over'))
plot.tss.data$promoter_motif = factor(plot.tss.data$promoter_motif, 
                                      levels=rev(c('other_motif', 'Cmyc', 'CTCF', 'E2F4', 'E2F6', 'ELF1', 'Gabp', 'Nrf1', 'Nrsf', 'PU1', 'SP1', 'Srf', 'Tr4', 'USF1', 'Yy1')))
fig_2B = ggplot(plot.tss.data, aes(x=ExpBin,y=NumBin)) + 
  geom_bar(aes(fill=promoter_motif),stat='identity',color='black') + theme_bw() +
  scale_fill_manual(values=vcols) + xlab('Outlier bin') + 
  ylab('Proportion with RV in promoter') +
  theme(panel.border = element_blank(),
        legend.key.height = unit(0.05, "in")) +
  gtex_v8_figure_theme() + theme(legend.title=element_blank(),
                                 legend.position="right")
  return(fig_2B)
}


#######################
# Command Line args
#######################
variant_position_enrichment_dir <- args[1]  # Input dir with positional enrichments
jxn_usage_nearby_altered_ss_enrichment_dir <- args[2]  # Second input dir with read count enrichments
variant_enrichment_dir <- args[3]  # Third input dir was enrichments of rare variants nearby splicing outliers
tissue_names_file <- args[4]  # Filename containing names of gtex tissues
tissue_colors_file <- args[5]  # Filename corresponding to colors of gtex tissues
visualize_variant_position_enrichment_dir <- args[6]  # Output dir
splice_site_cartoon <- args[7]
figure_2_ab_data <- args[8]  # File from Nicole with data for making figures 2a and 2b


options(bitmapType = 'cairo', device = 'pdf')

print(splice_site_cartoon)

########################
# Extract tissue color information
########################
# Extract vector tissue names
tissue_names <- as.character(unlist(read.table(tissue_names_file,header=FALSE), use.names=FALSE))
tissue_colors = read.table(tissue_colors_file, header = T, stringsAsFactors = F, sep = "\t")
# Get vector of hex colors in correct order
color_vector <- get_color_vector(tissue_colors, tissue_names)
# Add to coherent data frame
tissues_df <- data.frame(tissue_id=factor(tissue_names, levels=tissue_names), order=1:length(tissue_names))

#######################
# Make panel 2a (code from Nicole)
#########################
panel_2a <- make_panel_2a(figure_2_ab_data)
output_file <- paste0(visualize_variant_position_enrichment_dir, "panel_2a.pdf")
ggsave(panel_2a, file=output_file, width=7.2, height=3,units="in")


#######################
# Make panel 2b (code from Nicole)
#########################
panel_2b <- make_panel_2b(figure_2_ab_data)
output_file <- paste0(visualize_variant_position_enrichment_dir, "panel_2b.pdf")
ggsave(panel_2b, file=output_file, width=7.2, height=4.5,units="in")

#######################
# Make errorbar plot showing enrichment of rare variants nearby splice sites within splicing outliers
#########################
pvalue_threshold = '.00001'
distance = '6'
tissue_by_tissue_enrichment_file <- paste0(variant_enrichment_dir, "tbt_variant_outlier_enrichment_pvalue_", pvalue_threshold, "_distance_", distance, "_version_all.txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "tbt_variant_outlier_enrichment_pvalue_", pvalue_threshold, "_distance_", distance, "_errorbar.pdf")
tbt_variant_outlier_errorbar_plot <- tbt_variant_outlier_enrichment_errorbar_plot(tissue_by_tissue_enrichment_file, color_vector)
ggsave(tbt_variant_outlier_errorbar_plot, file=output_file, width=7.2, height=5, units="in")


#######################
# Make errorbar plot showing enrichment of rare variants nearby splice sites within heuristic approach outliers
#########################
distance = '8'
tissue_by_tissue_enrichment_file <- paste0(variant_enrichment_dir, "tbt_variant_heuristic_outlier_enrichment_distance_", distance, "_version_all.txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "tbt_variant_heuristic_outlier_enrichment_distance_", distance, "_errorbar.pdf")
tbt_variant_heuristic_outlier_errorbar_plot <- tbt_variant_outlier_enrichment_errorbar_plot(tissue_by_tissue_enrichment_file, color_vector)
ggsave(tbt_variant_heuristic_outlier_errorbar_plot, file=output_file, width=7.2, height=5, units="in")


#######################
# Make errorbar plot showing enrichment of rare variants nearby splice sites within multi-tissue-outliers
#########################
stem <- "cross_tissue_variant_outlier_enrichment_"
output_file <- paste0(visualize_variant_position_enrichment_dir, stem, "errorbar.pdf")
cross_tissue_variant_outlier_errorbar_plot <- cross_tissue_variant_outlier_enrichment_errorbar_plot(variant_enrichment_dir, stem)
ggsave(cross_tissue_variant_outlier_errorbar_plot, file=output_file, width=7.2, height=5, units="in")


#######################
# Make errorbar plot showing enrichment of rare variants nearby splice sites within multi-tissue-outliers
#########################
stem <- "cross_tissue_variant_outlier_enrichment_mutually_exclusive_"
output_file <- paste0(visualize_variant_position_enrichment_dir, stem, "errorbar.pdf")
cross_tissue_variant_outlier_mutually_exclusive_errorbar_plot <- cross_tissue_variant_outlier_enrichment_mutually_exclusive_errorbar_plot(variant_enrichment_dir, stem)
ggsave(cross_tissue_variant_outlier_mutually_exclusive_errorbar_plot, file=output_file, width=7.2, height=5, units="in")



#######################
# Make odds ratio enrichment plots showing pyrimidine/purine nature of ppt variants
#########################
distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
ss_type <- "acceptor"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "ppt_enrichment_", ss_type, "_", distance, "_", pvalue_threshold, "_purine_pyrimidine_enrichment.pdf")
ppt_start <- -35
ppt_end <- -5
purine_pyrimidine_variant_change_enrichment_errorbar <- purine_pyrimidine_variant_changes_in_specific_region(outlier_distance_file, inlier_distance_file, ss_type, ppt_start, ppt_end)
ggsave(purine_pyrimidine_variant_change_enrichment_errorbar, file=output_file, width=7.2, height=3.2,units="in")


#######################
#Make plot showing read count enrichments for rare variants in concensus sites
#########################
# Load in data
input_file <- paste0(jxn_usage_nearby_altered_ss_enrichment_dir, "tissue_by_tissue_outliers_with_rv_in_concensus_sites_outlier_individuals_1e-05_inlier_individuals_1e-05_with_read_counts.txt")
enrichment_data <- read.table(input_file, header=TRUE)
# Visualize enrichment distribution (data aggregrated) tissues
output_file <- paste0(visualize_variant_position_enrichment_dir, "tbt_concensus_jxn_usage_enrichment.pdf")
tbt_concensus_read_count_boxplot <- odds_ratio_concensus_boxplot_across_tissues(enrichment_data)
ggsave(tbt_concensus_read_count_boxplot, file=output_file, width=7.2, height=3,units="in")


#######################
#Make plot showing read count enrichments for rare variants in concensus sites
#########################
# Load in data
input_file <- paste0(jxn_usage_nearby_altered_ss_enrichment_dir, "cross_tissue_outliers_with_rv_in_concensus_sites_outlier_individuals_1e-05_inlier_individuals_1e-05_with_median_read_counts.txt")
enrichment_data <- read.table(input_file, header=TRUE)
# Visualize enrichment distribution (data aggregrated) tissues
output_file <- paste0(visualize_variant_position_enrichment_dir, "cross_tissue_concensus_median_jxn_usage_enrichment.pdf")
xt_median_concensus_read_count_boxplot <- odds_ratio_concensus_boxplot_across_tissues(enrichment_data)
ggsave(xt_median_concensus_read_count_boxplot, file=output_file, width=7.2, height=3,units="in")





#######################
#Make PWM plot showing mutation types across the annotated splice concensus sites
#########################
distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
input_file <- paste0(jxn_usage_nearby_altered_ss_enrichment_dir, "cross_tissue_outliers_with_rv_in_concensus_sites_outlier_individuals_1e-05_inlier_individuals_1e-05_with_median_read_counts.txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_type_across_annotated_splice_site_positions_pwm.pdf" )
mutation_type_annotated_pwm_plot_two_row <- mutation_type_pwm_plot_across_annotated_concensus_sites_two_row(input_file)
#ggsave(mutation_type_annotated_pwm_plot, file=output_file, width=7.2, height=4.0,units="in")

#######################
#Make PWM plot showing mutation types across the annotated splice concensus sites
#########################
distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
input_file <- paste0(jxn_usage_nearby_altered_ss_enrichment_dir, "cross_tissue_outliers_with_rv_in_concensus_sites_outlier_individuals_1e-05_inlier_individuals_1e-05_with_median_read_counts.txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_type_across_annotated_splice_site_positions_pwm.pdf" )
mutation_type_annotated_pwm_plot <- mutation_type_pwm_plot_across_annotated_concensus_sites(input_file)
#ggsave(mutation_type_annotated_pwm_plot, file=output_file, width=7.2, height=4.0,units="in")

#######################
#Make PWM plot showing mutation types across the novel splice concensus sites
#########################
distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
input_file <- paste0(jxn_usage_nearby_altered_ss_enrichment_dir, "cross_tissue_outliers_with_rv_in_concensus_sites_outlier_individuals_1e-05_inlier_individuals_1e-05_with_median_read_counts.txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_type_across_novel_splice_site_positions_pwm.pdf" )
mutation_type_novel_pwm_plot <- mutation_type_pwm_plot_across_novel_concensus_sites(input_file)
#ggsave(mutation_type_novel_pwm_plot, file=output_file, width=7.2, height=4.0,units="in")



#######################
#Make plot showing read count enrichments for rare variants in PPTs
#########################
# Load in data
input_file <- paste0(jxn_usage_nearby_altered_ss_enrichment_dir, "tissue_by_tissue_outliers_with_rv_in_ppt_sites_outlier_individuals_1e-05_inlier_individuals_1e-05_with_read_counts.txt")
enrichment_data <- read.table(input_file, header=TRUE)
# Visualize enrichment distribution (data aggregrated) tissues
output_file <- paste0(visualize_variant_position_enrichment_dir, "tbt_ppt_jxn_usage_enrichment.pdf")
tbt_ppt_read_count_enrichment_boxplot <- odds_ratio_ppt_boxplot_across_tissues(enrichment_data)
ggsave(tbt_ppt_read_count_enrichment_boxplot, file=output_file, width=7.2, height=3,units="in")


#######################
#Make plot showing read count enrichments for rare variants in PPTs
#########################
# Load in data
input_file <- paste0(jxn_usage_nearby_altered_ss_enrichment_dir, "cross_tissue_outliers_with_rv_in_ppt_sites_outlier_individuals_1e-05_inlier_individuals_1e-05_with_median_read_counts.txt")
enrichment_data <- read.table(input_file, header=TRUE)
# Visualize enrichment distribution (data aggregrated) tissues
output_file <- paste0(visualize_variant_position_enrichment_dir, "cross_tissue_ppt_median_jxn_usage_enrichment.pdf")
xt_median_ppt_read_count_enrichment_boxplot <- odds_ratio_ppt_boxplot_across_tissues(enrichment_data)
ggsave(xt_median_ppt_read_count_enrichment_boxplot, file=output_file, width=7.2, height=3,units="in")



#######################
# Positional odds ratio plots for RV seperated by donor vs acceptor splice sites
#########################
distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_ss_type_seperated_positional_odds_ratio.pdf")
panel_2d <- make_positional_odds_ratio_plot_seperated_by_ss_type(inlier_distance_file, outlier_distance_file)
ggsave(panel_2d, file=paste0(visualize_variant_position_enrichment_dir, "panel_2d.pdf"), width=7.2, height=3, units="in")

#######################
# Positional odds ratio plots for RV seperated by donor vs acceptor splice sites
#########################
distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
panel_2d_no_y_axis <- make_positional_odds_ratio_plot_seperated_by_ss_type_no_y_axis(inlier_distance_file, outlier_distance_file)
#ggsave(panel_2d, file=paste0(visualize_variant_position_enrichment_dir, "panel_2d.pdf"), width=7.2, height=3, units="in")


#######################
#Make plot showing read count enrichments for rare variants in PPTs
#########################
# Load in data
input_file <- paste0(jxn_usage_nearby_altered_ss_enrichment_dir, "cross_tissue_outliers_with_rv_in_ppt_sites_outlier_individuals_1e-05_inlier_individuals_1e-05_with_median_read_counts.txt")
enrichment_data <- read.table(input_file, header=TRUE)
# Visualize enrichment distribution (data aggregrated) tissues
output_file <- paste0(visualize_variant_position_enrichment_dir, "panel_2e.pdf")
panel_2e <- odds_ratio_ppt_boxplot_across_tissues(enrichment_data)
ggsave(panel_2e, file=output_file, width=7.2, height=3,units="in")


############
# Make figure 2 
#########################
first_row <- plot_grid(panel_2a, panel_2b, labels = c('A','B'), ncol=2)
second_row <- plot_grid(NULL, labels = c('C'), ncol=1)
third_row <- plot_grid(panel_2d, panel_2e, labels = c('D','E'), ncol=2, rel_widths=c(1.6,1))
fourth_row <- plot_grid(mutation_type_annotated_pwm_plot_two_row, labels = c('F'), ncol=1)
fig_2 = plot_grid(first_row,second_row,third_row, fourth_row, nrow = 4, align='v', rel_heights=c(1.3,.55,.9,1.4))
ggsave(fig_2, file=paste0(visualize_variant_position_enrichment_dir,'fig2.pdf'), width=7.2, height=8,units="in")


############
# Make figure 2 
#########################
first_row <- plot_grid(panel_2a, panel_2b, labels = c('A','B'), ncol=2)
second_row <- plot_grid(NULL, labels = c('C'), ncol=1)
third_row <- plot_grid(panel_2d, panel_2e, labels = c('D','F'), ncol=2, rel_widths=c(1.6,1))
fourth_row <- plot_grid(NULL, mutation_type_annotated_pwm_plot, labels = c('E', ''), ncol=2, rel_widths=c(.01,1))
fig_2 = plot_grid(first_row,second_row,third_row, fourth_row, nrow = 4, align='v', rel_heights=c(1.1,.55,.85,.5))
ggsave(fig_2, file=paste0(visualize_variant_position_enrichment_dir,'fig2_v2.pdf'), width=7.2, height=6.7,units="in")


############
# Make figure 2 
#########################

panel_2e_edited <- panel_2e + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.border = element_blank())
first_row <- plot_grid(panel_2a, panel_2b, labels = c('A','B'), ncol=2)
second_row <- plot_grid(NULL, labels = c('C'), ncol=1)
third_row <- plot_grid(NULL,NULL, labels = c('D', ''),rel_heights=c(.04,1), ncol=1)
fourth_row <- plot_grid(NULL, NULL, mutation_type_annotated_pwm_plot_two_row, NULL, labels = c("E","F", "", ""), ncol=2, rel_heights=c(.03,1), rel_widths=c(1,.7))
fig_2 = plot_grid(first_row,second_row,third_row, fourth_row, nrow = 4, align='v', rel_heights=c(.98,.45,.63,.78))
combined_fig2 <- ggdraw() + 
				draw_plot(fig_2,0,0,1,1) +
				draw_plot(panel_2e_edited,.59,-.019,.4,.255) +
				draw_plot(panel_2d_no_y_axis,-.0145,.242,1.00,.22) +
				draw_label("Relative risk", 0.007, .39, size=8,angle=90) 

ggsave(combined_fig2, file=paste0(visualize_variant_position_enrichment_dir,'fig2_v3.pdf'), width=7.2, height=7.0,units="in")


############
# Make figure 2 
#########################

panel_2e_edited <- panel_2e + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.border = element_blank())
first_row <- plot_grid(panel_2a, panel_2b, labels = c('A','B'), ncol=2)
second_row <- plot_grid(NULL, labels = c('C'), ncol=1)
third_row <- plot_grid(NULL,NULL, labels = c('D', ''),rel_heights=c(.04,1), ncol=1)
fourth_row <- plot_grid(NULL, NULL, mutation_type_annotated_pwm_plot_two_row, NULL, labels = c("E","F", "", ""), ncol=2, rel_heights=c(.03,1), rel_widths=c(1,.56))
fig_2 = plot_grid(first_row,second_row,third_row, fourth_row, nrow = 4, align='v', rel_heights=c(.98,.45,.63,.78))
combined_fig2 <- ggdraw() + 
				draw_plot(fig_2,0,0,1,1) +
				draw_plot(panel_2e_edited + labs(y=""),.632,-.019,.357,.255) +
				draw_plot(panel_2d_no_y_axis,-.0145,.242,1.00,.22) +
				draw_label("Relative risk", 0.007, .39, size=8,angle=90) +
				draw_label("Junction usage", 0.655, .145, size=8,angle=90) 

ggsave(combined_fig2, file=paste0(visualize_variant_position_enrichment_dir,'fig2_v4.pdf'), width=7.2, height=7.0,units="in")



print("DONE")



#######################
# Make TBT enrichment supplementary figure
#########################
#splice_concensus_site_enrichment_plot <- plot_grid(tbt_variant_outlier_errorbar_plot, plot_grid(cross_tissue_variant_outlier_mutually_exclusive_errorbar_plot, xt_median_concensus_read_count_boxplot, labels=c("B","C"), ncol=2, rel_widths=c(.6,.4)), ncol=1, rel_heights=c(.6,.44), labels=c("A",""))
#output_file <- paste0(visualize_variant_position_enrichment_dir, "splice_consensus_site_enrichment_supplementary_figure.pdf")
#ggsave(splice_concensus_site_enrichment_plot, file=output_file, width=7.2, height=6.0, units="in")
tbt_enrichment_plot <- plot_grid(tbt_variant_outlier_errorbar_plot, plot_grid(tbt_concensus_read_count_boxplot, tbt_ppt_read_count_enrichment_boxplot, labels=c("B","C"), ncol=2, rel_widths=c(.5,.5)), ncol=1, rel_heights=c(.6,.44), labels=c("A",""))
output_file <- paste0(visualize_variant_position_enrichment_dir, "tissue_by_tissue_enrichment_supplementary_figure.pdf")
ggsave(tbt_enrichment_plot, file=output_file, width=7.2, height=6.0, units="in")


#######################
# Make Positional enrichment supplementary figure
#########################
xt_enrichment_plot <- plot_grid(cross_tissue_variant_outlier_mutually_exclusive_errorbar_plot, xt_median_concensus_read_count_boxplot, labels=c("A","B"), ncol=2, rel_widths=c(.6,.4))
output_file <- paste0(visualize_variant_position_enrichment_dir, "cross_tissue_enrichment_supplementary_figure.pdf")
ggsave(xt_enrichment_plot, file=output_file, width=7.2, height=3.0, units="in")




########################################
# Supplementary figure for pwm for annotated and novel
########################################
#combined_pwm_plot <- plot_grid(mutation_type_annotated_pwm_plot, mutation_type_novel_pwm_plot,ncol=1, labels=c("A","B"))
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_variant_type_across_splice_site_positions_pwm_supplement.pdf" )
ggsave(mutation_type_novel_pwm_plot, file=output_file, width=7.2, height=3.0,units="in")

#######################
#Make plot showing broad mutation types across the splice concensus sites sepeated by novel and rare
#########################
distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
output_file <- paste0(visualize_variant_position_enrichment_dir, "distance_to_",version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, "_broad_variant_type_across_splice_site_positions_seperated_by_novel_and_rare.pdf" )
broad_mutation_type_bar_plot_seperated_by_novel_and_rare_v1 <- broad_mutation_type_bar_plot_across_concensus_sites_seperated_by_novel_and_rare_v1(inlier_distance_file, outlier_distance_file)
ggsave(broad_mutation_type_bar_plot_seperated_by_novel_and_rare_v1, file=output_file, width=7.2, height=3,units="in")


distance <- "1000"
version <- "observed_splice_site"
pvalue_threshold <- "1e-05"
ss_type <- "acceptor"
outlier_distance_file <- paste0(variant_position_enrichment_dir, "outlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
inlier_distance_file <- paste0(variant_position_enrichment_dir, "inlier_distance_to_", version, "_distance_", distance, "_pvalue_thresh_", pvalue_threshold, ".txt")
ppt_enrichment_acceptor(outlier_distance_file, inlier_distance_file, ss_type)


