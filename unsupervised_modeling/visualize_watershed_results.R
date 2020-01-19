args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(PRROC)
library(cowplot)
library(RColorBrewer)



get_tissue_names <- function(roc_object, number_of_dimensions) {
	tissue_names <- c()
	for (tissue_num in 1:number_of_dimensions) {
		tissue_name <- strsplit(roc_object$roc[[tissue_num]]$name, "_total_expression")[[1]][1]
		tissue_names <- c(tissue_names, tissue_name)
	}
	return(tissue_names)
}



make_theta_pair_heatmap <- function(theta_pair, number_of_dimensions,tissue_names, outlier_type) {
	# Convert theta_pair vector into matrix of number_of_tissuesXnumber_of_tissues
	theta_pair_mat = matrix(0, number_of_dimensions, number_of_dimensions)
	dimension_counter = 1
	for (dimension1 in 1:number_of_dimensions) {
		for (dimension2 in dimension1:number_of_dimensions) {
			if (dimension1 != dimension2) {
				theta_pair_mat[dimension1, dimension2] = theta_pair[1, dimension_counter]
				theta_pair_mat[dimension2, dimension1] = theta_pair[1, dimension_counter]
				dimension_counter = dimension_counter + 1
			}
		}
		theta_pair_mat[dimension1, dimension1] = NA
	}

	# Cluster tissues based on similarity of theta_pairs
	#order <- hclust( as.dist(1- abs(corr_mat)), method = "ward.D" )$order
	order <- hclust( dist(theta_pair_mat, method = "euclidean"), method = "ward.D" )$order

	rownames(theta_pair_mat) <- gsub("_", " ", tissue_names, fixed=TRUE)
	colnames(theta_pair_mat) <- gsub("_", " ", tissue_names, fixed=TRUE)
	#gsub("_"," ",ordered_tissue_names,fixed=TRUE)

	melted_mat <- melt(theta_pair_mat)
    # Axis labels are factors
    melted_mat$X1 <- factor(melted_mat$X1)
    melted_mat$X2 <- factor(melted_mat$X2)


    #  Use factors to represent covariate and pc name
    melted_mat$X1 <- factor(melted_mat$X1, levels = rownames(theta_pair_mat)[order])
    melted_mat$X2 <- factor(melted_mat$X2, levels = rownames(theta_pair_mat)[order])


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    heatmap <- heatmap + scale_fill_gradient2() #+ scale_fill_distiller(palette="RdBu")
    heatmap <- heatmap + theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5,size=6), axis.text.y = element_text(size=6))
    heatmap <- heatmap + gtex_v8_figure_theme() + theme(axis.text.x=element_blank())
    heatmap <- heatmap + labs(x = "Tissue", y = "Tissue",fill="Edge Weight", title=outlier_type)
    return(heatmap)

}

make.point.plot = function(tissuesdf, colors, vertical = TRUE){
    if (vertical) {
        p = ggplot(tissuesdf, aes(x = 1, y = order, label = tissue_id))
    } else {
        p = ggplot(tissuesdf, aes(x = order, y = 1))
    }
    p = p + geom_point(aes(colour = tissue_id), size = .4) +
        scale_colour_manual(values = colors) + guides(colour = FALSE) + 
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank())
    return(p)
}

make_fig_5_theta_pair_heatmap <- function(theta_pair, number_of_dimensions,tissue_names, tissue_colors) {
	# Convert theta_pair vector into matrix of number_of_tissuesXnumber_of_tissues
	theta_pair_mat = matrix(0, number_of_dimensions, number_of_dimensions)
	dimension_counter = 1
	for (dimension1 in 1:number_of_dimensions) {
		for (dimension2 in dimension1:number_of_dimensions) {
			if (dimension1 != dimension2) {
				theta_pair_mat[dimension1, dimension2] = theta_pair[1, dimension_counter]
				theta_pair_mat[dimension2, dimension1] = theta_pair[1, dimension_counter]
				dimension_counter = dimension_counter + 1
			}
		}
		theta_pair_mat[dimension1, dimension1] = NA
	}

	# Cluster tissues based on similarity of theta_pairs
	#order <- hclust( as.dist(1- abs(corr_mat)), method = "ward.D" )$order
	order <- hclust( dist(theta_pair_mat, method = "euclidean"), method = "ward.D" )$order

	#######
	tissue_colors$tissue_id = factor(tissue_colors$tissue_id, levels = as.character(tissue_names[order]))
	#print(tissue_colors$tissue_id)
	tissue_colors$order = as.numeric(tissue_colors$tissue_id)
	tissue_colors = tissue_colors[!is.na(tissue_colors$order),]
	colors = paste0("#",tissue_colors$tissue_color_hex)
	names(colors) = tissue_colors$tissue_id

	colors.vertical = make.point.plot(tissue_colors, colors)
	colors.horizontal = make.point.plot(tissue_colors, colors, vertical = FALSE)



	#######

	rownames(theta_pair_mat) <- tissue_names
	colnames(theta_pair_mat) <- tissue_names

	melted_mat <- melt(theta_pair_mat)
    # Axis labels are factors
    melted_mat$X1 <- factor(melted_mat$X1)
    melted_mat$X2 <- factor(melted_mat$X2)


    #  Use factors to represent covariate and pc name
    melted_mat$X1 <- factor(melted_mat$X1, levels = rownames(theta_pair_mat)[order])
    melted_mat$X2 <- factor(melted_mat$X2, levels = rownames(theta_pair_mat)[order])


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    heatmap <- heatmap + scale_fill_gradient2() #+ scale_fill_distiller(palette="RdBu")
    heatmap <- heatmap + theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5,size=7), axis.text.y = element_text(size=7))
    heatmap <- heatmap + gtex_v8_figure_theme() + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
    heatmap <- heatmap + labs(x = "Tissue", y = "Tissue",fill="Edge Weight")
    
   	combined = ggdraw() +
        draw_plot(heatmap, .05,.05,.95,0.95) +
        draw_plot(colors.vertical, -.321, .123, .903, .903) + 
        draw_plot(colors.horizontal,0.086,-.157, .69,.69)

    return(combined)

}


make_delta_auc_plot <- function(ase_watershed, ase_river, ase_tissue_names, splicing_watershed, splicing_river, splicing_tissue_names, te_watershed, te_river, te_tissue_names, tissue_names_subset) {
	delta_auc <- c()
	outlier_type <- c()
	# ASE
	for (tissue_num in 1:length(ase_tissue_names)) {
		tissue_name <- ase_tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			watershed_auc <- ase_watershed[[tissue_num]]$evaROC$watershed_pr_auc
			river_auc <- ase_river[[tissue_num]]$evaROC$watershed_pr_auc
			delta_auc <- c(delta_auc, watershed_auc-river_auc)
			outlier_type <- c(outlier_type, "ASE")
		}
	}
	# Splicing
	for (tissue_num in 1:length(splicing_tissue_names)) {
		tissue_name <- splicing_tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			watershed_auc <- splicing_watershed[[tissue_num]]$evaROC$watershed_pr_auc
			river_auc <- splicing_river[[tissue_num]]$evaROC$watershed_pr_auc
			delta_auc <- c(delta_auc, watershed_auc-river_auc)
			outlier_type <- c(outlier_type, "Splicing")
		}
	}
	# TE
	for (tissue_num in 1:length(te_tissue_names)) {
		tissue_name <- te_tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			watershed_auc <- te_watershed[[tissue_num]]$evaROC$watershed_pr_auc
			river_auc <- te_river[[tissue_num]]$evaROC$watershed_pr_auc
			delta_auc <- c(delta_auc, watershed_auc-river_auc)
			outlier_type <- c(outlier_type, "TE")
		}
	}
	# Put all into compact data frame
	df <- data.frame(outlier_type=factor(outlier_type), delta_auc=delta_auc)

	# Make boxplot
	boxplot <- ggplot(df, aes(x=outlier_type, y=delta_auc, fill=outlier_type)) + geom_boxplot() +
				gtex_v8_figure_theme() + 
				xlab("") +
               	ylab("Watershed AUC - River AUC") + 
               	theme(legend.position="none") +
               	scale_fill_manual(values=c("#7F5A83", "#0D324D", "#BFCDE0")) + 
               	geom_hline(yintercept = 0)

	return(boxplot)


}


make_tbt_auc_distribution_plot <- function(ase_watershed, ase_river, ase_tissue_names, splicing_watershed, splicing_river, splicing_tissue_names, te_watershed, te_river, te_tissue_names, tissue_names_subset) {
	auc <- c()
	outlier_type <- c()
	method <- c()
	# ASE
	for (tissue_num in 1:length(ase_tissue_names)) {
		tissue_name <- ase_tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			watershed_auc <- ase_watershed[[tissue_num]]$evaROC$watershed_pr_auc
			river_auc <- ase_river[[tissue_num]]$evaROC$GAM_pr_auc
			auc <- c(auc, watershed_auc, river_auc)
			outlier_type <- c(outlier_type, "ASE", "ASE")
			method <- c(method, "Watershed", "GAM")
			# delta_auc <- c(delta_auc, watershed_auc-river_auc)
			# outlier_type <- c(outlier_type, "ASE")
		}
	}
	# Splicing
	for (tissue_num in 1:length(splicing_tissue_names)) {
		tissue_name <- splicing_tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			watershed_auc <- splicing_watershed[[tissue_num]]$evaROC$watershed_pr_auc
			river_auc <- splicing_river[[tissue_num]]$evaROC$GAM_pr_auc
			auc <- c(auc, watershed_auc, river_auc)
			outlier_type <- c(outlier_type, "Splicing", "Splicing")
			method <- c(method, "watershed", "GAM")
		}
	}
	# TE
	for (tissue_num in 1:length(te_tissue_names)) {
		tissue_name <- te_tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			watershed_auc <- te_watershed[[tissue_num]]$evaROC$watershed_pr_auc
			river_auc <- te_river[[tissue_num]]$evaROC$GAM_pr_auc
			auc <- c(auc, watershed_auc, river_auc)
			outlier_type <- c(outlier_type, "TE", "TE")
			method <- c(method, "Watershed", "GAM")
		}
	}
	# Put all into compact data frame
	df <- data.frame(outlier_type=factor(outlier_type), auc=auc, method=factor(method,levels=c("GAM", "Watershed")))

	# Make boxplot
	boxplot <- ggplot(df, aes(x=method, y=auc, fill=outlier_type)) + geom_boxplot() +
				gtex_v8_figure_theme() + 
				xlab("") +
               	ylab("AUC(PR)") + 
               	theme(legend.position="none") +
               	scale_fill_manual(values=c("#7F5A83", "#0D324D", "#BFCDE0"))  
               	#geom_hline(yintercept = 0)

	return(boxplot)


}



make_tbt_auc_distribution_plot2 <- function(ase_watershed, ase_river, ase_tissue_names, splicing_watershed, splicing_river, splicing_tissue_names, te_watershed, te_river, te_tissue_names, tissue_names_subset) {
	auc <- c()
	outlier_type <- c()
	method <- c()
	# ASE
	for (tissue_num in 1:length(ase_tissue_names)) {
		tissue_name <- ase_tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			watershed_auc <- ase_watershed[[tissue_num]]$evaROC$watershed_pr_auc
			gam_auc <- ase_river[[tissue_num]]$evaROC$GAM_pr_auc
			river_auc <- ase_river[[tissue_num]]$evaROC$watershed_pr_auc
			auc <- c(auc, watershed_auc, river_auc, gam_auc)
			outlier_type <- c(outlier_type, "ASE","ASE", "ASE")
			method <- c(method, "Watershed", "RIVER", "GAM")
			# delta_auc <- c(delta_auc, watershed_auc-river_auc)
			# outlier_type <- c(outlier_type, "ASE")
		}
	}
	# Splicing
	for (tissue_num in 1:length(splicing_tissue_names)) {
		tissue_name <- splicing_tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			watershed_auc <- splicing_watershed[[tissue_num]]$evaROC$watershed_pr_auc
			gam_auc <- splicing_river[[tissue_num]]$evaROC$GAM_pr_auc
			river_auc <- splicing_river[[tissue_num]]$evaROC$watershed_pr_auc
			auc <- c(auc, watershed_auc, river_auc, gam_auc)
			outlier_type <- c(outlier_type, "Splicing", "Splicing", "Splicing")
			method <- c(method, "Watershed","RIVER", "GAM")
		}
	}
	# TE
	for (tissue_num in 1:length(te_tissue_names)) {
		tissue_name <- te_tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			watershed_auc <- te_watershed[[tissue_num]]$evaROC$watershed_pr_auc
			gam_auc <- te_river[[tissue_num]]$evaROC$GAM_pr_auc
			river_auc <- te_river[[tissue_num]]$evaROC$watershed_pr_auc
			auc <- c(auc, watershed_auc, river_auc, gam_auc)
			outlier_type <- c(outlier_type, "Expression", "Expression", "Expression")
			method <- c(method, "Watershed", "RIVER","GAM")
		}
	}
	# Put all into compact data frame
	df <- data.frame(outlier_type=factor(outlier_type, levels=c("ASE","Splicing", "Expression")), auc=auc, method=factor(method,levels=c("GAM", "RIVER", "Watershed")))

	# Make boxplot
	boxplot <- ggplot(df, aes(x=method, y=auc, fill=outlier_type)) + geom_boxplot() +
				gtex_v8_figure_theme() + 
				xlab("") +
               	ylab("AUC(PR)") + 
               	theme(legend.position="none") +
               	scale_fill_manual(values=c("#7F5A83", "#0D324D", "#BFCDE0"))  
               	#geom_hline(yintercept = 0)

	return(boxplot)


}

make_tbt_auc_lolipop_plot <- function(roc_object_vi, roc_object_exact, tissue_names, outlier_type, tissue_names_subset) {
	auc_vi <- c()
	auc_exact <- c()
	ordered_tissue_names <- c()
	for (tissue_num in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			tissue_auc_vi <- roc_object_vi[[tissue_num]]$evaROC$watershed_pr_auc
			tissue_auc_independent <- roc_object_exact[[tissue_num]]$evaROC$watershed_pr_auc
			auc_vi <- c(auc_vi, tissue_auc_vi)
			auc_exact <- c(auc_exact, tissue_auc_independent)
			ordered_tissue_names <- c(ordered_tissue_names, tissue_name)
		}
	}

	#"firebrick4", "dodgerblue3"
	#cols <- c( "c1" = rgb(0.2,0.7,0.1,0.5), "c2" = rgb(0.7,0.2,0.1,0.5) )
	cols <- c( "c1" = "steelblue3", "c2" = "firebrick4" )
	df <- data.frame(auc_watershed = auc_vi, auc_river=auc_exact, tissue=ordered_tissue_names, tissues_position=1:length(ordered_tissue_names))
	plotter <- ggplot(df) +
  			   geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_vi, yend=auc_exact), color="grey") +
  			   geom_point( aes(x=tissues_position, y=auc_vi, color="c1"), size=1.5) +
               geom_point( aes(x=tissues_position, y=auc_exact, color="c2"), size=1.5) +
               scale_color_manual(name="", breaks=c("c1","c2"), values=cols, labels=c("Watershed", "River")) +
               #ggtitle(outlier_type) +
                xlab("") +
               ylab("AUC (PR)") + 
               theme(legend.position="bottom") +
               theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=0.5)) +
               gtex_v8_figure_theme() +
               #scale_y_continuous(expand = c(0, 0), limits = c(0, 1.06)) + 
               scale_x_continuous(breaks=1:length(ordered_tissue_names),labels=gsub("_"," ",ordered_tissue_names,fixed=TRUE)) 
    return(plotter +draw_text(outlier_type, x=11, y=.68,size=8))
}
#gsub(".", " ", data, fixed=TRUE)

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}



make_tbt_auc_lolipop_plot_across_outliers_on_same_axis <- function(roc_watershed_splicing, roc_river_splicing, roc_watershed_te, roc_river_te, roc_watershed_ase, roc_river_ase, splicing_tissue_names, te_tissue_names, ase_tissue_names, output_file) {
	options(bitmapType = 'cairo', device = 'pdf')
	auc_watershed <- c()
	auc_river <- c()
	tissue_names <- c()
	outlier_type <- c()
	positions <- c()
	ase_num <- 1
	position_counter <- 1
	for (tissue_num in 1:length(splicing_tissue_names)) {
		tissue_name <- splicing_tissue_names[tissue_num]
		tissue_names <- c(tissue_names, tissue_name, tissue_name, tissue_name)
		outlier_type <- c(outlier_type, "ASE", "Splice", "TE")
		positions <- c(positions, position_counter, position_counter + 1, position_counter + 2)
		position_counter <- position_counter + 9

		if (tissue_name == "Kidney_Cortex") {
			auc_watershed <- c(auc_watershed, NaN, roc_watershed_splicing[[tissue_num]]$evaROC$watershed_pr_auc, roc_watershed_te[[tissue_num]]$evaROC$watershed_pr_auc)
			auc_river <- c(auc_river, NaN, roc_river_splicing[[tissue_num]]$evaROC$watershed_pr_auc, roc_river_te[[tissue_num]]$evaROC$watershed_pr_auc)
		} else {
			auc_watershed <- c(auc_watershed, roc_watershed_ase[[ase_num]]$evaROC$watershed_pr_auc, roc_watershed_splicing[[tissue_num]]$evaROC$watershed_pr_auc, roc_watershed_te[[tissue_num]]$evaROC$watershed_pr_auc)
			auc_river <- c(auc_river, roc_river_ase[[ase_num]]$evaROC$watershed_pr_auc, roc_river_splicing[[tissue_num]]$evaROC$watershed_pr_auc, roc_river_te[[tissue_num]]$evaROC$watershed_pr_auc)
			ase_num <- ase_num + 1
		}

	}
	df <- data.frame(auc_watershed = auc_watershed, auc_river=auc_river, tissue=tissue_names, tissues_position=positions, outlier_type=factor(outlier_type))
	cols <- c( "c1" = rgb(0.2,0.7,0.1,0.5), "c2" = rgb(0.1,0.2,0.6,0.4) )

	plotter <- ggplot(df) +
  			   geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_watershed, yend=auc_river), color="grey") +
  			   geom_point( aes(x=tissues_position, y=auc_watershed, color="c1", shape=outlier_type), size=1) +
               geom_point( aes(x=tissues_position, y=auc_river, color="c2", shape=outlier_type), size=1) +
               scale_color_manual(name="", breaks=c("c1","c2"), values=cols, labels=c("watershed-tbt", "river-median")) +
               #coord_flip()+
                xlab("Tissue") +
               ylab("AUC (PR)") + 
               theme(legend.position="bottom") +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
               gtex_v8_figure_theme()
               #scale_x_continuous(breaks=seq(1,length(tissue_names),3),labels=rep("", length(splicing_tissue_names)))
    ggsave(plotter, file=output_file, width=7.2, height=5.5, units="in")

}

make_tbt_auc_lolipop_plot_v2 <- function(roc_object_vi, roc_object_exact, tissue_names, outlier_type, tissue_names_subset) {
	auc_vi <- c()
	auc_exact <- c()
	ordered_tissue_names <- c()

	for (tissue_num in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			tissue_auc_vi <- roc_object_vi[[tissue_num]]$evaROC$watershed_pr_auc
			tissue_auc_independent <- roc_object_exact[[tissue_num]]$evaROC$watershed_pr_auc
			auc_vi <- c(auc_vi, tissue_auc_vi)
			auc_exact <- c(auc_exact, tissue_auc_independent)
			ordered_tissue_names <- c(ordered_tissue_names, tissue_name)
		}
	}

	cols <- c( "c1" = "steelblue3", "c2" = "firebrick4" )

	df <- data.frame(auc_watershed = auc_vi, auc_river=auc_exact, tissue=ordered_tissue_names, tissues_position=1:length(ordered_tissue_names))
	plotter <- ggplot(df) +
  			   geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_vi, yend=auc_exact), color="grey") +
  			   geom_point( aes(x=tissues_position, y=auc_vi, color="c1"), size=1.5) +
               geom_point( aes(x=tissues_position, y=auc_exact, color="c2"), size=1.5) +
               scale_color_manual(name="", breaks=c("c1","c2"), values=cols, labels=c("Watershed", "River")) +
               ggtitle(outlier_type) +
                xlab("") +
               ylab("AUC (PR)") + 
               theme(legend.position="bottom") +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
               gtex_v8_figure_theme() +
               #scale_y_continuous(expand = c(0, 0), limits = c(0, 1.06)) + 
               scale_x_continuous(breaks=1:length(ordered_tissue_names),labels=gsub("_"," ",ordered_tissue_names,fixed=TRUE)) 
    return(plotter)
}

make_tbt_auc_lolipop_plot_with_median_river <- function(roc_object, tissue_names, outlier_type, tissue_names_subset) {
	auc_vi <- c()
	auc_exact <- c()
	ordered_tissue_names <- c()

	for (tissue_num in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_num]
		if (tissue_name %in% tissue_names_subset) {
			tissue_auc_vi <- roc_object[[tissue_num]]$evaROC$watershed_pr_auc
			tissue_auc_independent <- roc_object[[tissue_num]]$evaROC$median_river_pr_auc
			auc_vi <- c(auc_vi, tissue_auc_vi)
			auc_exact <- c(auc_exact, tissue_auc_independent)
			ordered_tissue_names <- c(ordered_tissue_names, tissue_name)
		}
	}

	cols <- c( "c1" = "steelblue3", "c2" = rgb(0.2,0.7,0.1,0.5) )

	df <- data.frame(auc_watershed = auc_vi, auc_river=auc_exact, tissue=ordered_tissue_names, tissues_position=1:length(ordered_tissue_names))
	plotter <- ggplot(df) +
  			   geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_vi, yend=auc_exact), color="grey") +
  			   geom_point( aes(x=tissues_position, y=auc_vi, color="c1"), size=1.5) +
               geom_point( aes(x=tissues_position, y=auc_exact, color="c2"), size=1.5) +
               scale_color_manual(name="", breaks=c("c1","c2"), values=cols, labels=c("Watershed", "River")) +
               ggtitle(outlier_type) +
                xlab("") +
               ylab("AUC (PR)") + 
               theme(legend.position="bottom") +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
               gtex_v8_figure_theme() +
               #scale_y_continuous(expand = c(0, 0), limits = c(0, 1.06)) + 
               scale_x_continuous(breaks=1:length(ordered_tissue_names),labels=gsub("_"," ",ordered_tissue_names,fixed=TRUE)) 
    return(plotter)
}

make_conditional_phi_bar_plot <- function(phi, conditional_z, tissue_names, outlier_type) {
	ordered_tissue_names <- c()
	probability <- c()
	outlier_score <- c()
	for (tissue_num in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_num]
		probability <- c(probability, phi[tissue_num,1], phi[tissue_num,2], phi[tissue_num,3])
		ordered_tissue_names <- c(ordered_tissue_names, tissue_name, tissue_name, tissue_name)
		outlier_score <- c(outlier_score, 1, 2, 3)
		if (tissue_name == "Heart_Left_Ventricle" & outlier_type == "ase") {
			ordered_tissue_names <- c(ordered_tissue_names, "Kidney_Cortex", "Kidney_Cortex", "Kidney_Cortex")
			probability <- c(probability, NaN, NaN, NaN)
			outlier_score <- c(outlier_score, 1, 2, 3)
		}
	}
	df <- data.frame(probability=probability, tissue=ordered_tissue_names, outlier_score=factor(outlier_score))
	plotter <- ggplot() +
	           geom_bar(aes(y = probability, x = tissue, fill = outlier_score), data = df, stat="identity") +
	           ylab(paste0("P(E | Z=", conditional_z, ")")) +
	           xlab("") +
	           theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
	           gtex_v8_figure_theme() +
	           scale_fill_manual(values=c("steelblue3","mediumorchid3","tomato3"))
	return(plotter)
}

make_phi_bar_plot <- function(phi, tissue_names, outlier_type) {
	conditional_1_plot <- make_conditional_phi_bar_plot(phi$outlier_component, 1, tissue_names, outlier_type)
	conditional_0_plot <- make_conditional_phi_bar_plot(phi$inlier_component, 0, tissue_names, outlier_type)
	legend <- get_legend(conditional_0_plot + theme(legend.position="bottom"))

	combined <- plot_grid(conditional_0_plot + ggtitle(outlier_type) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), conditional_1_plot +theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none"), ncol=1, rel_heights=c(1,2))
	
	combined2 <- ggdraw() + draw_plot(combined,0,0,1,1) + draw_plot(legend,.37,-.465,1,1)

	return(combined2)
}

extract_betas <- function(eids, tissue_index, theta, version) {
	betas <- c()
	for (eid_index in 1:length(eids)) {
		eid <- eids[eid_index]
		anno_namer <- paste0('chrom_hmm_', eid, "_", version)
		anno_index <- which(rownames(theta) == anno_namer)
		if (length(anno_index) != 1) {
			print("ASSUMPTION ERROR")
		}
		betas <- c(betas, theta[anno_index[1], tissue_index])
	}
	return(betas)

}

make_tissue_specific_theta <- function(theta, anno_names, tissue_names, chrom_hmm_to_tissue, version) {
	# Need to also input tissue-names and castel file
	rownames(theta) = anno_names
	unique_tissues = unique(as.character(chrom_hmm_to_tissue$GTEx_tissue))
	unique_EIDs = unique(as.character(chrom_hmm_to_tissue$EID))

	betas <- c()
	type <- c()
	tissue_type <- c()

	for (tissue_num in 1:length(unique_tissues)) {
		tissue_name <- unique_tissues[tissue_num]
		indices <- which(as.character(chrom_hmm_to_tissue$GTEx_tissue) == tissue_name)
		matched_eids <- as.character(chrom_hmm_to_tissue$EID)[indices]
		unmatched_eids <- unique_EIDs[unique_EIDs != matched_eids]
		tissue_index <- which(tissue_names==tissue_name)
		tissue_matched_betas <- extract_betas(matched_eids, tissue_index, theta, version)
		tissue_unmatched_betas <- extract_betas(unmatched_eids, tissue_index, theta, version)

		betas <- c(betas, tissue_matched_betas, tissue_unmatched_betas)
		type <- c(type, rep("matched", length(tissue_matched_betas)), rep("unmatched", length(tissue_unmatched_betas)))
		tissue_type <- c(tissue_type, rep(tissue_name, length(tissue_unmatched_betas) + length(tissue_matched_betas)))
		print(tissue_name)
		print(mean(tissue_matched_betas))
		print(mean(tissue_unmatched_betas))
	}

	df <- data.frame(beta=betas, type=factor(type), tissue=factor(tissue_type))

	matched <- df$beta[as.character(type) == "matched"]
	unmatched <- df$beta[as.character(type) == "unmatched"]

	#print(mean(matched))
	#print(mean(unmatched))
	#print(wilcox.test(matched, unmatched))
	#for (tissue_num in 1:dim(chrom_hmm_to_tissue)[1]) {
	#	tissue_name <- as.character(chrom_hmm_to_tissue$GTEx_tissue[tissue_num])
	#	tissue_index = which(tissue_names==tissue_name)
	#	EID = as.character(chrom_hmm_to_tissue$EID[tissue_num])
	#}
}

plot_three_class_theta_pair_term <- function(theta_pair) {
	mat <- matrix(0, 3, 3)
	#mat[1,2] <- theta_pair[1,1]
	#mat[2,1] <- theta_pair[1,1]
	#mat[1,3] <- theta_pair[1,2]
	#mat[3,1] <- theta_pair[1,2]
	#mat[3,2] <- theta_pair[1,3]
	#mat[2,3] <- theta_pair[1,3]
	mat[1,2] <- theta_pair[1,2]
	mat[2,1] <- theta_pair[1,2]
	mat[1,3] <- theta_pair[1,3]
	mat[3,1] <- theta_pair[1,3]
	mat[3,2] <- theta_pair[1,1]
	mat[2,3] <- theta_pair[1,1]

	melted_corr <- melt(mat)

	
    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1)
    melted_corr$X2 <- factor(melted_corr$X2)

    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) #+ scale_fill_gradient(low="grey",high="plum2")

    heatmap <- heatmap + scale_fill_distiller(palette = "Blues", direction=1, na.value = "white")

    heatmap <- heatmap + theme(axis.text.x = element_text(angle = 0, vjust=.5),legend.position="bottom") 
    heatmap <- heatmap + gtex_v8_figure_theme()
    heatmap <- heatmap + labs(fill="Edge weight",x = "", y="")

    heatmap <- heatmap + scale_x_discrete(breaks=c("1", "2", "3"),labels=c("ASE", "Splicing", "Expression"))
    heatmap <- heatmap + scale_y_discrete(breaks=c("1", "2", "3"),labels=c("ASE", "Splicing", "Expression"))

    legend <- get_legend(heatmap)

    #combined <- plot_grid(heatmap+theme(legend.position="none"), legend, ncol=1,rel_heights=c(1,.15))
    combined <- ggdraw() + draw_plot(heatmap+theme(legend.position="none"),0,.15,1,.85) + draw_plot(legend,.2,-.31,1,1)

	return(combined)
}

#######################################
## Compute absolute risk for high watershed posterior variants
#######################################
absolute_risk_plot <- function(gam_predictions, watershed_predictions, data_input, pvalue) {
	#######################################
	# Load in all data (training and test)
	#######################################
	N2_pairs <- data_input$N2_pairs	

  	real_valued_outliers_test2 <- abs(rbind(data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),]))

  	watershed_threshold <- .9

  	thresh <- c()
  	absolute_risks <- c()
  	outlier_type <- c()
  	model_type <- c()


  		# splice
  		ws_variants <- watershed_predictions[,1] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,1] <= pvalue)/length(real_valued_outliers_test2[ws_variants,1])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "sOutlier")
  		model_type <- c(model_type, "Watershed")

  		# TE
  		ws_variants <- watershed_predictions[,2] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,2] <= pvalue)/length(real_valued_outliers_test2[ws_variants,2])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "eOutlier")
  		model_type <- c(model_type, "Watershed")

  		# ase
  		ws_variants <- watershed_predictions[,3] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,3] <= pvalue)/length(real_valued_outliers_test2[ws_variants,3])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "aseOutlier")
  		model_type <- c(model_type, "Watershed")

  		 # splice
  		nvar <- sum(watershed_predictions[,1] > watershed_threshold)
  		matched_thresh <- rev(sort(gam_predictions[,1]))[nvar+1]
  		ws_variants <- gam_predictions[,1] > matched_thresh
  		risk <- sum(real_valued_outliers_test2[ws_variants,1] <= pvalue)/length(real_valued_outliers_test2[ws_variants,1])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "sOutlier")
  		model_type <- c(model_type, "GAM")

  		# TE
  		nvar <- sum(watershed_predictions[,2] > watershed_threshold)
  		matched_thresh <- rev(sort(gam_predictions[,2]))[nvar+1]
  		ws_variants <- gam_predictions[,2] > matched_thresh
  		risk <- sum(real_valued_outliers_test2[ws_variants,2] <= pvalue)/length(real_valued_outliers_test2[ws_variants,2])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "eOutlier")
  		model_type <- c(model_type, "GAM")

  		# ase
   		nvar <- sum(watershed_predictions[,3] > watershed_threshold)
  		matched_thresh <- rev(sort(gam_predictions[,3]))[nvar+1]
  		ws_variants <- gam_predictions[,3] > matched_thresh
  		risk <- sum(real_valued_outliers_test2[ws_variants,3] <= pvalue)/length(real_valued_outliers_test2[ws_variants,3])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "aseOutlier")
  		model_type <- c(model_type, "GAM")


    df <- data.frame(watershed_threshold=thresh, absolute_risk=absolute_risks, outlier_type=factor(outlier_type, levels=c("aseOutlier","sOutlier","eOutlier")), model_type=factor(model_type, levels=c("GAM", "Watershed")))
	p <- ggplot(data=df, aes(x=model_type, y=absolute_risk, fill=outlier_type)) +
	geom_bar(stat="identity", color="black", position=position_dodge())+
  	gtex_v8_figure_theme() + 
  	labs(x="", y="Proportion of variants\nleading to outlier", fill="") + 
  	scale_fill_manual(values=c("#7F5A83", "#0D324D", "#BFCDE0"))

  	legend <- get_legend(p)

  	combined <- ggdraw() + draw_plot(p + theme(legend.position="none"),0, -.1, 1,1.1) + draw_plot(legend, .26, .33, 1,1)


	return(combined)
}


#######################################
## Compute absolute risk for high watershed posterior variants
#######################################
absolute_risk_plot_with_cadd <- function(gam_predictions, watershed_predictions, cadd_predictions, data_input, pvalue, watershed_threshold) {
	plot_1 <- absolute_risk_plot_with_cadd_helper(gam_predictions, watershed_predictions, cadd_predictions, data_input, pvalue, .5)
	plot_2 <- absolute_risk_plot_with_cadd_helper(gam_predictions, watershed_predictions, cadd_predictions, data_input, pvalue, .7)
	plot_3 <- absolute_risk_plot_with_cadd_helper(gam_predictions, watershed_predictions, cadd_predictions, data_input, pvalue, .9)
	legend <- get_legend(plot_1 + theme(legend.position="bottom"))
	combined <- plot_grid(plot_1 + theme(legend.position="none"), plot_2 + theme(legend.position="none"), plot_3 + theme(legend.position="none"), legend, labels=c("A","B","C"), rel_heights=c(1,1,1,.13), nrow=4)
	return(combined)
}
absolute_risk_plot_with_cadd_helper <- function(gam_predictions, watershed_predictions, cadd_predictions, data_input, pvalue, watershed_threshold) {
	#######################################
	# Load in all data (training and test)
	#######################################
	N2_pairs <- data_input$N2_pairs	

  	real_valued_outliers_test2 <- abs(rbind(data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),]))


  	thresh <- c()
  	absolute_risks <- c()
  	outlier_type <- c()
  	model_type <- c()


  		# splice
  		ws_variants <- watershed_predictions[,1] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,1] <= pvalue)/length(real_valued_outliers_test2[ws_variants,1])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "sOutlier")
  		model_type <- c(model_type, "Watershed")

  		# TE
  		ws_variants <- watershed_predictions[,2] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,2] <= pvalue)/length(real_valued_outliers_test2[ws_variants,2])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "eOutlier")
  		model_type <- c(model_type, "Watershed")

  		# ase
  		ws_variants <- watershed_predictions[,3] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,3] <= pvalue)/length(real_valued_outliers_test2[ws_variants,3])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "aseOutlier")
  		model_type <- c(model_type, "Watershed")

  		 # splice
  		nvar <- sum(watershed_predictions[,1] > watershed_threshold)
  		matched_thresh <- rev(sort(gam_predictions[,1]))[nvar+1]
  		ws_variants <- gam_predictions[,1] > matched_thresh
  		risk <- sum(real_valued_outliers_test2[ws_variants,1] <= pvalue)/length(real_valued_outliers_test2[ws_variants,1])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "sOutlier")
  		model_type <- c(model_type, "GAM")

  		# TE
  		nvar <- sum(watershed_predictions[,2] > watershed_threshold)
  		matched_thresh <- rev(sort(gam_predictions[,2]))[nvar+1]
  		ws_variants <- gam_predictions[,2] > matched_thresh
  		risk <- sum(real_valued_outliers_test2[ws_variants,2] <= pvalue)/length(real_valued_outliers_test2[ws_variants,2])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "eOutlier")
  		model_type <- c(model_type, "GAM")

  		# ase
   		nvar <- sum(watershed_predictions[,3] > watershed_threshold)
  		matched_thresh <- rev(sort(gam_predictions[,3]))[nvar+1]
  		ws_variants <- gam_predictions[,3] > matched_thresh
  		risk <- sum(real_valued_outliers_test2[ws_variants,3] <= pvalue)/length(real_valued_outliers_test2[ws_variants,3])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "aseOutlier")
  		model_type <- c(model_type, "GAM")


  		 # splice
  		nvar <- sum(watershed_predictions[,1] > watershed_threshold)
  		matched_thresh <- rev(sort(cadd_predictions))[nvar+1]
  		ws_variants <- cadd_predictions > matched_thresh
  		risk <- sum(real_valued_outliers_test2[ws_variants,1] <= pvalue)/length(real_valued_outliers_test2[ws_variants,1])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "sOutlier")
  		model_type <- c(model_type, "CADD")

  		# TE
  		nvar <- sum(watershed_predictions[,2] > watershed_threshold)
  		matched_thresh <- rev(sort(cadd_predictions))[nvar+1]
  		ws_variants <- cadd_predictions > matched_thresh
  		risk <- sum(real_valued_outliers_test2[ws_variants,2] <= pvalue)/length(real_valued_outliers_test2[ws_variants,2])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "eOutlier")
  		model_type <- c(model_type, "CADD")

  		# ase
   		nvar <- sum(watershed_predictions[,3] > watershed_threshold)
  		matched_thresh <- rev(sort(cadd_predictions))[nvar+1]
  		ws_variants <- cadd_predictions > matched_thresh
  		risk <- sum(real_valued_outliers_test2[ws_variants,3] <= pvalue)/length(real_valued_outliers_test2[ws_variants,3])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "aseOutlier")
  		model_type <- c(model_type, "CADD")


    df <- data.frame(watershed_threshold=thresh, absolute_risk=absolute_risks, outlier_type=factor(outlier_type, levels=c("aseOutlier","sOutlier","eOutlier")), model_type=factor(model_type, levels=c("CADD", "GAM", "Watershed")))
	print(watershed_threshold)
	print(df)
	p <- ggplot(data=df, aes(x=model_type, y=absolute_risk, fill=outlier_type)) +
	geom_bar(stat="identity", color="black", position=position_dodge())+
	ylim(0,.7) + 
  	gtex_v8_figure_theme() + 
  	labs(x="", y="Proportion of variants\nleading to outlier", fill="", title=paste0("Watershed posterior > ", watershed_threshold)) + 
  	scale_fill_manual(values=c("#7F5A83", "#0D324D", "#BFCDE0")) 
  	return(p)
}


plot_pr_gam_river_watershed_comparison_curve <- function(roc_object_exact, roc_object_ind, number_of_dimensions, output_file) {
	precision <- c()
	recall <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object_exact[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed (exact inference)
		precision <- c(precision, dimension_roc_object$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("Watershed", length(dimension_roc_object$evaROC$watershed_precision)))


		# GAM
		precision <- c(precision, dimension_roc_object$evaROC$GAM_precision)
		recall <- c(recall, dimension_roc_object$evaROC$GAM_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$GAM_precision)))
		prediction_type <- c(prediction_type, rep("GAM", length(dimension_roc_object$evaROC$GAM_precision)))


		# Indepdent
		precision <- c(precision, dimension_roc_object_ind$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object_ind$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("RIVER", length(dimension_roc_object_ind$evaROC$watershed_precision)))

	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("Watershed","RIVER", "GAM")))
  

	outlier_type <- "ase"
	outlier_name <- "ASE"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision", group="", linetype="", colour="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_linetype_manual(values=c("solid", "dotted", "solid")) +
                gtex_v8_figure_theme() + 
                draw_label(outlier_name,x=.5,y=.95,size=8)

	outlier_type <- "splicing"
	outlier_name <- "Splicing"
  	plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision",group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_linetype_manual(values=c("solid", "dotted", "solid")) +
                gtex_v8_figure_theme() +
                 draw_label(outlier_name,x=.5,y=.95,size=8)

	outlier_type <- "total_expression"
	outlier_name <- "Expression"
  	plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision", group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_linetype_manual(values=c("solid", "dotted", "solid")) +
                gtex_v8_figure_theme() + 
                draw_label(outlier_name,x=.5,y=.95,size=8)

    legend <- get_legend(plotter_ase + theme(legend.position="bottom"))
    combined_plots <- plot_grid(plotter_ase + theme(legend.position="none"), plotter_splice+ theme(legend.position="none"), plotter_te+ theme(legend.position="none"), rel_widths=c(1,1,1), nrow=1)

	combined <- ggdraw() + draw_plot(combined_plots,0,.07,1,.9) + draw_plot(legend,.38,-0.38,1,1)

	return(combined)
}

plot_pr_cadd_watershed_comparison_curve <- function(roc_object_exact, roc_object_ind, number_of_dimensions, output_file) {
	precision <- c()
	recall <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object_exact[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed (exact inference)
		precision <- c(precision, dimension_roc_object$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("Watershed", length(dimension_roc_object$evaROC$watershed_precision)))


		# GAM
		precision <- c(precision, dimension_roc_object$evaROC$CADD_precision)
		recall <- c(recall, dimension_roc_object$evaROC$CADD_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$CADD_precision)))
		prediction_type <- c(prediction_type, rep("CADD", length(dimension_roc_object$evaROC$CADD_precision)))


	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("Watershed","CADD")))
  

	outlier_type <- "ase"
	outlier_name <- "ASE"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision", group="", linetype="", colour="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                scale_linetype_manual(values=c("solid", "solid")) +
                gtex_v8_figure_theme() + 
                draw_label(outlier_name,x=.5,y=.95,size=8)

	outlier_type <- "splicing"
	outlier_name <- "Splicing"
  	plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision",group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                scale_linetype_manual(values=c("solid", "solid")) +
                gtex_v8_figure_theme() +
                 draw_label(outlier_name,x=.5,y=.95,size=8)

	outlier_type <- "total_expression"
	outlier_name <- "Expression"
  	plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision", group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                scale_linetype_manual(values=c("solid", "solid")) +
                gtex_v8_figure_theme() + 
                draw_label(outlier_name,x=.5,y=.95,size=8)

    legend <- get_legend(plotter_ase + theme(legend.position="bottom"))
    combined_plots <- plot_grid(plotter_ase + theme(legend.position="none"), plotter_splice+ theme(legend.position="none"), plotter_te+ theme(legend.position="none"), rel_widths=c(1,1,1), nrow=1)

	combined <- ggdraw() + draw_plot(combined_plots,0,.07,1,.9) + draw_plot(legend,.44,-0.42,1,1)

	return(combined)
}




get_tissues_with_many_n2_pairs <- function(splicing_obj, te_obj, ase_obj, splicing_tissue_names, te_tissue_names, ase_tissue_names, min_num_pos) {
	dicti = list()
	for (tissue_num in 1:length(splicing_tissue_names)) {
		tissue_name <- splicing_tissue_names[tissue_num]
		pos <- splicing_obj[[tissue_num]]$evaROC$num_positive_pairs
		neg <- splicing_obj[[tissue_num]]$evaROC$num_negative_pairs
		if (pos >= min_num_pos) {
			if (tissue_name %in% names(dicti)) {
				dicti[[tissue_name]] = as.numeric(dicti[[tissue_name]]) + 1
			} else {
				dicti[[tissue_name]] = 1
			}
		}
	}

	for (tissue_num in 1:length(te_tissue_names)) {
		tissue_name <- te_tissue_names[tissue_num]
		pos <- te_obj[[tissue_num]]$evaROC$num_positive_pairs
		neg <- te_obj[[tissue_num]]$evaROC$num_negative_pairs
		if (pos >= min_num_pos) {
			if (tissue_name %in% names(dicti)) {
				dicti[[tissue_name]] = as.numeric(dicti[[tissue_name]]) + 1
			} else {
				dicti[[tissue_name]] = 1
			}
		}
	}

	for (tissue_num in 1:length(ase_tissue_names)) {
		tissue_name <- ase_tissue_names[tissue_num]
		pos <- ase_obj[[tissue_num]]$evaROC$num_positive_pairs
		neg <- ase_obj[[tissue_num]]$evaROC$num_negative_pairs
		if (pos >= min_num_pos) {
			if (tissue_name %in% names(dicti)) {
				dicti[[tissue_name]] = as.numeric(dicti[[tissue_name]]) + 1
			} else {
				dicti[[tissue_name]] = 1
			}
		}
	}

	counter = 0
	tissue_subset <- c()
	for (tissue_num in 1:length(te_tissue_names)) {
		tissue_name <- te_tissue_names[tissue_num]
		if (tissue_name %in% names(dicti)) {
			count = as.numeric(dicti[[tissue_name]])
		} else {
			count = 0
		}
		if (count == 3) {
			tissue_subset <- c(tissue_subset, tissue_name)
		}
	}
	return(tissue_subset)

}



plot_pr_river_watershed_comparison_curve <- function(roc_object_exact, roc_object_vi, roc_object_ind, number_of_dimensions) {
	precision <- c()
	recall <- c()
	outlier_type <- c()
	prediction_type <- c()
	inference_method <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object_exact[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		dimension_roc_object_vi <- roc_object_vi[[dimension]]
		# Tied watershed (exact inference)
		precision <- c(precision, dimension_roc_object$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("Watershed", length(dimension_roc_object$evaROC$watershed_precision)))
		inference_method <- c(inference_method, rep("Exact", length(dimension_roc_object$evaROC$watershed_precision)))


		# Tied watershed (variational inference)
		precision <- c(precision, dimension_roc_object_vi$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object_vi$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_vi$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("Watershed", length(dimension_roc_object_vi$evaROC$watershed_precision)))
		inference_method <- c(inference_method, rep("Approximate", length(dimension_roc_object_vi$evaROC$watershed_precision)))


		# Indepdent
		precision <- c(precision, dimension_roc_object_ind$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object_ind$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("RIVER", length(dimension_roc_object_ind$evaROC$watershed_precision)))
		inference_method <- c(inference_method, rep("Exact", length(dimension_roc_object_ind$evaROC$watershed_precision)))

	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), inference=factor(inference_method, levels=c("Exact","Approximate")), prediction_type=factor(prediction_type, levels=c("Watershed","RIVER")))
  

	outlier_type <- "ase"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type, linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title="ASE") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                gtex_v8_figure_theme()
                #theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	outlier_type <- "splicing"
  	plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type, linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title="Splicing") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                gtex_v8_figure_theme()
                #theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	outlier_type <- "total_expression"
  	plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type,linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title="Expression") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                gtex_v8_figure_theme()
                #theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

    legend <- get_legend(plotter_ase)
    combined_plots <- plot_grid(plotter_ase + theme(legend.position="none"),plotter_splice+ theme(legend.position="none"), plotter_te + theme(legend.position="none"), nrow=1)


	return(plot_grid(combined_plots,legend,ncol=1,rel_heights=c(1,.06)))
}

plot_beta_difference_scatter_between_exact_and_vi <- function(model_params_exact, model_params_approximate) {
	exact_betas <- c()
	approximate_betas <- c()
	outlier_class <- c()

	exact_betas <- c(model_params_exact$theta[,1], model_params_exact$theta[,2], model_params_exact$theta[,3])
	approximate_betas <- c(model_params_approximate$theta[,1], model_params_approximate$theta[,2], model_params_approximate$theta[,3])
	outlier_class <- c(rep("Splicing", length(model_params_approximate$theta[,1])), rep("Expression", length(model_params_approximate$theta[,2])), rep("ASE", length(model_params_approximate$theta[,3])))

	df <- data.frame(exact_betas=exact_betas, approximate_betas=approximate_betas, outlier_class=factor(outlier_class, levels=c("ASE", "Splicing", "Expression")))

	plotter <- ggplot(df, aes(x=exact_betas, y=approximate_betas, colour=outlier_class)) + geom_point() +
			geom_abline() + 
			theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)) +
			labs(x="Beta (exact)", y="Beta (approximate)",colour="") + 
			scale_color_manual(values=c("#7F5A83", "#0D324D", "#BFCDE0")) +
			gtex_v8_figure_theme()
	return(plotter)
}

visualize_confusion_matrix <- function(confusion_matrix, titler) {
	confusion_matrix = confusion_matrix/rowSums(confusion_matrix)
    melted_corr <- melt(confusion_matrix)

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1, levels=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"))
    melted_corr$X2 <- factor(melted_corr$X2, levels=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"))
    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) #+ scale_fill_gradient(low="grey",high="plum2")

    heatmap <- heatmap + scale_fill_distiller(palette = "Blues", direction=1, limits=c(0,1))

    #heatmap <- heatmap + theme(text = element_text(size=12),axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11), axis.text.x = element_text(angle = 0, vjust=.5)) 
    heatmap <- heatmap + labs(x = "Observed Class", y = "Predicted Class",fill="",title=titler)
    heatmap <- heatmap + gtex_v8_figure_theme() + theme(axis.text.x = element_text(angle = 0, vjust=.5))
    heatmap <- heatmap + scale_x_discrete(breaks=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"),labels=c("0 0 0", "1 0 0", "0 1 0", "0 0 1", "1 0 1","0 1 1", "1 1 0", "1 1 1"))
    heatmap <- heatmap + scale_y_discrete(breaks=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"),labels=c("0 0 0", "1 0 0", "0 1 0", "0 0 1", "1 0 1","0 1 1", "1 1 0", "1 1 1"))

    return(heatmap)
}

visualize_river_and_watershed_confusion_matrices <- function(watershed_confusion_matrix, watershed_vi_confusion_matrix, river_confusion_matrix) {
	watershed_confusion_plot <- visualize_confusion_matrix(watershed_confusion_matrix, "Watershed (Exact)")
	watershed_vi_confusion_plot <- visualize_confusion_matrix(watershed_vi_confusion_matrix, "Watershed (Approximate)")

	river_confusion_plot <- visualize_confusion_matrix(river_confusion_matrix, "RIVER")


	combined_confusion_plots <- plot_grid(river_confusion_plot, watershed_confusion_plot, watershed_vi_confusion_plot, ncol=1)
	return(combined_confusion_plots)


}

extract_watershed_sample_names <- function(data_input) {
	feat_all <- data_input$feat
	N2_pairs <- data_input$N2_pairs
	feat_train <- feat_all[is.na(N2_pairs),]
	return(rownames(feat_train))
}

extract_sample_outlier_values <- function(data_input) {
	N2_pairs <- data_input$N2_pairs
	outliers_discrete_train <- data_input$outliers_discrete[is.na(N2_pairs),]
	return(outliers_discrete_train)
}

compare_watershed_posteriors_with_different_training_inputs <- function(roc_3_class_data_input, alt_roc_3_class_data_input, roc_object_exact, alt_roc_object_exact) {
	roc_3_class_sample_names <- extract_watershed_sample_names(roc_3_class_data_input)
	alt_roc_3_class_sample_names <- extract_watershed_sample_names(alt_roc_3_class_data_input)
	indices <- alt_roc_3_class_sample_names %in% roc_3_class_sample_names
	if (sum(roc_3_class_sample_names != alt_roc_3_class_sample_names[indices]) != 0) {
		print("FUNDAMENTAL ASSUMPTION ERROR in compare_watershed_posteriors_with_different_training_inputs")
	}

	roc_posteriors <- roc_object_exact$model_params$posterior
	alt_roc_posteriors <- alt_roc_object_exact$model_params$posterior[indices,]

	standard_posterior <- c()
	alt_posterior <- c()
	class <- c()
	num_samples <- length(roc_posteriors[,1])

	standard_posterior <- c(standard_posterior, roc_posteriors[,1])
	alt_posterior <- c(alt_posterior, alt_roc_posteriors[,1])
	class <- c(class, rep("Splicing", num_samples))

	standard_posterior <- c(standard_posterior, roc_posteriors[,2])
	alt_posterior <- c(alt_posterior, alt_roc_posteriors[,2])
	class <- c(class, rep("Expression", num_samples))

	standard_posterior <- c(standard_posterior, roc_posteriors[,3])
	alt_posterior <- c(alt_posterior, alt_roc_posteriors[,3])
	class <- c(class, rep("ASE", num_samples))

	df <- data.frame(standard_posterior=standard_posterior, alt_posterior=alt_posterior, outlier_class=factor(class,levels=c("ASE", "Splicing", "Expression")))

	plotter <- ggplot(df, aes(x=standard_posterior, y=alt_posterior, colour=outlier_class)) + geom_point(size=.1) +
			geom_abline() + 
			theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)) +
			labs(x="Standard Watershed Posterior", y="Alt Watershed Posterior)",colour="") + 
			scale_color_manual(values=c("#7F5A83", "#0D324D", "#BFCDE0")) +
			gtex_v8_figure_theme()
	return(plotter)

}

make_posterior_scatter_colored_by_outlier_class <- function(posterior, alt_posterior, outlier_status, title) {
	df <- data.frame(standard_posterior=posterior, alt_posterior=alt_posterior, outlier=factor(outlier_status))
	plotter <- ggplot(df, aes(x=standard_posterior, y=alt_posterior, colour=outlier)) + geom_point(size=.8) +
			geom_abline() + 
			theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)) +
			labs(x="Standard Watershed Posterior", y="Alt Watershed Posterior)",colour="",title=title) + 
			gtex_v8_figure_theme() +
			theme(legend.position="bottom")
	return(plotter)
}

compare_watershed_posteriors_seperated_by_class_with_different_training_inputs <- function(roc_3_class_data_input, alt_roc_3_class_data_input, roc_object_exact, alt_roc_object_exact) {
	roc_3_class_sample_names <- extract_watershed_sample_names(roc_3_class_data_input)
	alt_roc_3_class_sample_names <- extract_watershed_sample_names(alt_roc_3_class_data_input)
	indices <- alt_roc_3_class_sample_names %in% roc_3_class_sample_names
	if (sum(roc_3_class_sample_names != alt_roc_3_class_sample_names[indices]) != 0) {
		print("FUNDAMENTAL ASSUMPTION ERROR in compare_watershed_posteriors_with_different_training_inputs")
	}

	roc_posteriors <- roc_object_exact$model_params$posterior
	alt_roc_posteriors <- alt_roc_object_exact$model_params$posterior[indices,]
	
	outliers <- extract_sample_outlier_values(roc_3_class_data_input)
	alt_outliers <- extract_sample_outlier_values(alt_roc_3_class_data_input)
	alt_outliers <- alt_outliers[indices,]
	if (sum(outliers!=alt_outliers) != 0) {
		print("FUNDAMENTAL assumption error")
	}

	splicing_scatter_plot <- make_posterior_scatter_colored_by_outlier_class(roc_posteriors[,1], alt_roc_posteriors[,1], outliers[,1], "Splicing")
	expression_scatter_plot <- make_posterior_scatter_colored_by_outlier_class(roc_posteriors[,2], alt_roc_posteriors[,2], outliers[,2], "Expression")
	ase_scatter_plot <- make_posterior_scatter_colored_by_outlier_class(roc_posteriors[,3], alt_roc_posteriors[,3], outliers[,3], "ASE")

	combined_scatter_plots <- plot_grid(ase_scatter_plot, splicing_scatter_plot, expression_scatter_plot, ncol=3)

	return(combined_scatter_plots)

}

plot_three_class_auprc_bootstrap_distributions <- function(watershed_roc, river_roc, number_of_dimensions) {
	auprc <- c()
	outlier_type <- c()
	prediction_type <- c()

	for (dimension in 1:number_of_dimensions) {
		# Extract auprc data on this dimension
		watershed_roc_object <- watershed_roc[[dimension]]
		river_roc_object_ind <- river_roc[[dimension]]

		dimension_name <- watershed_roc_object$name

		# Tied watershed (exact inference)
		auprc <- c(auprc, watershed_roc_object$evaROC$watershed_pr_auc_bootstraps)
		outlier_type <- c(outlier_type, rep(dimension_name, length(watershed_roc_object$evaROC$watershed_pr_auc_bootstraps)))
		prediction_type <- c(prediction_type, rep("Watershed", length(watershed_roc_object$evaROC$watershed_pr_auc_bootstraps)))
		print(paste0(dimension_name, " Watershed"))
		full_pr_auc = watershed_roc_object$evaROC$watershed_pr_auc
		bootstrap_pr_aucs = watershed_roc_object$evaROC$watershed_pr_auc_bootstraps
		c_u = full_pr_auc - quantile(bootstrap_pr_aucs-full_pr_auc, .025)
		c_l = full_pr_auc - quantile(bootstrap_pr_aucs-full_pr_auc, .975)
		print(paste0("[", c_l, " , ", c_u, "]"))

		# RIVER
		auprc <- c(auprc, river_roc_object_ind$evaROC$watershed_pr_auc_bootstraps)
		outlier_type <- c(outlier_type, rep(dimension_name, length(river_roc_object_ind$evaROC$watershed_pr_auc_bootstraps)))
		prediction_type <- c(prediction_type, rep("RIVER", length(river_roc_object_ind$evaROC$watershed_pr_auc_bootstraps)))
		print(paste0(dimension_name, " RIVER"))
		full_pr_auc = river_roc_object_ind$evaROC$watershed_pr_auc
		bootstrap_pr_aucs = river_roc_object_ind$evaROC$watershed_pr_auc_bootstraps
		c_u = full_pr_auc - quantile(bootstrap_pr_aucs-full_pr_auc, .025)
		c_l = full_pr_auc - quantile(bootstrap_pr_aucs-full_pr_auc, .975)
		print(paste0("[", c_l, " , ", c_u, "]"))

		# GAM
		auprc <- c(auprc, river_roc_object_ind$evaROC$GAM_pr_auc_bootstraps)
		outlier_type <- c(outlier_type, rep(dimension_name, length(river_roc_object_ind$evaROC$GAM_pr_auc_bootstraps)))
		prediction_type <- c(prediction_type, rep("GAM", length(river_roc_object_ind$evaROC$GAM_pr_auc_bootstraps)))
		print(paste0(dimension_name, " GAM"))
		full_pr_auc = river_roc_object_ind$evaROC$GAM_pr_auc
		bootstrap_pr_aucs = river_roc_object_ind$evaROC$GAM_pr_auc_bootstraps
		c_u = full_pr_auc - quantile(bootstrap_pr_aucs-full_pr_auc, .025)
		c_l = full_pr_auc - quantile(bootstrap_pr_aucs-full_pr_auc, .975)
		print(paste0("[", c_l, " , ", c_u, "]"))
	}

	# Put into compact df
	df <- data.frame(auprc=auprc, model=factor(prediction_type,levels=c("Watershed","RIVER","GAM")), outlier_type=factor(outlier_type))

	outlier_type <- "ase"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=auprc, fill=model, colour=model)) + geom_histogram(position="identity",alpha=0.5) + 
                labs(x="AUC(PR)", y="Non-parametric bootstrap samples", colour="", fill="", title="ASE") +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_fill_manual(values=c("steelblue3", "white", "firebrick4")) +
                gtex_v8_figure_theme()
	outlier_type <- "splicing"
  	plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=auprc, fill=model, colour=model)) + geom_histogram(position="identity",alpha=0.5) + 
                labs(x="AUC(PR)", y="Non-parametric bootstrap samples", colour="", fill="", title="Splicing") +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_fill_manual(values=c("steelblue3", "white", "firebrick4")) +
                gtex_v8_figure_theme()

	outlier_type <- "total_expression"
  	plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=auprc, fill=model, colour=model)) + geom_histogram(position="identity",alpha=0.5) + 
                labs(x="AUC(PR)", y="Non-parametric bootstrap samples", colour="", fill="", title="Expression") +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_fill_manual(values=c("steelblue3", "white", "firebrick4")) +
                gtex_v8_figure_theme()


    legend <- get_legend(plotter_ase)
    combined_plots <- plot_grid(plotter_ase + theme(legend.position="none"),plotter_splice+ theme(legend.position="none"), plotter_te + theme(legend.position="none"), ncol=1)


	return(plot_grid(combined_plots,legend,ncol=1,rel_heights=c(1,.06)))

}

plot_three_class_delta_auprc_bootstrap_distributions <- function(watershed_roc, river_roc, number_of_dimensions) {
	delta_auprc <- c()
	outlier_type <- c()
	ci <- c()
	pvalues <- c()
	for (dimension in 1:number_of_dimensions) {
		# Extract auprc data on this dimension
		watershed_roc_object <- watershed_roc[[dimension]]
		river_roc_object_ind <- river_roc[[dimension]]

		dimension_name <- watershed_roc_object$name

		

		population_watershed_pr_auc = watershed_roc_object$evaROC$watershed_pr_auc
		population_river_pr_auc = river_roc_object_ind$evaROC$watershed_pr_auc
		population_delta_pr_auc = population_watershed_pr_auc - population_river_pr_auc


		# BootStrapped delta_auprc
		bootstrapped_delta_pr_auc = watershed_roc_object$evaROC$watershed_pr_auc_bootstraps - river_roc_object_ind$evaROC$watershed_pr_auc_bootstraps

		c_u = population_delta_pr_auc - quantile(bootstrapped_delta_pr_auc-population_delta_pr_auc, .025)
		c_l = population_delta_pr_auc - quantile(bootstrapped_delta_pr_auc-population_delta_pr_auc, .975)

		ci <- c(ci, paste0("[", c_l, " , ", c_u, "]"))

		fraction = sum(watershed_roc_object$evaROC$watershed_pr_auc_bootstraps > river_roc_object_ind$evaROC$watershed_pr_auc_bootstraps)/length(river_roc_object_ind$evaROC$watershed_pr_auc_bootstraps)
		print(paste0(dimension_name, ": ", 1-fraction))

		pvalues <- bootstrapped_delta_pr_auc
		delta_auprc <- c(delta_auprc, bootstrapped_delta_pr_auc)
		outlier_type <- c(outlier_type, rep(dimension_name, length(bootstrapped_delta_pr_auc)))

	}

	# Put into compact df
	df <- data.frame(delta_auprc=delta_auprc, outlier_type=factor(outlier_type))

	outlier_type <- "ase"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=delta_auprc)) + geom_histogram(position="identity", color="darkblue", fill="lightblue") + 
                labs(x="Watershed AUC(PR) - RIVER AUC(PR)", y="Non-parametric bootstrap samples", title=paste0("ASE (95% CI", ci[3], ")")) +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_vline(xintercept = 0.0, color='red') +
                gtex_v8_figure_theme()
	outlier_type <- "splicing"
  	plotter_splicing <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=delta_auprc)) + geom_histogram(position="identity", color="darkblue", fill="lightblue") + 
                labs(x="Watershed AUC(PR) - RIVER AUC(PR)", y="Non-parametric bootstrap samples", title=paste0("Splicing (95% CI", ci[1], ")")) +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_vline(xintercept = 0.0, color='red') +
                gtex_v8_figure_theme()

	outlier_type <- "total_expression"
  	plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=delta_auprc)) + geom_histogram(position="identity", color="darkblue", fill="lightblue") + 
                labs(x="Watershed AUC(PR) - RIVER AUC(PR)", y="Non-parametric bootstrap samples", title=paste0("Expression (95% CI", ci[2], ")")) +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_vline(xintercept = 0.0, color='red') +
                gtex_v8_figure_theme()


    combined_plots <- plot_grid(plotter_ase + theme(legend.position="none"),plotter_splicing+ theme(legend.position="none"), plotter_te + theme(legend.position="none"), ncol=1)


	return(combined_plots)
}


auprc_bootstrap_hypothesis_testing <- function(watershed_object, river_object, number_of_dimensions) {
	# Number of bootstap iterations
	B = 20000
	# Model predictions
	watershed_predictions <- watershed_object$watershed_predictions
	river_predictions <- river_object$watershed_predictions
	# Held out labels
	labels = watershed_object$held_out_labels

	num_instances = dim(labels)[1]


	auprc <- c()
	outlier_type <- c()
	pvalues <- c()
	population_delta_auprc_arr <- c()

	for (dimension in 1:number_of_dimensions) {
		# Get dimension name
		dimension_name <- watershed_object$roc[[dimension]]$name
		dimension_watershed_predictions <- watershed_predictions[,dimension]
		dimension_river_predictions <- river_predictions[,dimension]
		dimension_labels <- labels[,dimension]

		merged_predictions <- c(dimension_watershed_predictions, dimension_river_predictions)
		merged_labels <- c(dimension_labels, dimension_labels)

		# Compute delta AUPRC at population level
		watershed_population_auprc <- pr.curve(scores.class0=dimension_watershed_predictions[dimension_labels==1], scores.class1=dimension_watershed_predictions[dimension_labels==0])$auc.integral
		river_population_auprc <- pr.curve(scores.class0=dimension_river_predictions[dimension_labels==1], scores.class1=dimension_river_predictions[dimension_labels==0])$auc.integral
		delta_population_auprc <- watershed_population_auprc - river_population_auprc
		population_delta_auprc_arr <- c(population_delta_auprc_arr, delta_population_auprc)

		# Create emperical distribution on delta auprcs
		emperical_delta_auprc <- c()
		for (bootstrap_iteration in 1:B) {
			bootstrap_indices = sample(1:(2*num_instances),size=2*num_instances,replace=TRUE)
			bootstrap_predictions = merged_predictions[bootstrap_indices]
			bootstrap_labels = merged_labels[bootstrap_indices]
			
			bootstrap_watershed_predictions = bootstrap_predictions[1:num_instances]
			bootstrap_watershed_labels = bootstrap_labels[1:num_instances]
			bootstrap_watershed_auprc = pr.curve(scores.class0=bootstrap_watershed_predictions[bootstrap_watershed_labels==1], scores.class1=bootstrap_watershed_predictions[bootstrap_watershed_labels==0])$auc.integral

			bootstrap_river_predictions = bootstrap_predictions[(num_instances+1):length(bootstrap_predictions)]
			bootstrap_river_labels = bootstrap_labels[(num_instances+1):length(bootstrap_labels)]
			bootstrap_river_auprc = pr.curve(scores.class0=bootstrap_river_predictions[bootstrap_river_labels==1], scores.class1=bootstrap_river_predictions[bootstrap_river_labels==0])$auc.integral
			

			bootstrap_delta_auprc <- bootstrap_watershed_auprc - bootstrap_river_auprc

			emperical_delta_auprc <- c(emperical_delta_auprc, bootstrap_delta_auprc)
			#bootstrap_auprc <- pr.curve(scores.class0=bootstrap_predictions[bootstrap_labels==1], scores.class1=bootstrap_predictions[bootstrap_labels==0])$auc.integral
		}

		pvalue <- sum(emperical_delta_auprc > delta_population_auprc)/length(emperical_delta_auprc)
		pvalues <- c(pvalues, pvalue)

		auprc <- c(auprc, emperical_delta_auprc)
		outlier_type <- c(outlier_type, rep(dimension_name, length(emperical_delta_auprc)))
	}

	df <- data.frame(delta_auprc=auprc, outlier_type=factor(outlier_type))

	outlier_type <- "ase"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=delta_auprc)) + geom_histogram(position="identity", color="darkblue", fill="lightblue") + 
                labs(x="Watershed AUC(PR) - RIVER AUC(PR)", y="Non-parametric bootstrap samples", title=paste0("ASE (pvalue=", pvalues[3], ")")) +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_vline(xintercept = population_delta_auprc_arr[3], color='red') +
                gtex_v8_figure_theme()
	
	outlier_type <- "splicing"
  	plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=delta_auprc)) + geom_histogram(position="identity", color="darkblue", fill="lightblue") + 
                labs(x="Watershed AUC(PR) - RIVER AUC(PR)", y="Non-parametric bootstrap samples", title=paste0("Splicing (pvalue=", pvalues[1], ")")) +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_vline(xintercept = population_delta_auprc_arr[1], color='red') +
                gtex_v8_figure_theme()

	outlier_type <- "total_expression"
  	plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=delta_auprc)) + geom_histogram(position="identity", color="darkblue", fill="lightblue") + 
                labs(x="Watershed AUC(PR) - RIVER AUC(PR)", y="Non-parametric bootstrap samples", title=paste0("Expression (pvalue=", pvalues[2], ")")) +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_vline(xintercept = population_delta_auprc_arr[2], color='red') +
                gtex_v8_figure_theme()
	
	combined_plots <- plot_grid(plotter_ase + theme(legend.position="none"),plotter_splice+ theme(legend.position="none"), plotter_te + theme(legend.position="none"), ncol=1)

    return(combined_plots)
}

delta_auprc_bootstrap_hypothesis_testing <- function(watershed_object, river_object, number_of_dimensions) {
	# Number of bootstap iterations
	set.seed(5)
	B = 20000
	# Model predictions
	watershed_predictions <- watershed_object$watershed_predictions
	river_predictions <- river_object$watershed_predictions
	# Held out labels
	labels = watershed_object$held_out_labels

	num_instances = dim(labels)[1]


	auprc <- c()
	outlier_type <- c()
	pvalues <- c()
	population_delta_auprc_arr <- c()

	for (dimension in 1:number_of_dimensions) {
		# Get dimension name
		dimension_name <- watershed_object$roc[[dimension]]$name
		dimension_watershed_predictions <- watershed_predictions[,dimension]
		dimension_river_predictions <- river_predictions[,dimension]
		dimension_labels <- labels[,dimension]

		# Compute delta AUPRC at population level
		watershed_population_auprc <- pr.curve(scores.class0=dimension_watershed_predictions[dimension_labels==1], scores.class1=dimension_watershed_predictions[dimension_labels==0])$auc.integral
		river_population_auprc <- pr.curve(scores.class0=dimension_river_predictions[dimension_labels==1], scores.class1=dimension_river_predictions[dimension_labels==0])$auc.integral
		delta_population_auprc <- watershed_population_auprc - river_population_auprc
		population_delta_auprc_arr <- c(population_delta_auprc_arr, delta_population_auprc)

		# Create emperical distribution on delta auprcs
		emperical_delta_auprc <- c()
		for (bootstrap_iteration in 1:B) {
			bootstrap_indices = sample(1:num_instances,size=num_instances,replace=TRUE)

			bootstrap_watershed_predictions = dimension_watershed_predictions[bootstrap_indices]
			bootstrap_river_predictions = dimension_river_predictions[bootstrap_indices]
			bootstrap_labels = dimension_labels[bootstrap_indices]

			bootstrap_null_watershed_predictions <- bootstrap_watershed_predictions
			bootstrap_null_river_predictions <- bootstrap_watershed_predictions

			x1 <- rbinom(n=num_instances, size=1, prob=.5)
			x2 <- rbinom(n=num_instances, size=1, prob=.5)

			bootstrap_null_watershed_predictions[x1==1] <- bootstrap_river_predictions[x1==1]
			bootstrap_null_river_predictions[x2==1] <- bootstrap_river_predictions[x2==1]


			bootstrap_watershed_auprc = pr.curve(scores.class0=bootstrap_null_watershed_predictions[bootstrap_labels==1], scores.class1=bootstrap_null_watershed_predictions[bootstrap_labels==0])$auc.integral
			bootstrap_river_auprc = pr.curve(scores.class0=bootstrap_null_river_predictions[bootstrap_labels==1], scores.class1=bootstrap_null_river_predictions[bootstrap_labels==0])$auc.integral

			bootstrap_delta_auprc <- bootstrap_watershed_auprc - bootstrap_river_auprc

			emperical_delta_auprc <- c(emperical_delta_auprc, bootstrap_delta_auprc)
			#bootstrap_auprc <- pr.curve(scores.class0=bootstrap_predictions[bootstrap_labels==1], scores.class1=bootstrap_predictions[bootstrap_labels==0])$auc.integral
		}

		pvalue <- sum(emperical_delta_auprc > delta_population_auprc)/length(emperical_delta_auprc)
		print(dimension_name)
		print(pvalue)
		pvalues <- c(pvalues, pvalue)

		auprc <- c(auprc, emperical_delta_auprc)
		outlier_type <- c(outlier_type, rep(dimension_name, length(emperical_delta_auprc)))
	}

	df <- data.frame(delta_auprc=auprc, outlier_type=factor(outlier_type))

	outlier_type <- "ase"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=delta_auprc)) + geom_histogram(position="identity", color="darkblue", fill="lightblue") + 
                labs(x="Watershed AUC(PR) - RIVER AUC(PR)", y="Non-parametric bootstrap samples", title=paste0("ASE (pvalue=", pvalues[3], ")")) +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_vline(xintercept = population_delta_auprc_arr[3], color='red') +
                gtex_v8_figure_theme()
	
	outlier_type <- "splicing"
  	plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=delta_auprc)) + geom_histogram(position="identity", color="darkblue", fill="lightblue") + 
                labs(x="Watershed AUC(PR) - RIVER AUC(PR)", y="Non-parametric bootstrap samples", title=paste0("Splicing (pvalue=", pvalues[1], ")")) +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_vline(xintercept = population_delta_auprc_arr[1], color='red') +
                gtex_v8_figure_theme()

	outlier_type <- "total_expression"
  	plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=delta_auprc)) + geom_histogram(position="identity", color="darkblue", fill="lightblue") + 
                labs(x="Watershed AUC(PR) - RIVER AUC(PR)", y="Non-parametric bootstrap samples", title=paste0("Expression (pvalue=", pvalues[2], ")")) +
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                geom_vline(xintercept = population_delta_auprc_arr[2], color='red') +
                gtex_v8_figure_theme()
	
	combined_plots <- plot_grid(plotter_ase + theme(legend.position="none"),plotter_splice+ theme(legend.position="none"), plotter_te + theme(legend.position="none"), ncol=1)

    return(combined_plots)
}


options(bitmapType = 'cairo', device = 'pdf')
############################
# Command line args
############################
three_class_roc_dir <- args[1]
tbt_roc_dir <- args[2]
genomic_annotation_names_file <- args[3]
tissue_names_file <- args[4]
chrom_hmm_to_tissue_mapping_file <- args[5]
output_dir <- args[6]
tissue_colors_file <- args[7]


# Read in tissue colors and names
tissue_colors = read.table(tissue_colors_file, header = T, stringsAsFactors = F, sep = "\t")
# slight mislabeling
for (tiss_num in 1:length(tissue_colors$tissue_id)) {
	if (tissue_colors$tissue_id[tiss_num] == "Brain_Spinal_cord_cervical_c1") {
		tissue_colors$tissue_id[tiss_num] = "Brain_Spinal_cord_cervical_c.1"
	}
	if (tissue_colors$tissue_id[tiss_num] == "Cells_EBVtransformed_lymphocytes") {
		tissue_colors$tissue_id[tiss_num] = "Cells_EBV.transformed_lymphocytes"
	}
}


#############################
# Load in 3 outlier type data
############################
pseudocount <- 30
input_stem <- paste0(three_class_roc_dir, "fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_",pseudocount,"_seed_3")

####### input data
roc_3_class_data_input <- readRDS(paste0(input_stem, "_data_input.rds"))

####### Exact watershed
independent_variables = "false"
inference_method = "exact"
output_root <- paste0(input_stem,"_inference_", inference_method, "_independent_", independent_variables)
roc_object_exact <- readRDS(paste0(output_root, "_roc_object3.rds"))


####### Pseudolikelihood approximation to watershed
independent_variables = "false"
inference_method = "pseudolikelihood"
output_root <- paste0(input_stem,"_inference_", inference_method, "_independent_", independent_variables)
roc_object_pseudo <- readRDS(paste0(output_root, "_roc_object3.rds"))


####### Exact RIVER
independent_variables = "true"
inference_method = "exact"
output_root <- paste0(input_stem,"_inference_", inference_method, "_independent_", independent_variables)
roc_object_independent <- readRDS(paste0(output_root, "_roc_object3.rds"))



#######################################
## delta AUPRC bootstrap hypothesis testing
#######################################
if (FALSE) {
output_file <- paste0(output_dir, "watershed_3_class_delta_auprc_hypothesis_testing_bootstrap_distributions.pdf")
delta_auprc_bootstrap_hypothesis_testing_histogram <- delta_auprc_bootstrap_hypothesis_testing(roc_object_exact, roc_object_independent, 3)
ggsave(delta_auprc_bootstrap_hypothesis_testing_histogram, file=output_file, width=7.2, height=9, units="in")


#######################################
## AUPRC bootstrap hypothesis testing
#######################################
output_file <- paste0(output_dir, "watershed_3_class_auprc_hypothesis_testing_bootstrap_distributions.pdf")
auprc_bootstrap_hypothesis_testing_histogram <- auprc_bootstrap_hypothesis_testing(roc_object_exact, roc_object_independent, 3)
ggsave(auprc_bootstrap_hypothesis_testing_histogram, file=output_file, width=7.2, height=9, units="in")

#######################################
## Visualize bootstrap delta auprc distributions
#######################################
output_file <- paste0(output_dir, "watershed_3_class_delta_auprc_bootstrap_distributions.pdf")
roc_3_bootstrap_distributions <- plot_three_class_delta_auprc_bootstrap_distributions(roc_object_exact$roc, roc_object_independent$roc, 3)
ggsave(roc_3_bootstrap_distributions, file=output_file, width=7.2, height=9, units="in")

#######################################
## Visualize bootstrap auprc distributions
#######################################
output_file <- paste0(output_dir, "watershed_3_class_auprc_bootstrap_distributions.pdf")
roc_3_bootstrap_distributions <- plot_three_class_auprc_bootstrap_distributions(roc_object_exact$roc, roc_object_independent$roc, 3)
ggsave(roc_3_bootstrap_distributions, file=output_file, width=7.2, height=9, units="in")
}


#######################################
## Visualize theta pair terms for exact inference
#######################################
output_file <- paste0(output_dir, "watershed_3_class_roc_exact_inference_theta_pair_heatmap.pdf")
roc_3_theta_pair_exact <- plot_three_class_theta_pair_term(roc_object_exact$model_params$theta_pair)
ggsave(roc_3_theta_pair_exact, file=output_file, width=7.2, height=4, units="in")


#######################################
## Make watershed 3-class absolute risk plot
#######################################
pvalue= 0.0027
output_file <- paste0(output_dir, "watershed_gam_3_class_roc_absolute_risk.pdf")
roc_3_absolute_risk_bar_plot <- absolute_risk_plot(roc_object_exact$gam_predictions, roc_object_exact$watershed_predictions, roc_3_class_data_input, pvalue)
ggsave(roc_3_absolute_risk_bar_plot, file=output_file, width=7.2, height=4.6, units="in")


#######################################
## Make watershed 3-class absolute risk plot (including cadding)
#######################################
pvalue= 0.0027
output_file <- paste0(output_dir, "watershed_gam_cadd_3_class_roc_absolute_risk.pdf")
roc_3_absolute_risk_bar_plot_w_cadd <- absolute_risk_plot_with_cadd(roc_object_exact$gam_predictions, roc_object_exact$watershed_predictions, roc_object_exact$cadd_predictions, roc_3_class_data_input, pvalue)
ggsave(roc_3_absolute_risk_bar_plot_w_cadd, file=output_file, width=7.2, height=6, units="in")


#######################################
## Visualize precision-recall curves for river, GAM, watershed-exact comparison and all three outlier types (te, splice, ase)
#######################################
number_of_dimensions <- 3
output_file <- paste0(output_dir, "watershed_3_class_roc_exact_gam_river_watershed_precision_recall.pdf")
gam_river_watershed_3_class_pr_curves <- plot_pr_gam_river_watershed_comparison_curve(roc_object_exact$roc, roc_object_independent$roc, number_of_dimensions)
ggsave(gam_river_watershed_3_class_pr_curves, file=output_file, width=10, height=3.0, units="in")

#######################################
## Visualize precision-recall curves for CADD and  watershed-exact comparison and all three outlier types (te, splice, ase)
#######################################
number_of_dimensions <- 3
output_file <- paste0(output_dir, "watershed_3_class_roc_exact_cadd_watershed_precision_recall.pdf")
cadd_watershed_3_class_pr_curves <- plot_pr_cadd_watershed_comparison_curve(roc_object_exact$roc, roc_object_independent$roc, number_of_dimensions)
ggsave(cadd_watershed_3_class_pr_curves, file=output_file, width=10, height=3.0, units="in")


#######################################
## Visualize precision-recall curves for river, watershed-vi, watershed-exact comparison and all three outlier types (te, splice, ase)
#######################################
output_file <- paste0(output_dir, "watershed_3_class_roc_river_exact_vi_watershed_precision_recall.pdf")
gam_river_watershed_3_class_pr_curves_vi_cmp <- plot_pr_river_watershed_comparison_curve(roc_object_exact$roc, roc_object_pseudo$roc, roc_object_independent$roc, number_of_dimensions)
#ggsave(gam_river_watershed_3_class_pr_curves_vi_cmp, file=output_file, width=7.2, height=3.0, units="in")


#######################################
## Visualize beta differences for exact vs inference
#######################################
output_file <- paste0(output_dir, "watershed_exact_vs_approximate_beta_estimate_scatter.pdf")
beta_comparison_scatter <- plot_beta_difference_scatter_between_exact_and_vi(roc_object_exact$model_params, roc_object_pseudo$model_params)
#ggsave(beta_comparison_scatter, file=output_file, width=7.2, height=3.0, units="in")


#######################################
## Combined supplementary VI plot for 3class model
#######################################
output_file <- paste0(output_dir, "watershed_3_class_exact_vi_comparison_supplementary_figure.pdf")
combined <- plot_grid(beta_comparison_scatter, gam_river_watershed_3_class_pr_curves_vi_cmp, ncol=1,labels=c("A","B"))
ggsave(combined, file=output_file, width=7.2, height=5.0, units="in")


#######################################
## Visualize Confusion matrix for both RIVER and Watershed (exact and vi)
#######################################
output_file <- paste0(output_dir, "watershed_river_confusion_matrix_heatmap.pdf")
confusion_heatmap <- visualize_river_and_watershed_confusion_matrices(roc_object_exact$confusion, roc_object_pseudo$confusion, roc_object_independent$confusion)
ggsave(confusion_heatmap, file=output_file, width=7.2, height=7.0, units="in")



#######################################
## Compare Watershed to versions of Watershed trained on alternative training data sets
#######################################
gene_pval="0.05"
output_file <- paste0(output_dir, "compare_watershed_posterior_with_different_training_inputs_01_05_scatter.pdf")
alt_roc_3_class_data_input <- readRDS(paste0(three_class_roc_dir, "fully_observed_te_ase_splicing_outliers_gene_pvalue_", gene_pval, "_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_gene_theshold_comparison_data_input.rds"))
alt_roc_object_exact <- readRDS(paste0(three_class_roc_dir, "fully_observed_te_ase_splicing_outliers_gene_pvalue_", gene_pval, "_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_gene_theshold_comparison_inference_exact_independent_false_roc_object3.rds"))
scatter <- compare_watershed_posteriors_with_different_training_inputs(roc_3_class_data_input, alt_roc_3_class_data_input, roc_object_exact, alt_roc_object_exact)
ggsave(scatter, file=output_file, width=7.2, height=5.0, units="in")

gene_pval="0.05"
output_file <- paste0(output_dir, "compare_watershed_posterior_seperated_by_outlier_class_with_different_training_inputs_01_05_scatter.pdf")
alt_roc_3_class_data_input <- readRDS(paste0(three_class_roc_dir, "fully_observed_te_ase_splicing_outliers_gene_pvalue_", gene_pval, "_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_gene_theshold_comparison_data_input.rds"))
alt_roc_object_exact <- readRDS(paste0(three_class_roc_dir, "fully_observed_te_ase_splicing_outliers_gene_pvalue_", gene_pval, "_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_gene_theshold_comparison_inference_exact_independent_false_roc_object3.rds"))
scatter <- compare_watershed_posteriors_seperated_by_class_with_different_training_inputs(roc_3_class_data_input, alt_roc_3_class_data_input, roc_object_exact, alt_roc_object_exact)
ggsave(scatter, file=output_file, width=13.2, height=7.0, units="in")

#######################################
## Visualize precision-recall curves for river, GAM, watershed-exact comparison and all three outlier types (te, splice, ase)
#######################################
number_of_dimensions <- 3
output_file <- paste0(output_dir, "compare_gene_threshold_watershed_3_class_roc_exact_gam_river_watershed_precision_recall.pdf")
alt_roc_object_independent <- readRDS(paste0(three_class_roc_dir, "fully_observed_te_ase_splicing_outliers_gene_pvalue_", gene_pval, "_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_gene_theshold_comparison_inference_exact_independent_true_roc_object3.rds"))
gam_river_watershed_3_class_pr_curves <- plot_pr_gam_river_watershed_comparison_curve(alt_roc_object_exact$roc, alt_roc_object_independent$roc, number_of_dimensions)
ggsave(gam_river_watershed_3_class_pr_curves, file=output_file, width=10, height=3.0, units="in")

print("DONE")
############################
# Model hyperparameters
############################
pseudocount <- 10
phi_update_method <- "fixed"
number_of_dimensions <- 49


# Load in names of genomic annotations
anno_names <- as.character(read.table(genomic_annotation_names_file, header=FALSE)[,1])

# Load in tissue_names
tissue_names <- as.character(read.table(tissue_names_file)$V1)

# Load in data mapping tissue name to chromHMM cell type
chrom_hmm_to_tissue <- read.table(chrom_hmm_to_tissue_mapping_file, header=TRUE, sep="\t")

############################
# Load in splicing data
############################
outlier_type <- "splicing"
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
splicing_watershed_obj <- readRDS(paste0(tbt_roc_dir, stem, "_inference_pseudolikelihood_independent_false_roc_object2.rds"))
splicing_river_obj <- readRDS(paste0(tbt_roc_dir, stem,"_inference_exact_independent_true_roc_object2.rds"))
splicing_tissue_names <- get_tissue_names(splicing_watershed_obj, number_of_dimensions)

############################
# Load in TE data
############################
outlier_type <- "total_expression"
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
te_watershed_obj <- readRDS(paste0(tbt_roc_dir, stem, "_inference_pseudolikelihood_independent_false_roc_object2.rds"))
te_river_obj <- readRDS(paste0(tbt_roc_dir, stem,"_inference_exact_independent_true_roc_object2.rds"))
te_tissue_names <- get_tissue_names(te_watershed_obj, number_of_dimensions)
############################
# Load in ASE data
############################
outlier_type <- "ase"
number_of_dimensions <- 49
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
ase_watershed_obj <- readRDS(paste0(tbt_roc_dir, stem, "_inference_pseudolikelihood_independent_false_roc_object2.rds"))
ase_river_obj <- readRDS(paste0(tbt_roc_dir, stem,"_inference_exact_independent_true_roc_object2.rds"))
ase_tissue_names <- get_tissue_names(ase_watershed_obj, number_of_dimensions)




tissues_with_suffiecient_n2_pairs <- get_tissues_with_many_n2_pairs(splicing_watershed_obj$roc, te_watershed_obj$roc, ase_watershed_obj$roc, splicing_tissue_names, te_tissue_names, ase_tissue_names, 5)

######################################
# Visualize tissue specific thetas
######################################
#te_tissue_specific_theta_plot <- make_tissue_specific_theta(te_watershed_obj$model_params$theta, anno_names, tissue_names, chrom_hmm_to_tissue, "prom")

#te_tissue_specific_theta_plot <- make_tissue_specific_theta(te_watershed_obj$model_params$theta, anno_names, tissue_names, chrom_hmm_to_tissue, "enh")



###################################################
# Visualize splicing results
###################################################
outlier_type <- "splicing"
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
number_of_dimensions <- 49
######################################
# Visualize theta-pair heatmap to look at correlation structure across tissues
######################################
output_file <- paste0(output_dir, stem, "_tissue_by_tissue_theta_pair_heatmap.pdf")
splicing_theta_pair_heatmap <- make_theta_pair_heatmap(splicing_watershed_obj$model_params$theta_pair, number_of_dimensions, splicing_tissue_names, "Splicing")
ggsave(splicing_theta_pair_heatmap, file=output_file, width=7.2, height=4.6, units="in")

######################################
# Visualize P(E|Z)
######################################
output_file <- paste0(output_dir, stem, "_tissue_by_tissue_phi.pdf")
watershed_phi_plot <- make_phi_bar_plot(splicing_watershed_obj$model_params$phi, splicing_tissue_names, paste0("watershed-",outlier_type))
river_phi_plot <- make_phi_bar_plot(splicing_river_obj$model_params$phi, splicing_tissue_names, paste0("river-",outlier_type))
phi_plot <- plot_grid(river_phi_plot, watershed_phi_plot, ncol=2)
ggsave(phi_plot, file=output_file, width=8.0, height=4.6, units="in")


######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(output_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
splicing_tbt_auc_lolipop_plot <- make_tbt_auc_lolipop_plot(splicing_watershed_obj$roc, splicing_river_obj$roc, splicing_tissue_names, "Splicing", tissues_with_suffiecient_n2_pairs)
# ggsave(splicing_tbt_auc_lolipop_plot, file=output_file, width=7.2, height=4.0, units="in")

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(output_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
splicing_tbt_auc_lolipop_plot_v2 <- make_tbt_auc_lolipop_plot_v2(splicing_watershed_obj$roc, splicing_river_obj$roc, splicing_tissue_names, "Splicing", tissues_with_suffiecient_n2_pairs)
# ggsave(splicing_tbt_auc_lolipop_plot, file=output_file, width=7.2, height=4.0, units="in")



######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(output_dir, stem, "tissue_by_tissue_pr_auc_between_median_river_watershed_lolipop.pdf")
splicing_tbt_auc_lolipop_plot_with_median_river <- make_tbt_auc_lolipop_plot_with_median_river(splicing_watershed_obj$roc, splicing_tissue_names, "Splicing", tissues_with_suffiecient_n2_pairs)
# ggsave(splicing_tbt_auc_lolipop_plot_with_median_river, file=output_file, width=7.2, height=4.0, units="in")













###################################################
# Visualize total expression results
###################################################
outlier_type <- "total_expression"
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
number_of_dimensions <- 49
######################################
# Visualize theta-pair heatmap to look at correlation structure across tissues
######################################
output_file <- paste0(output_dir, stem, "_tissue_by_tissue_theta_pair_heatmap.pdf")
te_theta_pair_heatmap <- make_theta_pair_heatmap(te_watershed_obj$model_params$theta_pair, number_of_dimensions, te_tissue_names, "Expression")
ggsave(te_theta_pair_heatmap, file=output_file, width=7.2, height=4.6, units="in")

######################################
# Visualize theta-pair heatmap to look at correlation structure across tissues
######################################
output_file <- paste0(output_dir, stem, "_tissue_by_tissue_fig5_theta_pair_heatmap.pdf")
fig_5_te_theta_pair_heatmap <- make_fig_5_theta_pair_heatmap(te_watershed_obj$model_params$theta_pair, number_of_dimensions, te_tissue_names, tissue_colors)
ggsave(fig_5_te_theta_pair_heatmap, file=output_file, width=7.2, height=4.6, units="in")
######################################
# Visualize P(E|Z)
######################################
output_file <- paste0(output_dir, stem, "_tissue_by_tissue_phi.pdf")
watershed_phi_plot <- make_phi_bar_plot(te_watershed_obj$model_params$phi, te_tissue_names, paste0("watershed-",outlier_type))
river_phi_plot <- make_phi_bar_plot(te_river_obj$model_params$phi, te_tissue_names, paste0("river-",outlier_type))
phi_plot <- plot_grid(river_phi_plot, watershed_phi_plot, ncol=2)
ggsave(phi_plot, file=output_file, width=8.0, height=4.6, units="in")


######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(output_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
te_tbt_auc_lolipop_plot <- make_tbt_auc_lolipop_plot(te_watershed_obj$roc, te_river_obj$roc, te_tissue_names, "Expression", tissues_with_suffiecient_n2_pairs)
# ggsave(te_tbt_auc_lolipop_plot, file=output_file, width=7.2, height=4.0, units="in")

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(output_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
te_tbt_auc_lolipop_plot_v2 <- make_tbt_auc_lolipop_plot_v2(te_watershed_obj$roc, te_river_obj$roc, te_tissue_names, "Expression", tissues_with_suffiecient_n2_pairs)
# ggsave(splicing_tbt_auc_lolipop_plot, file=output_file, width=7.2, height=4.0, units="in")


######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(output_dir, stem, "tissue_by_tissue_pr_auc_between_median_river_watershed_lolipop.pdf")
te_tbt_auc_lolipop_plot_with_median_river <- make_tbt_auc_lolipop_plot_with_median_river(te_watershed_obj$roc, te_tissue_names, "Expression", tissues_with_suffiecient_n2_pairs)
# ggsave(te_tbt_auc_lolipop_plot_with_median_river, file=output_file, width=7.2, height=4.0, units="in")










###################################################
# Visualize ASE results
###################################################
outlier_type <- "ase"
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
number_of_dimensions <- 49
######################################
# Visualize theta-pair heatmap to look at correlation structure across tissues
######################################
output_file <- paste0(output_dir, stem, "_tissue_by_tissue_theta_pair_heatmap.pdf")
ase_theta_pair_heatmap <- make_theta_pair_heatmap(ase_watershed_obj$model_params$theta_pair, number_of_dimensions, ase_tissue_names, "ASE")
ggsave(ase_theta_pair_heatmap, file=output_file, width=7.2, height=4.6, units="in")
######################################
# Visualize P(E|Z)
######################################
output_file <- paste0(output_dir, stem, "_tissue_by_tissue_phi.pdf")
watershed_phi_plot <- make_phi_bar_plot(ase_watershed_obj$model_params$phi, ase_tissue_names, paste0("watershed-",outlier_type))
river_phi_plot <- make_phi_bar_plot(ase_river_obj$model_params$phi, ase_tissue_names, paste0("river-",outlier_type))
phi_plot <- plot_grid(river_phi_plot, watershed_phi_plot, ncol=2)
ggsave(phi_plot, file=output_file, width=8.0, height=4.6, units="in")

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(output_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
ase_tbt_auc_lolipop_plot <- make_tbt_auc_lolipop_plot(ase_watershed_obj$roc, ase_river_obj$roc, ase_tissue_names, "ASE", tissues_with_suffiecient_n2_pairs)
# ggsave(ase_tbt_auc_lolipop_plot, file=output_file, width=7.2, height=4.0, units="in")

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(output_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
ase_tbt_auc_lolipop_plot_v2 <- make_tbt_auc_lolipop_plot_v2(ase_watershed_obj$roc, ase_river_obj$roc, ase_tissue_names, "ASE", tissues_with_suffiecient_n2_pairs)
# ggsave(splicing_tbt_auc_lolipop_plot, file=output_file, width=7.2, height=4.0, units="in")


######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(output_dir, stem, "tissue_by_tissue_pr_auc_between_median_river_watershed_lolipop.pdf")
ase_tbt_auc_lolipop_plot_with_median_river <- make_tbt_auc_lolipop_plot_with_median_river(ase_watershed_obj$roc, ase_tissue_names, "ASE", tissues_with_suffiecient_n2_pairs)
# ggsave(ase_tbt_auc_lolipop_plot_with_median_river, file=output_file, width=7.2, height=4.0, units="in")










###################################################
# VISUALIZE CROSS OUTLIER TYPE RESULTS
###################################################

###########################
# Make delta auc distribution plot across tissues
##########################
output_file <- paste0(output_dir, "delta_auc_watershed_river.pdf")
delta_auc_plot <- make_delta_auc_plot(ase_watershed_obj$roc, ase_river_obj$roc, ase_tissue_names, splicing_watershed_obj$roc, splicing_river_obj$roc, splicing_tissue_names, te_watershed_obj$roc, te_river_obj$roc, te_tissue_names, tissues_with_suffiecient_n2_pairs)
ggsave(delta_auc_plot, file=output_file, width=7.2, height=4.0, units="in")

###########################
# Make tbt auc distribution plot across tissues
##########################
output_file <- paste0(output_dir, "tbt_auc_distribution.pdf")
tbt_auc_distribution_plot <- make_tbt_auc_distribution_plot(ase_watershed_obj$roc, ase_river_obj$roc, ase_tissue_names, splicing_watershed_obj$roc, splicing_river_obj$roc, splicing_tissue_names, te_watershed_obj$roc, te_river_obj$roc, te_tissue_names, tissues_with_suffiecient_n2_pairs)
ggsave(tbt_auc_distribution_plot, file=output_file, width=7.2, height=4.0, units="in")

###########################
# Make tbt auc distribution plot across tissues
##########################
output_file <- paste0(output_dir, "tbt_auc_distribution2.pdf")
tbt_auc_distribution_plot2 <- make_tbt_auc_distribution_plot2(ase_watershed_obj$roc, ase_river_obj$roc, ase_tissue_names, splicing_watershed_obj$roc, splicing_river_obj$roc, splicing_tissue_names, te_watershed_obj$roc, te_river_obj$roc, te_tissue_names, tissues_with_suffiecient_n2_pairs)
ggsave(tbt_auc_distribution_plot2, file=output_file, width=7.2, height=4.0, units="in")

######################################
# Visualize tbt pr auc curves between watershed and river across outlier tyeps all on same plot
######################################
output_file <- paste0(output_dir, "tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001_cross_outlier_tissue_by_tissue_pr_auc_between_river_and_watershed_lolipop_by_row.pdf")
auc_tbt_legend <- get_legend(ase_tbt_auc_lolipop_plot)
#combined <- plot_grid(ase_tbt_auc_lolipop_plot + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), splicing_tbt_auc_lolipop_plot + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), te_tbt_auc_lolipop_plot +theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), legend.position="none",plot.margin=grid::unit(c(0,0,0,0), "mm")), rel_heights=c(1,1,2.44), ncol=1)
combined <- plot_grid(ase_tbt_auc_lolipop_plot + ylab("        ")  + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), splicing_tbt_auc_lolipop_plot + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), te_tbt_auc_lolipop_plot+ ylab("        ") +theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), legend.position="none",plot.margin=grid::unit(c(0,0,0,0), "mm")), rel_heights=c(1.24,1.24,3), ncol=1)

tbt_auc_plot <- ggdraw() + draw_plot(combined,0,0,1,1) + draw_plot(auc_tbt_legend,0,-.46,1,1)
ggsave(tbt_auc_plot,file=output_file,width=7.2, height=8.2, units="in")

######################################
# Visualize tbt pr auc curves between watershed and median river across outlier tyeps all on same plot
######################################
output_file <- paste0(output_dir, "tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001_cross_outlier_tissue_by_tissue_pr_auc_between_median_river_and_watershed_lolipop_by_row.pdf")
auc_tbt_legend <- get_legend(ase_tbt_auc_lolipop_plot_with_median_river)
combined <- plot_grid(ase_tbt_auc_lolipop_plot_with_median_river + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), splicing_tbt_auc_lolipop_plot_with_median_river + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), te_tbt_auc_lolipop_plot_with_median_river +theme(legend.position="none",plot.margin=grid::unit(c(0,0,0,0), "mm")), rel_heights=c(1,1,2), ncol=1)
combined2 <- ggdraw() + draw_plot(combined,0,0,1,1) + draw_plot(auc_tbt_legend,.7,.49,1,1)
ggsave(combined2,file=output_file,width=7.2, height=7.2, units="in")

######################################
# Visualize tbt pr auc curves between watershed and median river across outlier tyeps all on same plot
######################################
output_file <- paste0(output_dir, "tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001_cross_outlier_tissue_by_tissue_pr_auc_between_river_and_watershed_lolipop_by_row_v2.pdf")
auc_tbt_legend <- get_legend(ase_tbt_auc_lolipop_plot_v2)
combined <- plot_grid(ase_tbt_auc_lolipop_plot_v2 + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), splicing_tbt_auc_lolipop_plot_v2 + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), te_tbt_auc_lolipop_plot_v2 +theme(legend.position="none",plot.margin=grid::unit(c(0,0,0,0), "mm")), rel_heights=c(1,1,2), ncol=1)
combined2 <- ggdraw() + draw_plot(combined,0,0,1,1) + draw_plot(auc_tbt_legend,.7,.49,1,1)
ggsave(combined2,file=output_file,width=7.2, height=7.2, units="in")



######################################
# Combined heatmap
######################################
output_file <- paste0(output_dir, "combined_tbt_heatmap.pdf")
combined_tbt_heatmap <- plot_grid(ase_theta_pair_heatmap, splicing_theta_pair_heatmap, te_theta_pair_heatmap, ncol=1)
ggsave(combined_tbt_heatmap,file=output_file,width=7.2, height=11, units="in")







############################################
# MAKE FIGURE 4
############################################

row_1 <- plot_grid(NULL, roc_3_theta_pair_exact, roc_3_absolute_risk_bar_plot, ncol=3, labels=c("A","B","C"))
row_2 <- plot_grid(gam_river_watershed_3_class_pr_curves, labels=c("D"))
row_3 <- plot_grid(tbt_auc_plot, labels=c("E"))
fig_4 <- plot_grid(row_1, row_2, row_3, ncol=1, rel_heights=c(.55,.5,1.2))
ggsave(fig_4, file=paste0(output_dir, "fig4.pdf"), width=7.2, height=8, units="in")


############################################
# MAKE FIGURE 4 v2
############################################

row_1 <- plot_grid(NULL, roc_3_theta_pair_exact, roc_3_absolute_risk_bar_plot, ncol=3, labels=c("A","B","C"))
row_2 <- plot_grid(gam_river_watershed_3_class_pr_curves, labels=c("D"))
row_3 <- plot_grid(fig_5_te_theta_pair_heatmap, delta_auc_plot,ncol=2, rel_widths=c(1,.75), labels=c("E", "F"))
fig_4 <- plot_grid(row_1, row_2, row_3, ncol=1, rel_heights=c(.55,.5,.7))
ggsave(fig_4, file=paste0(output_dir, "fig4_v2.pdf"), width=7.2, height=6, units="in")


############################################
# MAKE FIGURE 4 v3
############################################

row_1 <- plot_grid(NULL, roc_3_theta_pair_exact, roc_3_absolute_risk_bar_plot, ncol=3, labels=c("A","B","C"))
row_2 <- plot_grid(gam_river_watershed_3_class_pr_curves, labels=c("D"))
row_3 <- plot_grid(fig_5_te_theta_pair_heatmap, tbt_auc_distribution_plot,ncol=2, rel_widths=c(1,.75), labels=c("E", "F"))
fig_4 <- plot_grid(row_1, row_2, row_3, ncol=1, rel_heights=c(.55,.5,.7))
ggsave(fig_4, file=paste0(output_dir, "fig4_v3.pdf"), width=7.2, height=6, units="in")

############################################
# MAKE FIGURE 4 v4
############################################

row_1 <- plot_grid(NULL, roc_3_theta_pair_exact, roc_3_absolute_risk_bar_plot, ncol=3, labels=c("A","B","C"), rel_widths=c(1,1.1,1))
row_2 <- plot_grid(gam_river_watershed_3_class_pr_curves, labels=c("D"))
row_3 <- plot_grid(fig_5_te_theta_pair_heatmap, tbt_auc_distribution_plot2,ncol=2, rel_widths=c(1,.75), labels=c("E", "F"))
fig_4 <- plot_grid(row_1, row_2, row_3, ncol=1, rel_heights=c(.55,.5,.7))
ggsave(fig_4, file=paste0(output_dir, "fig4_v4.pdf"), width=7.2, height=6, units="in")

