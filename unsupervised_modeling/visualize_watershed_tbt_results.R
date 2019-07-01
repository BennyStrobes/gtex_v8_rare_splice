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
    heatmap <- heatmap + gtex_v8_figure_theme() + theme(axis.text.x=element_blank())
    heatmap <- heatmap + labs(x = "Tissue", y = "Tissue",fill="Theta pair", title=outlier_type)
    return(heatmap)

}



make_tbt_auc_lolipop_plot <- function(roc_object_vi, roc_object_exact, tissue_names, outlier_type) {
	auc_vi <- c()
	auc_exact <- c()
	ordered_tissue_names <- c()
	for (tissue_num in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_num]
		tissue_auc_vi <- roc_object_vi[[tissue_num]]$evaROC$watershed_pr_auc
		tissue_auc_independent <- roc_object_exact[[tissue_num]]$evaROC$watershed_pr_auc
		auc_vi <- c(auc_vi, tissue_auc_vi)
		auc_exact <- c(auc_exact, tissue_auc_independent)
		ordered_tissue_names <- c(ordered_tissue_names, tissue_name)
		if (tissue_name == "Heart_Left_Ventricle" & outlier_type == "ase") {
			auc_vi <- c(auc_vi, NaN)
			auc_exact <- c(auc_exact, NaN)
			ordered_tissue_names <- c(ordered_tissue_names, "Kidney_Cortex")
		}
	}

	cols <- c( "c1" = rgb(0.2,0.7,0.1,0.5), "c2" = rgb(0.7,0.2,0.1,0.5) )

	df <- data.frame(auc_watershed = auc_vi, auc_river=auc_exact, tissue=ordered_tissue_names, tissues_position=1:length(ordered_tissue_names))
	plotter <- ggplot(df) +
  			   geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_vi, yend=auc_exact), color="grey") +
  			   geom_point( aes(x=tissues_position, y=auc_vi, color="c1"), size=1.5) +
               geom_point( aes(x=tissues_position, y=auc_exact, color="c2"), size=1.5) +
               scale_color_manual(name="", breaks=c("c1","c2"), values=cols, labels=c("watershed-tbt", "river-tbt")) +
               ggtitle(outlier_type) +
                xlab("") +
               ylab("AUC (PR)") + 
               theme(legend.position="bottom") +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
               gtex_v8_figure_theme() +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 1.06)) + 
               scale_x_continuous(breaks=1:length(ordered_tissue_names),labels=(ordered_tissue_names)) 
    return(plotter)
}

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

make_tbt_auc_lolipop_plot_with_median_river <- function(roc_object, tissue_names, outlier_type) {
	auc_vi <- c()
	auc_exact <- c()
	ordered_tissue_names <- c()

	for (tissue_num in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_num]
		tissue_auc_vi <- roc_object[[tissue_num]]$evaROC$watershed_pr_auc
		tissue_auc_independent <- roc_object[[tissue_num]]$evaROC$median_river_pr_auc
		auc_vi <- c(auc_vi, tissue_auc_vi)
		auc_exact <- c(auc_exact, tissue_auc_independent)
		ordered_tissue_names <- c(ordered_tissue_names, tissue_name)
		if (tissue_name == "Heart_Left_Ventricle" & outlier_type == "ase") {
			auc_vi <- c(auc_vi, NaN)
			auc_exact <- c(auc_exact, NaN)
			ordered_tissue_names <- c(ordered_tissue_names, "Kidney_Cortex")
		}
	}

	cols <- c( "c1" = rgb(0.2,0.7,0.1,0.5), "c2" = rgb(0.1,0.2,0.6,0.4) )

	df <- data.frame(auc_watershed = auc_vi, auc_river=auc_exact, tissue=ordered_tissue_names, tissues_position=1:length(ordered_tissue_names))
	plotter <- ggplot(df) +
  			   geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_vi, yend=auc_exact), color="grey") +
  			   geom_point( aes(x=tissues_position, y=auc_vi, color="c1"), size=1.5) +
               geom_point( aes(x=tissues_position, y=auc_exact, color="c2"), size=1.5) +
               scale_color_manual(name="", breaks=c("c1","c2"), values=cols, labels=c("watershed-tbt", "river-median")) +
               ggtitle(outlier_type) +
                xlab("") +
               ylab("AUC (PR)") + 
               theme(legend.position="bottom") +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
               gtex_v8_figure_theme() +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 1.06)) + 
               scale_x_continuous(breaks=1:length(ordered_tissue_names),labels=(ordered_tissue_names)) 
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

options(bitmapType = 'cairo', device = 'pdf')
############################
# Command line args
############################
input_dir <- args[1]


############################
# Model hyperparameters
############################
pseudocount <- 30
phi_update_method <- "fixed"
number_of_dimensions <- 49




############################
# Load in splicing data
############################
outlier_type <- "splicing"
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
splicing_watershed_obj <- readRDS(paste0(input_dir, stem, "_inference_pseudolikelihood_independent_false_roc_object.rds"))
splicing_river_obj <- readRDS(paste0(input_dir, stem,"_inference_exact_independent_true_roc_object.rds"))
splicing_tissue_names <- get_tissue_names(splicing_watershed_obj, number_of_dimensions)


############################
# Load in TE data
############################
outlier_type <- "total_expression"
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
te_watershed_obj <- readRDS(paste0(input_dir, stem, "_inference_pseudolikelihood_independent_false_roc_object.rds"))
te_river_obj <- readRDS(paste0(input_dir, stem,"_inference_exact_independent_true_roc_object.rds"))
te_tissue_names <- get_tissue_names(te_watershed_obj, number_of_dimensions)

############################
# Load in ASE data
############################
outlier_type <- "ase"
number_of_dimensions <- 48
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
ase_watershed_obj <- readRDS(paste0(input_dir, stem, "_inference_pseudolikelihood_independent_false_roc_object.rds"))
ase_river_obj <- readRDS(paste0(input_dir, stem,"_inference_exact_independent_true_roc_object.rds"))
ase_tissue_names <- get_tissue_names(ase_watershed_obj, number_of_dimensions)








###################################################
# Visualize splicing results
###################################################
outlier_type <- "splicing"
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
number_of_dimensions <- 49
######################################
# Visualize theta-pair heatmap to look at correlation structure across tissues
######################################
output_file <- paste0(input_dir, stem, "_tissue_by_tissue_theta_pair_heatmap.pdf")
splicing_theta_pair_heatmap <- make_theta_pair_heatmap(splicing_watershed_obj$model_params$theta_pair, number_of_dimensions, splicing_tissue_names, outlier_type)
ggsave(splicing_theta_pair_heatmap, file=output_file, width=7.2, height=4.6, units="in")

######################################
# Visualize P(E|Z)
######################################
output_file <- paste0(input_dir, stem, "_tissue_by_tissue_phi.pdf")
watershed_phi_plot <- make_phi_bar_plot(splicing_watershed_obj$model_params$phi, splicing_tissue_names, paste0("watershed-",outlier_type))
river_phi_plot <- make_phi_bar_plot(splicing_river_obj$model_params$phi, splicing_tissue_names, paste0("river-",outlier_type))
phi_plot <- plot_grid(river_phi_plot, watershed_phi_plot, ncol=2)
ggsave(phi_plot, file=output_file, width=8.0, height=4.6, units="in")


######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
splicing_tbt_auc_lolipop_plot <- make_tbt_auc_lolipop_plot(splicing_watershed_obj$roc, splicing_river_obj$roc, splicing_tissue_names, outlier_type)
# ggsave(splicing_tbt_auc_lolipop_plot, file=output_file, width=7.2, height=4.0, units="in")


######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_median_river_watershed_lolipop.pdf")
splicing_tbt_auc_lolipop_plot_with_median_river <- make_tbt_auc_lolipop_plot_with_median_river(splicing_watershed_obj$roc, splicing_tissue_names, outlier_type)
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
output_file <- paste0(input_dir, stem, "_tissue_by_tissue_theta_pair_heatmap.pdf")
te_theta_pair_heatmap <- make_theta_pair_heatmap(te_watershed_obj$model_params$theta_pair, number_of_dimensions, te_tissue_names, outlier_type)
ggsave(te_theta_pair_heatmap, file=output_file, width=7.2, height=4.6, units="in")
######################################
# Visualize P(E|Z)
######################################
output_file <- paste0(input_dir, stem, "_tissue_by_tissue_phi.pdf")
watershed_phi_plot <- make_phi_bar_plot(te_watershed_obj$model_params$phi, te_tissue_names, paste0("watershed-",outlier_type))
river_phi_plot <- make_phi_bar_plot(te_river_obj$model_params$phi, te_tissue_names, paste0("river-",outlier_type))
phi_plot <- plot_grid(river_phi_plot, watershed_phi_plot, ncol=2)
ggsave(phi_plot, file=output_file, width=8.0, height=4.6, units="in")


######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
te_tbt_auc_lolipop_plot <- make_tbt_auc_lolipop_plot(te_watershed_obj$roc, te_river_obj$roc, te_tissue_names, outlier_type)
# ggsave(te_tbt_auc_lolipop_plot, file=output_file, width=7.2, height=4.0, units="in")


######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_median_river_watershed_lolipop.pdf")
te_tbt_auc_lolipop_plot_with_median_river <- make_tbt_auc_lolipop_plot_with_median_river(te_watershed_obj$roc, te_tissue_names, outlier_type)
# ggsave(te_tbt_auc_lolipop_plot_with_median_river, file=output_file, width=7.2, height=4.0, units="in")











###################################################
# Visualize ASE results
###################################################
outlier_type <- "ase"
stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
number_of_dimensions <- 48
######################################
# Visualize theta-pair heatmap to look at correlation structure across tissues
######################################
output_file <- paste0(input_dir, stem, "_tissue_by_tissue_theta_pair_heatmap.pdf")
ase_theta_pair_heatmap <- make_theta_pair_heatmap(ase_watershed_obj$model_params$theta_pair, number_of_dimensions, ase_tissue_names, outlier_type)
ggsave(ase_theta_pair_heatmap, file=output_file, width=7.2, height=4.6, units="in")
######################################
# Visualize P(E|Z)
######################################
output_file <- paste0(input_dir, stem, "_tissue_by_tissue_phi.pdf")
watershed_phi_plot <- make_phi_bar_plot(ase_watershed_obj$model_params$phi, ase_tissue_names, paste0("watershed-",outlier_type))
river_phi_plot <- make_phi_bar_plot(ase_river_obj$model_params$phi, ase_tissue_names, paste0("river-",outlier_type))
phi_plot <- plot_grid(river_phi_plot, watershed_phi_plot, ncol=2)
ggsave(phi_plot, file=output_file, width=8.0, height=4.6, units="in")

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
ase_tbt_auc_lolipop_plot <- make_tbt_auc_lolipop_plot(ase_watershed_obj$roc, ase_river_obj$roc, ase_tissue_names, outlier_type)
# ggsave(ase_tbt_auc_lolipop_plot, file=output_file, width=7.2, height=4.0, units="in")


######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_median_river_watershed_lolipop.pdf")
ase_tbt_auc_lolipop_plot_with_median_river <- make_tbt_auc_lolipop_plot_with_median_river(ase_watershed_obj$roc, ase_tissue_names, outlier_type)
# ggsave(ase_tbt_auc_lolipop_plot_with_median_river, file=output_file, width=7.2, height=4.0, units="in")














###################################################
# VISUALIZE CROSS OUTLIER TYPE RESULTS
###################################################
######################################
# Visualize tbt pr auc curves between watershed and river across outlier tyeps all on same plot
######################################
output_file <- paste0(input_dir, "tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001_cross_outlier_tissue_by_tissue_pr_auc_between_river_and_watershed_lolipop_on_same_axis.pdf")
make_tbt_auc_lolipop_plot_across_outliers_on_same_axis(splicing_watershed_obj$roc, splicing_river_obj$roc, te_watershed_obj$roc, te_river_obj$roc, ase_watershed_obj$roc, ase_river_obj$roc, splicing_tissue_names, te_tissue_names, ase_tissue_names, output_file)

######################################
# Visualize tbt pr auc curves between watershed and river across outlier tyeps all on same plot
######################################
output_file <- paste0(input_dir, "tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001_cross_outlier_tissue_by_tissue_pr_auc_between_river_and_watershed_lolipop_by_row.pdf")
auc_tbt_legend <- get_legend(ase_tbt_auc_lolipop_plot)
combined <- plot_grid(ase_tbt_auc_lolipop_plot + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), splicing_tbt_auc_lolipop_plot + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), te_tbt_auc_lolipop_plot +theme(legend.position="none",plot.margin=grid::unit(c(0,0,0,0), "mm")), rel_heights=c(1,1,2), ncol=1)
combined2 <- ggdraw() + draw_plot(combined,0,0,1,1) + draw_plot(auc_tbt_legend,.75,.49,1,1)
ggsave(combined2,file=output_file,width=7.2, height=7.2, units="in")

######################################
# Visualize tbt pr auc curves between watershed and median river across outlier tyeps all on same plot
######################################
output_file <- paste0(input_dir, "tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001_cross_outlier_tissue_by_tissue_pr_auc_between_median_river_and_watershed_lolipop_by_row.pdf")
auc_tbt_legend <- get_legend(ase_tbt_auc_lolipop_plot_with_median_river)
combined <- plot_grid(ase_tbt_auc_lolipop_plot_with_median_river + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), splicing_tbt_auc_lolipop_plot_with_median_river + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), te_tbt_auc_lolipop_plot_with_median_river +theme(legend.position="none",plot.margin=grid::unit(c(0,0,0,0), "mm")), rel_heights=c(1,1,2), ncol=1)
combined2 <- ggdraw() + draw_plot(combined,0,0,1,1) + draw_plot(auc_tbt_legend,.7,.49,1,1)
ggsave(combined2,file=output_file,width=7.2, height=7.2, units="in")
