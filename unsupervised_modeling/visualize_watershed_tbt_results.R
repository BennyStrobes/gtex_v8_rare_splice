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



make_theta_pair_heatmap <- function(theta_pair, number_of_dimensions,tissue_names, output_file, outlier_type) {
	options(bitmapType = 'cairo', device = 'pdf')

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
	corr_mat = matrix(0, number_of_dimensions, number_of_dimensions)
	dimension_counter = 1
	for (dimension1 in 1:number_of_dimensions) {
		for (dimension2 in dimension1:number_of_dimensions) {
			if (dimension1 != dimension2) {
				corr_mat[dimension1, dimension2] = cor(theta_pair_mat[dimension1,], theta_pair_mat[dimension2,], use="pairwise.complete.obs", method="spearman")
				corr_mat[dimension2, dimension1] = cor(theta_pair_mat[dimension1,], theta_pair_mat[dimension2,], use="pairwise.complete.obs", method="spearman")
				dimension_counter = dimension_counter + 1
			}
		}
		corr_mat[dimension1, dimension1] = NA
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
    heatmap <- heatmap + theme(legend.text = element_text(size=7),legend.title = element_text(size=7),text = element_text(size=7),plot.title = element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    heatmap <- heatmap + labs(x = "Tissue", y = "Tissue",fill="Theta pair", title=outlier_type)

    ggsave(heatmap, file=output_file,width=7.2, height=6.0, units="in")
}



make_tbt_auc_lolipop_plot <- function(roc_object_vi, roc_object_exact, tissue_names, output_file, outlier_type) {
	options(bitmapType = 'cairo', device = 'pdf')
	auc_vi <- c()
	auc_exact <- c()
	for (tissue_num in length(tissue_names):1) {
		tissue_name <- tissue_names[tissue_num]
		tissue_auc_vi <- roc_object_vi[[tissue_num]]$evaROC$watershed_pr_auc
		tissue_auc_independent <- roc_object_exact[[tissue_num]]$evaROC$watershed_pr_auc
		auc_vi <- c(auc_vi, tissue_auc_vi)
		auc_exact <- c(auc_exact, tissue_auc_independent)
	}

	cols <- c( "c1" = rgb(0.2,0.7,0.1,0.5), "c2" = rgb(0.7,0.2,0.1,0.5) )

	df <- data.frame(auc_watershed = auc_vi, auc_river=auc_exact, tissue=tissue_names, tissues_position=1:length(tissue_names))
	plotter <- ggplot(df) +
  			   geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_vi, yend=auc_exact), color="grey") +
  			   geom_point( aes(x=tissues_position, y=auc_vi, color="c1"), size=3) +
               geom_point( aes(x=tissues_position, y=auc_exact, color="c2"), size=3) +
               scale_color_manual(name="", breaks=c("c1","c2"), values=cols, labels=c("watershed-tbt", "river-tbt")) +
               coord_flip()+
               ggtitle(outlier_type) +
                xlab("") +
               ylab("AUC (PR)") + 
               theme(legend.position="bottom") +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
               scale_x_continuous(breaks=1:length(tissue_names),labels=rev(tissue_names)) 
    ggsave(plotter, file=output_file, width=10, height=7.0, units="in")
}

make_tbt_auc_lolipop_plot_with_median_river <- function(roc_object, tissue_names, output_file, outlier_type) {
	options(bitmapType = 'cairo', device = 'pdf')
	auc_vi <- c()
	auc_exact <- c()
	for (tissue_num in length(tissue_names):1) {
		tissue_name <- tissue_names[tissue_num]
		tissue_auc_vi <- roc_object[[tissue_num]]$evaROC$watershed_pr_auc
		tissue_auc_independent <- roc_object[[tissue_num]]$evaROC$median_river_pr_auc
		auc_vi <- c(auc_vi, tissue_auc_vi)
		auc_exact <- c(auc_exact, tissue_auc_independent)
	}

	cols <- c( "c1" = rgb(0.2,0.7,0.1,0.5), "c2" = rgb(0.1,0.2,0.6,0.4) )

	df <- data.frame(auc_watershed = auc_vi, auc_river=auc_exact, tissue=tissue_names, tissues_position=1:length(tissue_names))
	plotter <- ggplot(df) +
  			   geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_vi, yend=auc_exact), color="grey") +
  			   geom_point( aes(x=tissues_position, y=auc_vi, color="c1"), size=3) +
               geom_point( aes(x=tissues_position, y=auc_exact, color="c2"), size=3) +
               scale_color_manual(name="", breaks=c("c1","c2"), values=cols, labels=c("watershed-tbt", "river-median")) +
               coord_flip()+
                xlab("") +
                ggtitle(outlier_type) +
               ylab("AUC (PR)") + 
               theme(legend.position="bottom") +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
               scale_x_continuous(breaks=1:length(tissue_names),labels=rev(tissue_names))
    ggsave(plotter, file=output_file, width=10, height=7.0, units="in")
}




input_dir <- args[1]






outlier_type <- "splicing"
pseudocount <- 30
phi_update_method <- "sample_size"
number_of_dimensions <- 49

stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
watershed_obj <- readRDS(paste0(input_dir, stem, "_inference_pseudolikelihood_independent_false_roc_object.rds"))
river_obj <- readRDS(paste0(input_dir, stem,"_inference_exact_independent_true_roc_object.rds"))

tissue_names <- get_tissue_names(watershed_obj, number_of_dimensions)

######################################
# Visualize theta-pair heatmap to look at correlation structure across tissues
######################################
output_file <- paste0(input_dir, stem, "_tissue_by_tissue_theta_pair_heatmap.pdf")
make_theta_pair_heatmap(watershed_obj$model_params$theta_pair, number_of_dimensions, tissue_names, output_file, outlier_type)

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
make_tbt_auc_lolipop_plot(watershed_obj$roc, river_obj$roc, tissue_names, output_file, outlier_type)

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_median_river_watershed_lolipop.pdf")
make_tbt_auc_lolipop_plot_with_median_river(watershed_obj$roc, tissue_names, output_file, outlier_type)



outlier_type <- "total_expression"
pseudocount <- 30
phi_update_method <- "sample_size"
number_of_dimensions <- 49

stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
watershed_obj <- readRDS(paste0(input_dir, stem, "_inference_pseudolikelihood_independent_false_roc_object.rds"))
river_obj <- readRDS(paste0(input_dir, stem,"_inference_exact_independent_true_roc_object.rds"))

tissue_names <- get_tissue_names(watershed_obj, number_of_dimensions)

######################################
# Visualize theta-pair heatmap to look at correlation structure across tissues
######################################
output_file <- paste0(input_dir, stem, "_tissue_by_tissue_theta_pair_heatmap.pdf")
make_theta_pair_heatmap(watershed_obj$model_params$theta_pair, number_of_dimensions, tissue_names, output_file, outlier_type)

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
make_tbt_auc_lolipop_plot(watershed_obj$roc, river_obj$roc, tissue_names, output_file, outlier_type)

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_median_river_watershed_lolipop.pdf")
make_tbt_auc_lolipop_plot_with_median_river(watershed_obj$roc, tissue_names, output_file, outlier_type)


outlier_type <- "ase"
pseudocount <- 30
phi_update_method <- "sample_size"
number_of_dimensions <- 48

stem <- paste0(outlier_type,"_tbt_intersect_te_ase_splicing_out_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_", pseudocount, "_", phi_update_method, "_.001_.001")
watershed_obj <- readRDS(paste0(input_dir, stem, "_inference_pseudolikelihood_independent_false_roc_object.rds"))
river_obj <- readRDS(paste0(input_dir, stem,"_inference_exact_independent_true_roc_object.rds"))

tissue_names <- get_tissue_names(watershed_obj, number_of_dimensions)

######################################
# Visualize theta-pair heatmap to look at correlation structure across tissues
######################################
output_file <- paste0(input_dir, stem, "_tissue_by_tissue_theta_pair_heatmap.pdf")
make_theta_pair_heatmap(watershed_obj$model_params$theta_pair, number_of_dimensions, tissue_names, output_file, outlier_type)

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
make_tbt_auc_lolipop_plot(watershed_obj$roc, river_obj$roc, tissue_names, output_file, outlier_type)

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
output_file <- paste0(input_dir, stem, "tissue_by_tissue_pr_auc_between_median_river_watershed_lolipop.pdf")
make_tbt_auc_lolipop_plot_with_median_river(watershed_obj$roc, tissue_names, output_file, outlier_type)



