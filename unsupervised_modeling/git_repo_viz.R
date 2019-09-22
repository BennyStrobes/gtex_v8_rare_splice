options(bitmapType = 'cairo', device = 'pdf')
args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(PRROC)
library(cowplot)
library(RColorBrewer)





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
gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

###################

plot_pr_gam_river_watershed_comparison_curve <- function(roc_object_exact, roc_object_ind, number_of_dimensions, output_file) {
	precision <- c()
	recall <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object_exact[[dimension]]
		dimension_name <- dimension_roc_object$name
		print(dimension_name)
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
  	print(summary(df))

	outlier_type <- "pheno_1"
	outlier_name <- "ASpheno_1E"
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

	outlier_type <- "pheno_2"
	outlier_name <- "pheno_2"
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

	outlier_type <- "pheno_2"
	outlier_name <- "pheno_2"
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
  

	outlier_type <- "pheno_1"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type, linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title="ASE") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                gtex_v8_figure_theme()
                #theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	outlier_type <- "pheno_2"
  	plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type, linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title="Splicing") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                gtex_v8_figure_theme()
                #theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	outlier_type <- "pheno_3"
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
	outlier_class <- c(rep("pheno_1", length(model_params_approximate$theta[,1])), rep("pheno_2", length(model_params_approximate$theta[,2])), rep("pheno_3", length(model_params_approximate$theta[,3])))

	df <- data.frame(exact_betas=exact_betas, approximate_betas=approximate_betas, outlier_class=factor(outlier_class, levels=c("pheno_1", "pheno_2", "pheno_3")))

	plotter <- ggplot(df, aes(x=exact_betas, y=approximate_betas, colour=outlier_class)) + geom_point() +
			geom_abline() + 
			theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)) +
			labs(x="Beta (exact)", y="Beta (approximate)",colour="") + 
			scale_color_manual(values=c("#7F5A83", "#0D324D", "#BFCDE0")) +
			gtex_v8_figure_theme()
	return(plotter)
}



input_stem = args[1]


####### input data
roc_3_class_data_input <- readRDS(paste0(input_stem, "_data_input.rds"))

####### Exact watershed
independent_variables = "false"
inference_method = "exact"
output_root <- paste0(input_stem,"_inference_", inference_method, "_independent_", independent_variables)
roc_object_exact <- readRDS(paste0(output_root, "_roc_object.rds"))


####### Pseudolikelihood approximation to watershed
independent_variables = "false"
inference_method = "pseudolikelihood"
output_root <- paste0(input_stem,"_inference_", inference_method, "_independent_", independent_variables)
roc_object_pseudo <- readRDS(paste0(output_root, "_roc_object.rds"))


####### Exact RIVER
independent_variables = "true"
inference_method = "exact"
output_root <- paste0(input_stem,"_inference_", inference_method, "_independent_", independent_variables)
roc_object_independent <- readRDS(paste0(output_root, "_roc_object.rds"))


#######################################
## Visualize theta pair terms for exact inference
#######################################
output_file <- paste0(input_stem, "watershed_3_class_roc_exact_inference_theta_pair_heatmap.pdf")
roc_3_theta_pair_exact <- plot_three_class_theta_pair_term(roc_object_exact$model_params$theta_pair)
ggsave(roc_3_theta_pair_exact, file=output_file, width=7.2, height=4, units="in")

#######################################
## Visualize precision-recall curves for river, GAM, watershed-exact comparison and all three outlier types (te, splice, ase)
#######################################
number_of_dimensions <- 3
output_file <- paste0(input_stem, "watershed_3_class_roc_exact_gam_river_watershed_precision_recall.pdf")
gam_river_watershed_3_class_pr_curves <- plot_pr_gam_river_watershed_comparison_curve(roc_object_exact$roc, roc_object_independent$roc, number_of_dimensions)
ggsave(gam_river_watershed_3_class_pr_curves, file=output_file, width=10, height=3.0, units="in")




#######################################
## Visualize precision-recall curves for river, watershed-vi, watershed-exact comparison and all three outlier types (te, splice, ase)
#######################################
gam_river_watershed_3_class_pr_curves_vi_cmp <- plot_pr_river_watershed_comparison_curve(roc_object_exact$roc, roc_object_pseudo$roc, roc_object_independent$roc, number_of_dimensions)
#ggsave(gam_river_watershed_3_class_pr_curves_vi_cmp, file=output_file, width=7.2, height=3.0, units="in")


#######################################
## Visualize beta differences for exact vs inference
#######################################
beta_comparison_scatter <- plot_beta_difference_scatter_between_exact_and_vi(roc_object_exact$model_params, roc_object_pseudo$model_params)
#ggsave(beta_comparison_scatter, file=output_file, width=7.2, height=3.0, units="in")


#######################################
## Combined supplementary VI plot for 3class model
#######################################
output_file <- paste0(input_stem, "watershed_3_class_exact_vi_comparison_supplementary_figure.pdf")
combined <- plot_grid(beta_comparison_scatter, gam_river_watershed_3_class_pr_curves_vi_cmp, ncol=1,labels=c("A","B"))
ggsave(combined, file=output_file, width=7.2, height=5.0, units="in")


