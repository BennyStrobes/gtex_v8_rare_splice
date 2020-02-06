args = commandArgs(trailingOnly=TRUE)
library(glmnet)
library(methods)
library(stats)
library(utils)
library(Biobase)
library(pROC)
library(ggplot2)
library(sigmoid)
library(Rcpp)
library(numDeriv)
library(lbfgs)
library(reshape)
library(grid)
library(PRROC)
library(cowplot)
library(RColorBrewer)
source("watershed.R")







# x is 1-specificity (false positive rate)
# y is sensitivity  (true positive rate)
plot_pr_all_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	precision <- c()
	recall <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		precision <- c(precision, dimension_roc_object$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_precision)))
		# independent
		precision <- c(precision, dimension_roc_object_ind$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object_ind$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("river", length(dimension_roc_object_ind$evaROC$watershed_precision)))
		# rna only
		precision <- c(precision, dimension_roc_object$evaROC$rna_only_precision)
		recall <- c(recall, dimension_roc_object$evaROC$rna_only_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$rna_only_precision)))
		prediction_type <- c(prediction_type, rep("rna only", length(dimension_roc_object$evaROC$rna_only_precision)))
		# GAM
		precision <- c(precision, dimension_roc_object$evaROC$GAM_precision)
		recall <- c(recall, dimension_roc_object$evaROC$GAM_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$GAM_precision)))
		prediction_type <- c(prediction_type, rep("GAM", length(dimension_roc_object$evaROC$GAM_precision)))
	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type))
  

  	plotter <- ggplot(data=df, aes(x=recall, y=precision, colour=prediction_type)) + geom_line() + facet_wrap( ~ outlier_type, ncol=3) +
                labs(x="Recall", y="Precision", colour="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10))

	ggsave(plotter, file=output_file,width = 26,height=11,units="cm")
}



# x is 1-specificity (false positive rate)
# y is sensitivity  (true positive rate)
plot_roc_all_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	tpr <- c()
	fpr <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		tpr <- c(tpr, dimension_roc_object$evaROC$watershed_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$watershed_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_sens)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_sens)))
		# independent
		tpr <- c(tpr, dimension_roc_object_ind$evaROC$watershed_sens)
		fpr <- c(fpr, 1-dimension_roc_object_ind$evaROC$watershed_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_sens)))
		prediction_type <- c(prediction_type, rep("river", length(dimension_roc_object_ind$evaROC$watershed_sens)))
		# rna only
		tpr <- c(tpr, dimension_roc_object$evaROC$rna_only_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$rna_only_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$rna_only_sens)))
		prediction_type <- c(prediction_type, rep("rna only", length(dimension_roc_object$evaROC$rna_only_sens)))
		# GAM
		tpr <- c(tpr, dimension_roc_object$evaROC$GAM_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$GAM_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$GAM_sens)))
		prediction_type <- c(prediction_type, rep("GAM", length(dimension_roc_object$evaROC$GAM_sens)))
	}
	df <- data.frame(tpr=tpr, fpr=fpr, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type))
  

  	plotter <- ggplot(data=df, aes(x=fpr, y=tpr, colour=prediction_type)) + geom_line() + facet_wrap( ~ outlier_type, ncol=3) +
  				geom_abline(slope=1) +
                labs(x="False positive rate", y="True positive rate", colour="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10))

	ggsave(plotter, file=output_file,width = 26,height=11,units="cm")
}



# x is 1-specificity (false positive rate)
# y is sensitivity  (true positive rate)
plot_roc_gam_watershed_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	tpr <- c()
	fpr <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		tpr <- c(tpr, dimension_roc_object$evaROC$watershed_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$watershed_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_sens)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_sens)))
		# GAM
		tpr <- c(tpr, dimension_roc_object$evaROC$GAM_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$GAM_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$GAM_sens)))
		prediction_type <- c(prediction_type, rep("GAM", length(dimension_roc_object$evaROC$GAM_sens)))
	}
	df <- data.frame(tpr=tpr, fpr=fpr, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("watershed","GAM")))
  

  	plotter <- ggplot(data=df, aes(x=fpr, y=tpr, colour=outlier_type, linetype=prediction_type)) + geom_line() +
  				geom_abline(slope=1) +
                labs(x="False positive rate", y="True positive rate", colour="",linetype="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="right") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	ggsave(plotter, file=output_file,width = 15,height=11,units="cm")
}

# x is 1-specificity (false positive rate)
# y is sensitivity  (true positive rate)
plot_roc_river_watershed_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	tpr <- c()
	fpr <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		tpr <- c(tpr, dimension_roc_object$evaROC$watershed_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$watershed_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_sens)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_sens)))
		# GAM
		tpr <- c(tpr, dimension_roc_object_ind$evaROC$watershed_sens)
		fpr <- c(fpr, 1-dimension_roc_object_ind$evaROC$watershed_sens)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_sens)))
		prediction_type <- c(prediction_type, rep("river", length(dimension_roc_object_ind$evaROC$watershed_sens)))
	}
	df <- data.frame(tpr=tpr, fpr=fpr, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("watershed","river")))
  

  	plotter <- ggplot(data=df, aes(x=fpr, y=tpr, colour=outlier_type, linetype=prediction_type)) + geom_line() +
  				geom_abline(slope=1) +
                labs(x="False positive rate", y="True positive rate", colour="",linetype="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="right") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	ggsave(plotter, file=output_file,width = 15,height=11,units="cm")
}

plot_pr_gam_watershed_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	precision <- c()
	recall <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		precision <- c(precision, dimension_roc_object$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_precision)))

		# GAM
		precision <- c(precision, dimension_roc_object$evaROC$GAM_precision)
		recall <- c(recall, dimension_roc_object$evaROC$GAM_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$GAM_precision)))
		prediction_type <- c(prediction_type, rep("GAM", length(dimension_roc_object$evaROC$GAM_precision)))
	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("watershed","GAM")))
  

  	plotter <- ggplot(data=df, aes(x=recall, y=precision, colour=outlier_type, linetype=prediction_type)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="",linetype="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="right") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	ggsave(plotter, file=output_file,width = 19,height=11,units="cm")
}

plot_pr_river_watershed_comparison_curve <- function(roc_object_exact, roc_object_vi, roc_object_ind, number_of_dimensions, output_file) {
	options(bitmapType = 'cairo', device = 'pdf')
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
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_precision)))
		inference_method <- c(inference_method, rep("exact", length(dimension_roc_object$evaROC$watershed_precision)))


		# Tied watershed (variational inference)
		precision <- c(precision, dimension_roc_object_vi$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object_vi$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_vi$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object_vi$evaROC$watershed_precision)))
		inference_method <- c(inference_method, rep("vi", length(dimension_roc_object_vi$evaROC$watershed_precision)))


		# Indepdent
		precision <- c(precision, dimension_roc_object_ind$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object_ind$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("river", length(dimension_roc_object_ind$evaROC$watershed_precision)))
		inference_method <- c(inference_method, rep("exact", length(dimension_roc_object_ind$evaROC$watershed_precision)))

	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), inference=factor(inference_method, levels=c("exact","vi")), prediction_type=factor(prediction_type, levels=c("watershed","river")))
  

	outlier_type <- "ase"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type, linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title=outlier_type) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	outlier_type <- "splicing"
  	plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type, linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title=outlier_type) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	outlier_type <- "total_expression"
  	plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type,linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title=outlier_type) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

    legend <- get_legend(plotter_ase)
    combined_plots <- plot_grid(plotter_ase + theme(legend.position="none"), plotter_splice+ theme(legend.position="none"), plotter_te+ theme(legend.position="none"), nrow=1)


	ggsave(plot_grid(combined_plots, legend,ncol=1, rel_heights=c(1,.1)), file=output_file,width = 34,height=11,units="cm")
}


plot_pr_gam_river_watershed_comparison_curve <- function(roc_object_exact, roc_object_ind, number_of_dimensions, output_file) {
	options(bitmapType = 'cairo', device = 'pdf')
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
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_precision)))


		# GAM
		precision <- c(precision, dimension_roc_object$evaROC$GAM_precision)
		recall <- c(recall, dimension_roc_object$evaROC$GAM_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$GAM_precision)))
		prediction_type <- c(prediction_type, rep("GAM", length(dimension_roc_object$evaROC$GAM_precision)))


		# Indepdent
		precision <- c(precision, dimension_roc_object_ind$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object_ind$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("river", length(dimension_roc_object_ind$evaROC$watershed_precision)))

	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("watershed","river", "GAM")))
  

	outlier_type <- "ase"
	outlier_name <- "ASE"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", title=outlier_name) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	outlier_type <- "splicing"
	outlier_name <- "Splice"
  	plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", title=outlier_name) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	outlier_type <- "total_expression"
	outlier_name <- "Total expression"
  	plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", title=outlier_name) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

    legend <- get_legend(plotter_ase + theme(legend.position="right"))
    combined_plots <- plot_grid(plotter_te+ theme(legend.position="none"),plotter_ase + theme(legend.position="none"), plotter_splice+ theme(legend.position="none"), legend, rel_widths=c(1,1,1,.26), nrow=1)


	ggsave(combined_plots, file=output_file,width = 40,height=9,units="cm")
}

plot_pr_watershed_watershed_edge_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	options(bitmapType = 'cairo', device = 'pdf')

	precision <- c()
	recall <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		precision <- c(precision, dimension_roc_object$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_precision)))

		# Indepdent
		precision <- c(precision, dimension_roc_object_ind$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object_ind$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("watershed_edge", length(dimension_roc_object_ind$evaROC$watershed_precision)))
	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("watershed","watershed_edge")))
  

  	plotter <- ggplot(data=df, aes(x=recall, y=precision, colour=outlier_type, linetype=prediction_type)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="",linetype="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="right") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	ggsave(plotter, file=output_file,width = 19,height=11,units="cm")
}

# Helper function to remove NAs from vector
remove_na <- function(x) {
	return(x[!is.na(x)])
}

compute_confidence_interval_on_auprc_with_non_parametric_bootstrap <- function(predictions, labels, num_bootstrap_samples) {
	set.seed(2)
	# Number of samples we have
	num_samples <- length(predictions)
	# Initialize vector to keep track of auprc across bootstrap iterations
	auprc_arr <- c()
	# Compute aucprc across whole vector
	population_auprc <- pr.curve(scores.class0=predictions[labels==1], scores.class1=predictions[labels==0])$auc.integral
	# Loop through bootstrap samples
	for (bootstrap_sample in 1:num_bootstrap_samples) {
		bootstrap_indices = sample(1:num_samples,size=num_samples,replace=TRUE)
		bootstrap_predictions = predictions[bootstrap_indices]
		bootstrap_labels = labels[bootstrap_indices]
		bootstrap_auprc <- pr.curve(scores.class0=bootstrap_predictions[bootstrap_labels==1], scores.class1=bootstrap_predictions[bootstrap_labels==0])$auc.integral
		auprc_arr <- c(auprc_arr, bootstrap_auprc)
	}
	#print(quantile(auprc_arr,.025))
	#print(quantile(auprc_arr, .975))
	return(auprc_arr)
}


#######################################
# Extract ROC curves and precision recall curves for test set (in each dimension seperately) using:
#### 1. Watershed
#### 2. GAM
#### 3. RNA-only
#######################################
compute_roc_across_dimensions <- function(number_of_dimensions, dimension_labels, posterior_prob_test, real_valued_outliers_test1, gam_posteriors, binary_outliers_test2, median_river_prob, outlier_type) {
	num_non_parametric_bootstrap_samples <- 20000
	roc_object_across_dimensions <- list()
	pos_list <- c()
	neg_list <- c()
	if (outlier_type == "splicing") {
		outlier_dimension <- 1
	} else if (outlier_type == "total_expression") {
		outlier_dimension <- 2
	} else if (outlier_type == "ase") {
		outlier_dimension <- 3
	}
  	# Loop through dimensions
  	for (dimension in 1:number_of_dimensions) {
  		print(dimension)
  		# Name of dimension
  		dimension_name <- strsplit(dimension_labels[dimension],"_pval")[[1]][1]
  		# Pseudo gold standard
  		test_outlier_status <- binary_outliers_test2[,dimension]
  		# river predictions
  		# roc_obj <- roc(test_outlier_status, posterior_prob_test[,dimension])
  		roc_obj <- roc.curve(scores.class0 = remove_na(posterior_prob_test[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(posterior_prob_test[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
  		pr_obj <- pr.curve(scores.class0 = remove_na(posterior_prob_test[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(posterior_prob_test[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
  		
  		pos_list <- c(pos_list, remove_na(posterior_prob_test[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]))
  		neg_list <- c(neg_list, remove_na(posterior_prob_test[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]))
  	
  		valid_indices <- !is.na(real_valued_outliers_test1[,dimension]) & !is.na(test_outlier_status)
  		watershed_auc_prc_ci <- compute_confidence_interval_on_auprc_with_non_parametric_bootstrap(posterior_prob_test[valid_indices, dimension], test_outlier_status[valid_indices], num_non_parametric_bootstrap_samples)

  		# Predictions with only RNA
  		#rna_only_roc_obj <- roc(test_outlier_status, real_valued_outliers_test1[,dimension])
  		rna_only_roc_obj <- roc.curve(scores.class0 = remove_na(real_valued_outliers_test1[,dimension][test_outlier_status==1]), scores.class1 = remove_na(real_valued_outliers_test1[,dimension][test_outlier_status==0]), curve = T)
  		rna_only_pr_obj <- pr.curve(scores.class0 = remove_na(real_valued_outliers_test1[,dimension][test_outlier_status==1]), scores.class1 = remove_na(real_valued_outliers_test1[,dimension][test_outlier_status==0]), curve = T)

  		# predictions with only genomic annotations
  		#gam_roc_obj <- roc(test_outlier_status, gam_posteriors[,dimension])
   		gam_roc_obj <- roc.curve(scores.class0 = remove_na(gam_posteriors[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(gam_posteriors[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
   		gam_pr_obj <- pr.curve(scores.class0 = remove_na(gam_posteriors[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(gam_posteriors[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)

   		valid_indices <- !is.na(real_valued_outliers_test1[,dimension]) & !is.na(test_outlier_status)
  		gam_auc_prc_ci <- compute_confidence_interval_on_auprc_with_non_parametric_bootstrap(gam_posteriors[valid_indices, dimension], test_outlier_status[valid_indices], num_non_parametric_bootstrap_samples)

   		median_river_roc_obj <- roc.curve(scores.class0 = remove_na(median_river_prob[,outlier_dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(median_river_prob[,outlier_dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
   		median_river_pr_obj <- pr.curve(scores.class0 = remove_na(median_river_prob[,outlier_dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(median_river_prob[,outlier_dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)

   		valid_indices <- !is.na(real_valued_outliers_test1[,dimension]) & !is.na(test_outlier_status)
  		median_river_auc_prc_ci <- compute_confidence_interval_on_auprc_with_non_parametric_bootstrap(median_river_prob[valid_indices, outlier_dimension], test_outlier_status[valid_indices], num_non_parametric_bootstrap_samples)


		evaROC <-	
		 list(watershed_sens=roc_obj$curve[,2],
              watershed_spec=1-roc_obj$curve[,1],
         	  watershed_auc=roc_obj$auc,
         	  watershed_pr_auc=pr_obj$auc.integral,
         	  watershed_recall=pr_obj$curve[,1],
         	  watershed_precision=pr_obj$curve[,2],
         	  watershed_pr_auc_bootstraps=watershed_auc_prc_ci,
         	  GAM_sens=gam_roc_obj$curve[,2],
              GAM_spec=1-gam_roc_obj$curve[,1],
              GAM_auc=gam_roc_obj$auc,
         	  GAM_pr_auc=gam_pr_obj$auc.integral,
         	  GAM_recall=gam_pr_obj$curve[,1],
         	  GAM_precision=gam_pr_obj$curve[,2],
         	  GAM_pr_auc_bootstraps=gam_auc_prc_ci,
         	  rna_only_pr_auc=rna_only_pr_obj$auc.integral,
         	  rna_only_recall=rna_only_pr_obj$curve[,1],
         	  rna_only_precision=rna_only_pr_obj$curve[,2],
              rna_only_sens=rna_only_roc_obj$curve[,2],
              rna_only_spec=1-rna_only_roc_obj$curve[,1],
              rna_only_auc=rna_only_roc_obj$auc,
              median_river_sens=median_river_roc_obj$curve[,2],
              median_river_spec=1-median_river_roc_obj$curve[,1],
              median_river_auc=median_river_roc_obj$auc,
              median_river_pr_auc=median_river_pr_obj$auc.integral,
              median_river_recall=median_river_pr_obj$curve[,1],
              median_river_precision=median_river_pr_obj$curve[,2],
              median_river_pr_auc_bootstraps=median_river_auc_prc_ci,
              num_positive_pairs=length(remove_na(posterior_prob_test[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])])),
              num_negative_pairs=length(remove_na(posterior_prob_test[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])])))


		 roc_object_across_dimensions[[dimension]] <- list(name=dimension_name, evaROC=evaROC)
	}

		pr_obj <- pr.curve(scores.class0=pos_list, scores.class1=neg_list, curve = T)

		evaROC <-	
		 list(watershed_pr_auc=pr_obj$auc.integral,
         	  watershed_recall=pr_obj$curve[,1],
         	  watershed_precision=pr_obj$curve[,2])

		roc_object_across_dimensions[[(number_of_dimensions + 1)]] <- list(name="joint", evaROC=evaROC)

	return(roc_object_across_dimensions)
}



roc_analysis <- function(data_input, number_of_dimensions, lambda_costs, pseudoc, inference_method, independent_variables, vi_step_size, vi_threshold, phi_method, lambda_init, lambda_pair_init, outlier_type, file_name) {
	watershed_data <- readRDS(file_name)
	watershed_model <- watershed_data$model_params
	gam_data <- watershed_data$gam_params
	#######################################
	# Load in all data (training and test)
	#######################################
	feat_all <- data_input$feat
	discrete_outliers_all <- data_input$outliers_discrete
	binary_outliers_all <- data_input$outliers_binary
	fraction_binary_outliers_all <- data_input$fraction_outliers_binary
	N2_pairs <- data_input$N2_pairs

	#######################################
	# Extract training data
	#######################################
	feat_train <- feat_all[is.na(N2_pairs),]
	discrete_outliers_train <- discrete_outliers_all[is.na(N2_pairs),]
	binary_outliers_train <-  binary_outliers_all[is.na(N2_pairs),]

	#######################################
	# Extract Test data
	#######################################
  	feat_test <- rbind(feat_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], feat_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	discrete_outliers_test1 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	discrete_outliers_test2 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])
  	binary_outliers_test1 <- rbind(binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	binary_outliers_test2 <- rbind(fraction_binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], fraction_binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])
  	# Absolute pvalues from test prediction data set (to be used for RNA-only analysis)
  	real_valued_outliers_test1 <- -log10(abs(rbind(data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])) + 1e-7)

	#######################################
	## Standardize Genomic Annotations (features)
	#######################################
	mean_feat <- apply(feat_all, 2, mean)
	sd_feat <- apply(feat_all, 2, sd)
 	feat_all <- scale(feat_all, center=mean_feat, scale=sd_feat)
 	feat_train <- scale(feat_train, center=mean_feat, scale=sd_feat)
 	feat_test <- scale(feat_test, center=mean_feat, scale=sd_feat)

 	#######################################
	## Fit Genomic Annotation Model (GAM)
	#######################################
	nfolds <- 5

	#gam_data <- logistic_regression_genomic_annotation_model_cv(feat_train, binary_outliers_train, nfolds, lambda_costs, lambda_init)
	#saveRDS(gam_data,"gam.RDS")
	#gam_data <- readRDS("gam.RDS")
	print(paste0(nfolds,"-fold cross validation on GAM yielded optimal lambda of ", gam_data$lambda))

	gam_posterior_test_obj <- update_independent_marginal_probabilities_exact_inference_cpp(feat_test, binary_outliers_test1, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta, matrix(0,2,2), matrix(0,2,2), number_of_dimensions, choose(number_of_dimensions, 2), FALSE)
	gam_test_posteriors <- gam_posterior_test_obj$probability


	#######################################
	### Initialize phi using genomic annotation model
	#######################################
	gam_posterior_train_obj <- update_independent_marginal_probabilities_exact_inference_cpp(feat_train, binary_outliers_test1, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta, matrix(0,2,2), matrix(0,2,2), number_of_dimensions, choose(number_of_dimensions, 2), FALSE)
	gam_train_posteriors <- gam_posterior_train_obj$probability

	if (phi_method == "sample_size") {
		pseudoc <- .001*as.vector(colSums(!is.na(discrete_outliers_train)))
	} else if (phi_method == "marginal") {
		pseudoc <- NA
	}
	#phi_init <- map_phi_initialization(discrete_outliers_train, gam_train_posteriors, number_of_dimensions, pseudoc)

	phi_init <- map_phi_initialization(discrete_outliers_train, gam_train_posteriors, number_of_dimensions, pseudoc)
	
 	#######################################
	## Fit Watershed Model (using training data)
	#######################################
	lambda_singleton <- 0
  	lambda_pair <- lambda_pair_init
  	lambda <- lambda_init

  	# watershed_model <- integratedEM(feat_train, discrete_outliers_train, phi_init, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta, pseudoc, lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method, independent_variables, vi_step_size, vi_threshold)

 	#######################################
	## Get test data watershed posterior probabilities
	#######################################
 	posterior_info_test <- update_marginal_posterior_probabilities(feat_test, discrete_outliers_test1, watershed_model)
  	posterior_prob_test <- posterior_info_test$probability  # Marginal posteriors
  	posterior_pairwise_prob_test <- posterior_info_test$probability_pairwise  # Pairwise posteriors



 	#######################################
	# Extract ROC curves and precision recall curves for test set (in each dimension seperately) using:
	#### 1. Watershed
	#### 2. GAM
	#### 3. RNA-only
	#######################################
	dimension_labels <- colnames(data_input$outliers_binary)
	median_tissue_river_model <- readRDS("/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/watershed_three_class_roc/fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_30_inference_exact_independent_true_roc_object2.rds")

	roc_object_across_dimensions <- compute_roc_across_dimensions(number_of_dimensions, dimension_labels, posterior_prob_test, real_valued_outliers_test1, gam_test_posteriors, binary_outliers_test2,median_tissue_river_model$watershed_predictions, outlier_type)

  	return(list(gam_params=gam_data, model_params=watershed_model,roc=roc_object_across_dimensions))

}



make_theta_pair_heatmap <- function(theta_pair, number_of_dimensions,tissue_names, output_file) {
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
    heatmap <- heatmap + labs(x = "Tissue", y = "Tissue",fill="Theta pair")

    ggsave(heatmap, file=output_file,width=7.2, height=6.0, units="in")
}


make_expression_overlap_heatmap <- function(pvalues, number_of_dimensions,tissue_names, output_file) {
	options(bitmapType = 'cairo', device = 'pdf')

	# Convert theta_pair vector into matrix of number_of_tissuesXnumber_of_tissues
	theta_pair_mat = matrix(0, number_of_dimensions, number_of_dimensions)
	dimension_counter = 1
	for (dimension1 in 1:number_of_dimensions) {
		for (dimension2 in dimension1:number_of_dimensions) {
			if (dimension1 != dimension2) {
				theta_pair_mat[dimension1, dimension2] = sum(!is.na(pvalues[,dimension1]) & !is.na(pvalues[,dimension2]))
				theta_pair_mat[dimension2, dimension1] = sum(!is.na(pvalues[,dimension1]) & !is.na(pvalues[,dimension2]))
				dimension_counter = dimension_counter + 1
			}
		}
	}
	# Cluster tissues based on similarity of theta_pairs
	order <- hclust( dist(scale(theta_pair_mat), method = "euclidean"), method = "ward.D" )$order

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
    heatmap <- heatmap + labs(x = "Tissue", y = "Tissue",fill="Overlap Counts")

    ggsave(heatmap, file=output_file,width=7.2, height=6.0, units="in")
}



get_tissue_names <- function(roc_object, number_of_dimensions) {
	tissue_names <- c()
	for (tissue_num in 1:number_of_dimensions) {
		tissue_name <- strsplit(roc_object$roc[[tissue_num]]$name, "_total_expression")[[1]][1]
		tissue_names <- c(tissue_names, tissue_name)
	}
	return(tissue_names)
}

pr_curve_in_one_tissue <- function(tissue_name, roc_vi, roc_independent, output_file) {
	options(bitmapType = 'cairo', device = 'pdf')

	precision <- c()
	recall <- c()
	prediction_type <- c()

	precision <- c(precision, roc_vi$watershed_precision)
	recall <- c(recall, roc_vi$watershed_recall)
	prediction_type <- c(prediction_type, rep("watershed", length(roc_vi$watershed_precision)))

	precision <- c(precision, roc_independent$watershed_precision)
	recall <- c(recall, roc_independent$watershed_recall)
	prediction_type <- c(prediction_type, rep("river", length(roc_independent$watershed_precision)))	

	df <- data.frame(precision, recall, prediction_type=factor(prediction_type, levels=c("watershed","river")))
  

  	plotter <- ggplot(data=df, aes(x=recall, y=precision, colour=prediction_type)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="",linetype="",title=tissue_name) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="right") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	ggsave(plotter, file=output_file,width = 19,height=11,units="cm")


}

make_tbt_auc_lolipop_plot <- function(roc_object_vi, roc_object_exact, tissue_names, output_file) {
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

	df <- data.frame(auc_watershed = auc_vi, auc_river=auc_exact, tissue=tissue_names, tissues_position=1:length(tissue_names))
	plotter <- ggplot(df) +
  			   geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_vi, yend=auc_exact), color="grey") +
  			   geom_point( aes(x=tissues_position, y=auc_vi), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
               geom_point( aes(x=tissues_position, y=auc_exact), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
               coord_flip()+
                xlab("") +
               ylab("AUC (PR)") + 
               scale_x_continuous(breaks=1:49,labels=rev(tissue_names))
    ggsave(plotter, file=output_file, width=8.8, height=7.0, units="in")
}

make_tbt_auc_lolipop_plot_with_median_river <- function(roc_object, tissue_names, output_file) {
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

	df <- data.frame(auc_watershed = auc_vi, auc_river=auc_exact, tissue=tissue_names, tissues_position=1:length(tissue_names))
	plotter <- ggplot(df) +
  			   geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_vi, yend=auc_exact), color="grey") +
  			   geom_point( aes(x=tissues_position, y=auc_vi), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
               geom_point( aes(x=tissues_position, y=auc_exact), color="dodgerblue4", size=3 ) +
               coord_flip()+
                xlab("") +
               ylab("AUC (PR)") + 
               scale_x_continuous(breaks=1:49,labels=rev(tissue_names))
    ggsave(plotter, file=output_file, width=8.8, height=7.0, units="in")
}



#########################################
# Command line arguments
#########################################
input_file <- args[1]  # Watershed input file
output_stem <- args[2]  # Stem to save all output files
number_of_dimensions <- as.numeric(args[3])  # Dimensionality of space
pseudoc <- as.numeric(args[4])  # Prior specification for P(E|Z)
n2_pair_pvalue_fraction <- as.numeric(args[5])  # For N2 Pair cross-validation, pick pvalue threshold for each outlier dimension such that this fraction of cases are positive examples
binary_pvalue_threshold <- as.numeric(args[6])  # Pvalue threshold to call binary outliers for genomic annotation model
phi_method <- args[7]
lambda_init <- as.numeric(args[8])
lambda_pair_init <- as.numeric(args[9])
independent_variables <- args[10]
inference_method <- args[11]
outlier_type <- args[12]


#####################
# Parameters
#####################
#lambda_costs <- c(.1,.01,1e-3, 1e-4)
lambda_costs <- c(.1,.01,1e-3,1e-4)
vi_step_size=.8
vi_threshold=1e-8

#######################################
## Load in data
#######################################
data_input <- load_watershed_data(input_file, number_of_dimensions, n2_pair_pvalue_fraction, binary_pvalue_threshold)





#######################################
## Run Watershed model
output_root <- paste0(output_stem,"_inference_", inference_method, "_independent_", independent_variables)
roc_object <- roc_analysis(data_input, number_of_dimensions, lambda_costs, pseudoc, inference_method, independent_variables, vi_step_size, vi_threshold, phi_method, lambda_init, lambda_pair_init, outlier_type, paste0(output_root, "_roc_object2.rds"))
saveRDS(roc_object, paste0(output_root, "_roc_object3.rds"))



#tissue_names <- get_tissue_names(roc_object_vi, number_of_dimensions)



######################################
# Visualize theta-pair heatmap to look at correlation structure across tissues
######################################
#output_file <- paste0(output_stem, "tissue_by_tissue_theta_pair_heatmap.pdf")
#make_theta_pair_heatmap(roc_object_vi$model_params$theta_pair, number_of_dimensions, tissue_names, output_file)

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
#output_file <- paste0(output_stem, "tissue_by_tissue_pr_auc_between_river_watershed_lolipop.pdf")
#make_tbt_auc_lolipop_plot(roc_object_vi$roc, roc_object_independent$roc, tissue_names, output_file)

######################################
# Visualize tbt pr auc curves between watershed and river
######################################
#output_file <- paste0(output_stem, "tissue_by_tissue_pr_auc_between_median_river_watershed_lolipop.pdf")
#make_tbt_auc_lolipop_plot_with_median_river(roc_object_vi$roc, tissue_names, output_file)

######################################
# Visualize num overlap heatmap to look at missingness structure across tissues
######################################
#output_file <- paste0(output_stem, "tissue_by_tissue_expression_overlap_heatmap.pdf")
#make_expression_overlap_heatmap(data_input$outlier_pvalues, number_of_dimensions, tissue_names, output_file)

if (FALSE) {
######################################
# Make Precision recall for each tissue
######################################
for (tissue_num in 1:(number_of_dimensions+1)) {
	tissue_name <- tissue_names[tissue_num]
	tissue_roc_vi <- roc_object_vi$roc[[tissue_num]]$evaROC
	tissue_roc_independent <- roc_object_independent$roc[[tissue_num]]$evaROC
	auc_vi <- tissue_roc_vi$watershed_pr_auc
	auc_independent <- tissue_roc_independent$watershed_pr_auc
	print(paste0(tissue_name, ": ", auc_vi, " ", auc_independent))
	#output_file <- paste0(output_stem, tissue_name, "_precision_recall_curves.pdf")
	#pr_curve_in_one_tissue(tissue_name, tissue_roc_vi, tissue_roc_independent, output_file)
}
}
#######################################
## Visualize precision-recall curves for river, GAM, watershed-exact comparison and all three outlier types (te, splice, ase)
#######################################
# plot_pr_gam_river_watershed_comparison_curve(roc_object_exact$roc, roc_object_independent$roc, number_of_dimensions, paste0(output_stem, "_watershed_exact_river_gam_comparison_pr.pdf"))
# print("DONE")

#######################################
## Visualize precision-recall curves for river, watershed-vi, watershed-exact comparison and all three outlier types (te, splice, ase)
#######################################
# plot_pr_river_watershed_comparison_curve(roc_object_exact$roc, roc_object_vi$roc, roc_object_independent$roc, number_of_dimensions, paste0(output_stem, "_watershed_exact_watershed_vi_river_comparison_pr.pdf"))



#######################################
## Visualize Confusion matrix for both RIVER and Watershed (exact and vi)
#######################################
# visualize_river_and_watershed_confusion_matrices(roc_object_exact$confusion, roc_object_vi$confusion, roc_object_independent$confusion, paste0(output_stem, "_confusion_heatmap_river_watershed_exact_watershed_vi_comparison.pdf"))









































