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
options(bitmapType = 'cairo', device = 'pdf')

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}




make_gtex_posterior_vs_amish_pvalue_scatter <- function(merged_data, outlier_type) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_watershed_posterior

	corry <- cor.test(posterior, log_pval, method="spearman")
	pvalue = corry[3]
	correlation = corry[4]

	df <- data.frame(pvalue=log_pval, posterior=posterior)
	scatter <- ggplot(df, aes(x=pvalue, y=posterior)) +
 			 geom_point() +
 			 gtex_v8_figure_theme() + 
 			 labs(x=paste0("median -log10(pvalue) ", outlier_type, " in Amish cohort"),y=paste0("median GTEx Watershed ", outlier_type, " posterior"), title=paste0("Spearman rho: ", correlation, " / Spearman pvalue: ", pvalue))
 	return(scatter)

}

make_gtex_gam_posterior_vs_amish_pvalue_scatter <- function(merged_data, outlier_type) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_gam_posterior

	corry <- cor.test(posterior, log_pval, method="spearman")
	pvalue = corry[3]
	correlation = corry[4]

	df <- data.frame(pvalue=log_pval, posterior=posterior)
	scatter <- ggplot(df, aes(x=pvalue, y=posterior)) +
 			 geom_point() +
 			 gtex_v8_figure_theme() + 
 			 labs(x=paste0("median -log10(pvalue) ", outlier_type, " in Amish cohort"),y=paste0("median GTEx GAM ", outlier_type, " posterior"), title=paste0("Spearman rho: ", correlation, " / Spearman pvalue: ", pvalue))
 	return(scatter)

}

make_gtex_directional_watershed_posterior_vs_amish_pvalue_boxplot <- function(merged_data, outlier_type) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_watershed_posterior
	direction = merged_data$mode_gtex_direction


	high_watershed_indices = posterior > .8
	low_watershed_indices = posterior < .01

	high_watershed_pvalz = log_pval[high_watershed_indices]
	low_watershed_pvalz = log_pval[low_watershed_indices]

	high_watershed_direction = direction[high_watershed_indices]
	low_watershed_direction = direction[low_watershed_indices]


	pval <- c(high_watershed_pvalz, low_watershed_pvalz)
	type <- c(rep("Median GTEx Watershed posterior > .8", length(high_watershed_pvalz)), rep("Median GTEx Watershed posterior < .01", length(low_watershed_pvalz)))
	direction_arr <- c(as.character(high_watershed_direction), as.character(low_watershed_direction))


	df <- data.frame(pval=pval, type=factor(type, levels=c("Median GTEx Watershed posterior < .01", "Median GTEx Watershed posterior > .8")), direction=factor(direction_arr, levels=c("+", "-")))
	p <- ggplot(df, aes(x=type, y=pval, fill=direction)) + 
  		geom_boxplot() +
  		gtex_v8_figure_theme() + 
  		labs(title=outlier_type,fill="Expression outlier direction",x="", y=paste0(outlier_type, " median -log10(pvalue) in Amish cohort"))

  	return(p)

}

make_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot <- function(merged_data, outlier_type) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_watershed_posterior

	high_watershed_indices = posterior > .8
	# med_watershed_indices = posterior > .4 & posterior <=.8
	low_watershed_indices = posterior < .01

	high_watershed_pvalz = log_pval[high_watershed_indices]
	# med_watershed_pvalz = log_pval[med_watershed_indices]
	low_watershed_pvalz = log_pval[low_watershed_indices]

	# get gam thresholds
	num_high = length(high_watershed_pvalz)
	num_low = length(low_watershed_pvalz)

	high_threshold = sort(merged_data$median_gam_posterior)[(length(merged_data$median_gam_posterior) - num_high)]
	low_threshold = sort(merged_data$median_gam_posterior)[(num_low + 1)]
	high_gam_indices = merged_data$median_gam_posterior > high_threshold
	low_gam_indices = merged_data$median_gam_posterior < low_threshold

	high_gam_pvalz = log_pval[high_gam_indices]
	low_gam_pvalz = log_pval[low_gam_indices]


	pval <- c(high_watershed_pvalz, low_watershed_pvalz, high_gam_pvalz, low_gam_pvalz)
	type <- c(rep("Median GTEx posterior > .8", length(high_watershed_pvalz)), rep("Median GTEx posterior < .01", length(low_watershed_pvalz)), rep("Median GTEx posterior > .8", length(high_gam_pvalz)), rep("Median GTEx posterior < .01", length(low_gam_pvalz)))
	model <- c(rep("Watershed", length(high_watershed_pvalz)), rep("Watershed", length(low_watershed_pvalz)), rep("GAM", length(high_gam_pvalz)), rep("GAM", length(low_gam_pvalz)))

	#print(outlier_type)
	#print(wilcox.test(high_watershed_pvalz, low_watershed_pvalz))
	#print(wilcox.test(high_watershed_pvalz, high_gam_pvalz))
	#print(wilcox.test(high_gam_pvalz, low_gam_pvalz))

	df <- data.frame(pval=pval, type=factor(type, levels=c("Median GTEx posterior < .01", "Median GTEx posterior > .8")), model=factor(model, levels=c("Watershed", "GAM")))

	p <- ggplot(df, aes(x=type, y=pval, fill=model)) + 
  		geom_boxplot() +
  		scale_fill_manual(values=c("steelblue3", "firebrick4")) + 
  		gtex_v8_figure_theme() + 
  		labs(title=outlier_type,fill="",x="", y=paste0(outlier_type, " median -log10(p-value)\n in Amish cohort"))
  	return(p)

}
make_gtex_watershed_gam_posterior_vs_amish_pvalue_three_threshold_dotplot_boxplot <- function(merged_data, outlier_type) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_watershed_posterior

	high_watershed_indices = posterior > .8
	med_watershed_indices = posterior > .5 & posterior <=.8
	low_watershed_indices = posterior < .01

	high_watershed_pvalz = log_pval[high_watershed_indices]
	med_watershed_pvalz = log_pval[med_watershed_indices]
	low_watershed_pvalz = log_pval[low_watershed_indices]

	# get gam thresholds
	num_high = length(high_watershed_pvalz)
	num_med = length(med_watershed_pvalz)
	num_low = length(low_watershed_pvalz)



	high_threshold = sort(merged_data$median_gam_posterior)[(length(merged_data$median_gam_posterior) - num_high)]
	med_threshold = sort(merged_data$median_gam_posterior)[(length(merged_data$median_gam_posterior) - num_high - num_med)]
	low_threshold = sort(merged_data$median_gam_posterior)[(num_low + 1)]
	high_gam_indices = merged_data$median_gam_posterior > high_threshold
	med_gam_indices = merged_data$median_gam_posterior <= high_threshold & merged_data$median_gam_posterior > med_threshold
	low_gam_indices = merged_data$median_gam_posterior < low_threshold



	high_gam_pvalz = log_pval[high_gam_indices]
	med_gam_pvalz = log_pval[med_gam_indices]
	low_gam_pvalz = log_pval[low_gam_indices]


	print(outlier_type)
	print(wilcox.test(high_watershed_pvalz, high_gam_pvalz))

	high_watershed_label <- paste0("GTEx Watershed posterior > .8\nn=", num_high)
	med_watershed_label <- paste0(".5 > GTEx Watershed posterior <= .8\nn=", num_med)
	low_watershed_label <- paste0("GTEx Watershed posterior < .01\nn=", num_low)

	pval <- c(high_watershed_pvalz, med_watershed_pvalz, low_watershed_pvalz, high_gam_pvalz, med_gam_pvalz, low_gam_pvalz)
	type <- c(rep(high_watershed_label, length(high_watershed_pvalz)), rep(med_watershed_label, length(med_watershed_pvalz)), rep(low_watershed_label, length(low_watershed_pvalz)), rep(high_watershed_label, length(high_gam_pvalz)), rep(med_watershed_label, length(med_gam_pvalz)), rep(low_watershed_label, length(low_gam_pvalz)))
	model <- c(rep("Watershed", length(high_watershed_pvalz)),rep("Watershed", length(med_watershed_pvalz)), rep("Watershed", length(low_watershed_pvalz)), rep("GAM", length(high_gam_pvalz)),rep("GAM", length(med_gam_pvalz)), rep("GAM", length(low_gam_pvalz)))

	df <- data.frame(pval=pval, type=factor(type, levels=c(low_watershed_label, med_watershed_label, high_watershed_label)), model=factor(model, levels=c("Watershed", "GAM")))



	pval <- c(high_watershed_pvalz, med_watershed_pvalz, high_gam_pvalz, med_gam_pvalz)
	type <- c(rep(high_watershed_label, length(high_watershed_pvalz)), rep(med_watershed_label, length(med_watershed_pvalz)), rep(high_watershed_label, length(high_gam_pvalz)), rep(med_watershed_label, length(med_gam_pvalz)))
	model <- c(rep("Watershed", length(high_watershed_pvalz)),rep("Watershed", length(med_watershed_pvalz)), rep("GAM", length(high_gam_pvalz)),rep("GAM", length(med_gam_pvalz)))

	df2 <- data.frame(pval=pval, type=factor(type, levels=c(med_watershed_label, high_watershed_label)), model=factor(model, levels=c("Watershed", "GAM")))



	p <- ggplot() + 
  		geom_boxplot(data=df, aes(x=type, y=pval,fill=model), outlier.size=.1, alpha=0.22) +
  		geom_dotplot(data=df2, aes(x=type, y=pval, fill=model), binaxis='y', stackdir='center', position=position_dodge(width=.75)) +
  		scale_fill_manual(values=c("steelblue3", "firebrick4")) + 
  		scale_colour_manual(values=c("steelblue3", "firebrick4")) + 
  		gtex_v8_figure_theme() + 
  		labs(title=outlier_type,fill="", colour="",x="", y=paste0(outlier_type, " median -log10(p-value)\n in ASMAD cohort"))
  	return(p)	
}
make_gtex_watershed_gam_posterior_vs_amish_pvalue_three_threshold_boxplot <- function(merged_data, outlier_type) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_watershed_posterior

	high_watershed_indices = posterior > .8
	med_watershed_indices = posterior > .5 & posterior <=.8
	low_watershed_indices = posterior < .01

	high_watershed_pvalz = log_pval[high_watershed_indices]
	med_watershed_pvalz = log_pval[med_watershed_indices]
	low_watershed_pvalz = log_pval[low_watershed_indices]

	# get gam thresholds
	num_high = length(high_watershed_pvalz)
	num_med = length(med_watershed_pvalz)
	num_low = length(low_watershed_pvalz)



	high_threshold = sort(merged_data$median_gam_posterior)[(length(merged_data$median_gam_posterior) - num_high)]
	med_threshold = sort(merged_data$median_gam_posterior)[(length(merged_data$median_gam_posterior) - num_high - num_med)]
	low_threshold = sort(merged_data$median_gam_posterior)[(num_low + 1)]
	high_gam_indices = merged_data$median_gam_posterior > high_threshold
	med_gam_indices = merged_data$median_gam_posterior <= high_threshold & merged_data$median_gam_posterior > med_threshold
	low_gam_indices = merged_data$median_gam_posterior < low_threshold



	high_gam_pvalz = log_pval[high_gam_indices]
	med_gam_pvalz = log_pval[med_gam_indices]
	low_gam_pvalz = log_pval[low_gam_indices]


	high_watershed_label <- paste0("GTEx Watershed posterior > .8\nn=", num_high)
	med_watershed_label <- paste0(".4 > GTEx Watershed posterior <= .8\nn=", num_med)
	low_watershed_label <- paste0("GTEx Watershed posterior < .01\nn=", num_low)

	pval <- c(high_watershed_pvalz, med_watershed_pvalz, low_watershed_pvalz, high_gam_pvalz, med_gam_pvalz, low_gam_pvalz)
	type <- c(rep(high_watershed_label, length(high_watershed_pvalz)), rep(med_watershed_label, length(med_watershed_pvalz)), rep(low_watershed_label, length(low_watershed_pvalz)), rep(high_watershed_label, length(high_gam_pvalz)), rep(med_watershed_label, length(med_gam_pvalz)), rep(low_watershed_label, length(low_gam_pvalz)))
	model <- c(rep("Watershed", length(high_watershed_pvalz)),rep("Watershed", length(med_watershed_pvalz)), rep("Watershed", length(low_watershed_pvalz)), rep("GAM", length(high_gam_pvalz)),rep("GAM", length(med_gam_pvalz)), rep("GAM", length(low_gam_pvalz)))

	#print(outlier_type)
	#print(wilcox.test(med_watershed_pvalz, low_watershed_pvalz, alternative="greater"))
	#print(wilcox.test(high_watershed_pvalz, med_watershed_pvalz, alternative="greater"))
	#print(wilcox.test(high_watershed_pvalz, high_gam_pvalz))
	#print(wilcox.test(high_gam_pvalz, low_gam_pvalz))

	df <- data.frame(pval=pval, type=factor(type, levels=c(low_watershed_label, med_watershed_label, high_watershed_label)), model=factor(model, levels=c("Watershed", "GAM")))

	p <- ggplot(df, aes(x=type, y=pval, fill=model)) + 
  		geom_boxplot() +
  		scale_fill_manual(values=c("steelblue3", "firebrick4")) + 
  		gtex_v8_figure_theme() + 
  		labs(title=outlier_type,fill="",x="", y=paste0(outlier_type, " median -log10(p-value)\n in Amish cohort"))
  	return(p)

}

make_gtex_watershed_posterior_vs_amish_pvalue_boxplot_across_outlier_types <- function(splicing_merged_data, ase_merged_data, expression_merged_data) {
	# Load in Data
	splice_posterior = splicing_merged_data$median_watershed_posterior
	splice_log_pval = -log10(splicing_merged_data$median_amish_pvalue + 1e-6)
	ase_posterior = ase_merged_data$median_watershed_posterior
	ase_log_pval = -log10(ase_merged_data$median_amish_pvalue + 1e-6)
	te_posterior = expression_merged_data$median_watershed_posterior
	te_log_pval = -log10(expression_merged_data$median_amish_pvalue + 1e-6)


	# Indices corresponding to watershed bins
	high_watershed_splice_indices = splice_posterior > .8
	med_watershed_splice_indices = splice_posterior > .5 & splice_posterior <=.8
	low_watershed_splice_indices = splice_posterior < .01
	high_watershed_ase_indices = ase_posterior > .8
	med_watershed_ase_indices = ase_posterior > .5 & ase_posterior <=.8
	low_watershed_ase_indices = ase_posterior < .01
	high_watershed_te_indices = te_posterior > .8
	med_watershed_te_indices = te_posterior > .5 & te_posterior <=.8
	low_watershed_te_indices = te_posterior < .01

	# pvalues corresponding to above bins
	high_watershed_splice_pvalz = splice_log_pval[high_watershed_splice_indices]
	med_watershed_splice_pvalz = splice_log_pval[med_watershed_splice_indices]
	low_watershed_splice_pvalz = splice_log_pval[low_watershed_splice_indices]
	high_watershed_ase_pvalz = ase_log_pval[high_watershed_ase_indices]
	med_watershed_ase_pvalz = ase_log_pval[med_watershed_ase_indices]
	low_watershed_ase_pvalz = ase_log_pval[low_watershed_ase_indices]
	high_watershed_te_pvalz = te_log_pval[high_watershed_te_indices]
	med_watershed_te_pvalz = te_log_pval[med_watershed_te_indices]
	low_watershed_te_pvalz = te_log_pval[low_watershed_te_indices]

	# Concatenate into one vector
	pval <- c(high_watershed_splice_pvalz, med_watershed_splice_pvalz, low_watershed_splice_pvalz, high_watershed_ase_pvalz, med_watershed_ase_pvalz, low_watershed_ase_pvalz, high_watershed_te_pvalz, med_watershed_te_pvalz, low_watershed_te_pvalz)
	outlier_type <- c(rep("Splicing", length(high_watershed_splice_pvalz) + length(med_watershed_splice_pvalz) + length(low_watershed_splice_pvalz)), rep("ASE", length(high_watershed_ase_pvalz) + length(med_watershed_ase_pvalz) + length(low_watershed_ase_pvalz)), rep("Expression", length(high_watershed_te_pvalz) + length(med_watershed_te_pvalz) + length(low_watershed_te_pvalz)))
	type <- c()
	type <- c(type, rep("Median GTEx Watershed posterior > .8", l))

}


make_gtex_watershed_river_gam_posterior_vs_amish_pvalue_boxplot <- function(merged_data, outlier_type) {
	log_pval = -log10(merged_data$median_amish_pvalue + 1e-6)
	posterior = merged_data$median_watershed_posterior

	high_watershed_indices = posterior > .8
	# med_watershed_indices = posterior > .4 & posterior <=.8
	low_watershed_indices = posterior < .01

	high_watershed_pvalz = log_pval[high_watershed_indices]
	# med_watershed_pvalz = log_pval[med_watershed_indices]
	low_watershed_pvalz = log_pval[low_watershed_indices]

	# get gam thresholds
	num_high = length(high_watershed_pvalz)
	num_low = length(low_watershed_pvalz)

	high_gam_threshold = sort(merged_data$median_gam_posterior)[(length(merged_data$median_gam_posterior) - num_high)]
	low_gam_threshold = sort(merged_data$median_gam_posterior)[(num_low + 1)]
	high_gam_indices = merged_data$median_gam_posterior > high_gam_threshold
	low_gam_indices = merged_data$median_gam_posterior < low_gam_threshold

	high_gam_pvalz = log_pval[high_gam_indices]
	low_gam_pvalz = log_pval[low_gam_indices]


	high_river_threshold = sort(merged_data$median_river_posterior)[(length(merged_data$median_river_posterior) - num_high)]
	low_river_threshold = sort(merged_data$median_river_posterior)[(num_low + 1)]
	high_river_indices = merged_data$median_river_posterior > high_river_threshold
	low_river_indices = merged_data$median_river_posterior < low_river_threshold



	high_river_pvalz = log_pval[high_river_indices]
	#high_river_pvalz = log_pval[ merged_data$median_river_posterior > .8]
	low_river_pvalz = log_pval[low_river_indices]


	pval <- c(high_watershed_pvalz, low_watershed_pvalz, high_gam_pvalz, low_gam_pvalz, high_river_pvalz, low_river_pvalz)
	type <- c(rep("Median GTEx posterior > .8", length(high_watershed_pvalz)), rep("Median GTEx posterior < .01", length(low_watershed_pvalz)), rep("Median GTEx posterior > .8", length(high_gam_pvalz)), rep("Median GTEx posterior < .01", length(low_gam_pvalz)), rep("Median GTEx posterior > .8", length(high_river_pvalz)), rep("Median GTEx posterior < .01", length(low_river_pvalz)))
	model <- c(rep("Watershed", length(high_watershed_pvalz)), rep("Watershed", length(low_watershed_pvalz)), rep("GAM", length(high_gam_pvalz)), rep("GAM", length(low_gam_pvalz)), rep("RIVER", length(high_river_pvalz)), rep("RIVER", length(low_river_pvalz)))

	#print(wilcox.test(high_watershed_pvalz, low_watershed_pvalz))
	#print(wilcox.test(high_watershed_pvalz, high_gam_pvalz))
	#print(wilcox.test(high_gam_pvalz, low_gam_pvalz))

	df <- data.frame(pval=pval, type=factor(type, levels=c("Median GTEx posterior < .01", "Median GTEx posterior > .8")), model=factor(model, levels=c("Watershed", "RIVER", "GAM")))

	p <- ggplot(df, aes(x=type, y=pval, fill=model,color=model)) + 
  		geom_boxplot() +
  		scale_fill_manual(values=c("steelblue3", "White", "firebrick4")) + 
  		scale_color_manual(values=c("Black", "steelblue3", "Black")) + 
  		gtex_v8_figure_theme() + 
  		labs(title=outlier_type,fill="",color="",x="", y=paste0(outlier_type, " median -log10(pvalue) in Amish cohort"))
  	return(p)

}

make_gtex_watershed_river_gam_posterior_vs_amish_pvalue_pr_curve <- function(merged_data, outlier_type, pr_pvalue_threshold) {
	# Extract relevent data
	outlier_status = 1.0*(merged_data$median_amish_pvalue <= pr_pvalue_threshold)
	watershed_posteriors = merged_data$median_watershed_posterior
	river_posteriors = merged_data$median_river_posterior
	gam_posteriors = merged_data$median_gam_posterior

	# Make Precision-recall curves
	watershed_pr_object = pr.curve(scores.class0=watershed_posteriors[outlier_status==1.0], scores.class1=watershed_posteriors[outlier_status==0.0], curve=T)
	river_pr_object = pr.curve(scores.class0=river_posteriors[outlier_status==1.0], scores.class1=river_posteriors[outlier_status==0.0], curve=T)
	gam_pr_object = pr.curve(scores.class0=gam_posteriors[outlier_status==1.0], scores.class1=gam_posteriors[outlier_status==0.0], curve=T)

	print(watershed_pr_object$auc.integral)
	print(river_pr_object$auc.integral)
	print(gam_pr_object$auc.integral)
	precision = c(watershed_pr_object$curve[,2], river_pr_object$curve[,2], gam_pr_object$curve[,2])
	recall = c(watershed_pr_object$curve[,1], river_pr_object$curve[,1], gam_pr_object$curve[,1])
	model_type = c(rep("Watershed", length(watershed_pr_object$curve[,1])), rep("RIVER", length(river_pr_object$curve[,1])), rep("GAM", length(gam_pr_object$curve[,1])))

	df = data.frame(precision=precision, recall=recall, model_type=factor(model_type, levels=c("Watershed", "RIVER", "GAM")))
  	
  	plotter <- ggplot(data=df, aes(x=recall, y=precision, group=model_type)) + geom_line(aes(linetype=model_type, colour=model_type)) + 
                labs(x="Recall", y="Precision", group="", linetype="", colour="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_linetype_manual(values=c("solid", "dotted", "solid")) +
                gtex_v8_figure_theme() + 
                draw_label(outlier_type,x=.5,y=.95,size=8)
    return(plotter)
	#watershed_pr_auc=pr_obj$auc.integral,
	#watershed_recall=pr_obj$curve[,1],
	#watershed_precision=pr_obj$curve[,2],
	#pr.curve(scores.class0 = remove_na(gam_posteriors[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(gam_posteriors[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
}


merged_splicing_data_set_file <- args[1]
merged_ase_data_set_file <- args[2]
merged_expression_data_set_file <- args[3]
merged_rare_expression_data_set_file <- args[4]
output_dir <- args[5]


splicing_merged_data <- read.table(merged_splicing_data_set_file,header=TRUE)
ase_merged_data <- read.table(merged_ase_data_set_file,header=TRUE)
expression_merged_data <- read.table(merged_expression_data_set_file,header=TRUE)
#rare_expression_merged_data <- read.table(merged_rare_expression_data_set_file,header=TRUE)





###########################
# Splicing
############################

# make boxplot of amish splicing pvalues at different watershed thresholds
output_file <- paste0(output_dir, "splicing_gtex_watershed_gam_posterior_vs_amish_pvalue_dotplot_boxplot_3_thresholds.pdf")
splicing_3_threshold_dotplot_boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_three_threshold_dotplot_boxplot(splicing_merged_data, "Splicing")
ggsave(splicing_3_threshold_dotplot_boxplot, file=output_file, width=7.2, height=4,units="in")

# make boxplot of amish splicing pvalues at different watershed thresholds
output_file <- paste0(output_dir, "splicing_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot.pdf")
splicing_boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot(splicing_merged_data, "Splicing")
ggsave(splicing_boxplot, file=output_file, width=7.2, height=4,units="in")

# make boxplot of amish splicing pvalues at different watershed thresholds
output_file <- paste0(output_dir, "splicing_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot_3_thresholds.pdf")
splicing_3_threshold_boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_three_threshold_boxplot(splicing_merged_data, "Splicing")
ggsave(splicing_3_threshold_boxplot, file=output_file, width=7.2, height=4,units="in")


###########################
# ASE
############################


# make boxplot of amish splicing pvalues at different watershed thresholds
output_file <- paste0(output_dir, "ase_gtex_watershed_gam_posterior_vs_amish_pvalue_dotplot_boxplot_3_thresholds.pdf")
ase_3_threshold_dotplot_boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_three_threshold_dotplot_boxplot(ase_merged_data, "ASE")
ggsave(ase_3_threshold_dotplot_boxplot, file=output_file, width=7.2, height=4,units="in")

# make boxplot of amish splicing pvalues at different watershed thresholds
output_file <- paste0(output_dir, "ase_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot.pdf")
ase_boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot(ase_merged_data, "ASE")
ggsave(ase_boxplot, file=output_file, width=7.2, height=4,units="in")

# make boxplot of amish splicing pvalues at different watershed thresholds
output_file <- paste0(output_dir, "ase_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot_3_thresholds.pdf")
ase_3_threshold_boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_three_threshold_boxplot(ase_merged_data, "ASE")
ggsave(ase_3_threshold_boxplot, file=output_file, width=7.2, height=4,units="in")

###########################
# Expression
############################

# make boxplot of amish splicing pvalues at different watershed thresholds
output_file <- paste0(output_dir, "expression_gtex_watershed_gam_posterior_vs_amish_pvalue_dotplot_boxplot_3_thresholds.pdf")
te_3_threshold_dotplot_boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_three_threshold_dotplot_boxplot(expression_merged_data, "Expression")
ggsave(te_3_threshold_dotplot_boxplot, file=output_file, width=7.2, height=4,units="in")

# make boxplot of amish expression pvalues at different watershed thresholds
output_file <- paste0(output_dir, "expression_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot.pdf")
te_boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot(expression_merged_data, "Expression")
ggsave(te_boxplot, file=output_file, width=7.2, height=4,units="in")


# make boxplot of amish splicing pvalues at different watershed thresholds
output_file <- paste0(output_dir, "expression_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot_3_thresholds.pdf")
te_3_threshold_boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_three_threshold_boxplot(expression_merged_data, "Expression")
ggsave(te_3_threshold_boxplot, file=output_file, width=7.2, height=4,units="in")

# make boxplot of amish expression pvalues at different watershed thresholds seperated by direction
output_file <- paste0(output_dir, "expression_gtex_watershed_directional_posterior_vs_amish_pvalue_boxplot.pdf")
te_directional_boxplot <- make_gtex_directional_watershed_posterior_vs_amish_pvalue_boxplot(expression_merged_data, "Expression")
ggsave(te_directional_boxplot, file=output_file, width=7.2, height=4,units="in")




###########################
# Make Boxplot across three outlier types
###########################


output_file <- paste0(output_dir, "gtex_watershed_gam_posterior_vs_amish_pvalue_combined_boxplot.pdf")
legend <- get_legend(splicing_boxplot + theme(legend.position="bottom"))
combined_plots <- plot_grid(ase_boxplot + theme(legend.position="none"), splicing_boxplot+ theme(legend.position="none"), te_boxplot + theme(legend.position="none"), ncol=1)
boxplot_combined <- ggdraw() + draw_plot(combined_plots,0,.07,1,.9) + draw_plot(legend,.42,-0.43,1,1)
ggsave(boxplot_combined, file=output_file, width=7.2, height=7,units="in")


output_file <- paste0(output_dir, "gtex_watershed_gam_posterior_vs_amish_pvalue_combined_boxplot_3_thresholds.pdf")
legend <- get_legend(splicing_3_threshold_boxplot + theme(legend.position="bottom"))
combined_plots <- plot_grid(ase_3_threshold_boxplot + theme(legend.position="none"), splicing_3_threshold_boxplot+ theme(legend.position="none"), te_3_threshold_boxplot + theme(legend.position="none"), ncol=1)
boxplot_combined <- ggdraw() + draw_plot(combined_plots,0,.07,1,.9) + draw_plot(legend,.42,-0.43,1,1)
ggsave(boxplot_combined, file=output_file, width=7.2, height=7,units="in")



output_file <- paste0(output_dir, "gtex_watershed_gam_posterior_vs_amish_pvalue_combined_dotplot_boxplot_3_thresholds.pdf")
legend <- get_legend(splicing_3_threshold_dotplot_boxplot + theme(legend.position="bottom"))
combined_plots <- plot_grid(te_3_threshold_dotplot_boxplot + theme(legend.position="none"), ase_3_threshold_dotplot_boxplot + theme(legend.position="none"), splicing_3_threshold_dotplot_boxplot+ theme(legend.position="none"), ncol=1)
boxplot_combined <- ggdraw() + draw_plot(combined_plots,0,.07,1,.9) + draw_plot(legend,.42,-0.43,1,1)
ggsave(boxplot_combined, file=output_file, width=7.2, height=7,units="in")























#############
## OLD: NO longer used
########################


if (FALSE) {
pr_pvalue_threshold = .05

###########################
# make boxplot of amish rare-expression pvalues at different watershed thresholds
###########################
output_file <- paste0(output_dir, "rare_expression_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot(rare_expression_merged_data, "Expression")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")

output_file <- paste0(output_dir, "rare_expression_gtex_watershed_river_gam_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_watershed_river_gam_posterior_vs_amish_pvalue_boxplot(rare_expression_merged_data, "Expression")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")








###########################
# Splicing
############################

# make scatterplot of watershed posteriors vs amish splicing pvalue
###########################
output_file <- paste0(output_dir, "splicing_gtex_posterior_vs_amish_pvalue_scatter.pdf")
#scatter <- make_gtex_posterior_vs_amish_pvalue_scatter(splicing_merged_data, "Splicing")
#ggsave(scatter, file=output_file, width=7.2, height=4,units="in")

output_file <- paste0(output_dir, "splicing_gtex_gam_posterior_vs_amish_pvalue_scatter.pdf")
#scatter <- make_gtex_gam_posterior_vs_amish_pvalue_scatter(splicing_merged_data, "Splicing")
#ggsave(scatter, file=output_file, width=7.2, height=4,units="in")

###########################
# make precision recall curve
###########################
pr_curve_splice <- make_gtex_watershed_river_gam_posterior_vs_amish_pvalue_pr_curve(splicing_merged_data, "Splicing", pr_pvalue_threshold)


###########################
# make boxplot of amish splicing pvalues at different watershed thresholds
###########################
output_file <- paste0(output_dir, "splicing_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot(splicing_merged_data, "Splicing")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")

output_file <- paste0(output_dir, "splicing_gtex_watershed_river_gam_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_watershed_river_gam_posterior_vs_amish_pvalue_boxplot(splicing_merged_data, "Splicing")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")

###########################
# ASE
############################

###########################
# make scatterplot of watershed posteriors vs amish splicing pvalue
###########################
output_file <- paste0(output_dir, "ase_gtex_posterior_vs_amish_pvalue_scatter.pdf")
#scatter <- make_gtex_posterior_vs_amish_pvalue_scatter(ase_merged_data, "ASE")
#ggsave(scatter, file=output_file, width=7.2, height=4,units="in")

output_file <- paste0(output_dir, "ase_gtex_gam_posterior_vs_amish_pvalue_scatter.pdf")
#scatter <- make_gtex_gam_posterior_vs_amish_pvalue_scatter(ase_merged_data, "ASE")
#ggsave(scatter, file=output_file, width=7.2, height=4,units="in")
###########################
# make boxplot of amish splicing pvalues at different watershed thresholds
###########################
output_file <- paste0(output_dir, "ase_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot(ase_merged_data, "ASE")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")

output_file <- paste0(output_dir, "ase_gtex_watershed_river_gam_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_watershed_river_gam_posterior_vs_amish_pvalue_boxplot(ase_merged_data, "ASE")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")

###########################
# make precision recall curve
###########################
pr_curve_ase <- make_gtex_watershed_river_gam_posterior_vs_amish_pvalue_pr_curve(ase_merged_data, "ASE", pr_pvalue_threshold)


###########################
# Expression
############################

###########################
# make scatterplot of watershed posteriors vs amish splicing pvalue
###########################
output_file <- paste0(output_dir, "expression_gtex_posterior_vs_amish_pvalue_scatter.pdf")
#scatter <- make_gtex_posterior_vs_amish_pvalue_scatter(expression_merged_data, "Expression")
#ggsave(scatter, file=output_file, width=7.2, height=4,units="in")

output_file <- paste0(output_dir, "expression_gtex_gam_posterior_vs_amish_pvalue_scatter.pdf")
#scatter <- make_gtex_gam_posterior_vs_amish_pvalue_scatter(expression_merged_data, "Expression")
#ggsave(scatter, file=output_file, width=7.2, height=4,units="in")
###########################
# make boxplot of amish expression pvalues at different watershed thresholds
###########################
output_file <- paste0(output_dir, "expression_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_watershed_gam_posterior_vs_amish_pvalue_boxplot(expression_merged_data, "Expression")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")

output_file <- paste0(output_dir, "expression_gtex_watershed_river_gam_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_watershed_river_gam_posterior_vs_amish_pvalue_boxplot(expression_merged_data, "Expression")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")

###########################
# make boxplot of amish expression pvalues at different watershed thresholds seperated by direction
###########################
output_file <- paste0(output_dir, "expression_gtex_watershed_directional_posterior_vs_amish_pvalue_boxplot.pdf")
boxplot <- make_gtex_directional_watershed_posterior_vs_amish_pvalue_boxplot(expression_merged_data, "Expression")
ggsave(boxplot, file=output_file, width=7.2, height=4,units="in")


###########################
# make precision recall curve
###########################
pr_curve_te <- make_gtex_watershed_river_gam_posterior_vs_amish_pvalue_pr_curve(expression_merged_data, "Expression", pr_pvalue_threshold)











###########################
# Make combined precision recall curve across three outlier types (with cowplot)
###########################
output_file <- paste0(output_dir, "gtex_watershed_river_gam_posterior_vs_amish_pvalue_pr_curve.pdf")
legend <- get_legend(pr_curve_ase + theme(legend.position="bottom"))
combined_plots <- plot_grid(pr_curve_ase + theme(legend.position="none"), pr_curve_splice+ theme(legend.position="none"), pr_curve_te+ theme(legend.position="none"), legend, rel_heights=c(1,1,1,.03), nrow=1)
ggsave(combined_plots, file=output_file, width=11, height=3,units="in")















}
