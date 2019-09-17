args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(RColorBrewer)
library(ggplot2)
options(bitmapType = 'cairo', device = 'pdf')


gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

make_outlier_type_specific_gam_posterior_comparison_histogram <- function(gene_level_posteriors, variant_level_posteriors, title) {
	posteriors <- c((gene_level_posteriors), (variant_level_posteriors))
	typer <- c(rep("gene", length(gene_level_posteriors)), rep("variant", length(variant_level_posteriors)))

	df <- data.frame(posterior=posteriors, type=factor(typer))

	plotter <- ggplot(df, aes(x=posterior, color=type, fill=type)) +
		geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
		scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
		scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
		scale_x_continuous(trans = "log") +
		labs(title=title,x="GAM Posterior", y = "Density", fill="", color="")

	return(plotter)
}

gam_posterior_comparison_histogram <- function(gene_level_posteriors, variant_level_posteriors) {
	splicing_histogram <- make_outlier_type_specific_gam_posterior_comparison_histogram(gene_level_posteriors$splicing_gam_posterior, variant_level_posteriors$splicing_gam_posterior, "Splicing")
	ase_histogram <- make_outlier_type_specific_gam_posterior_comparison_histogram(gene_level_posteriors$ase_gam_posterior, variant_level_posteriors$ase_gam_posterior, "ASE")
	te_histogram <- make_outlier_type_specific_gam_posterior_comparison_histogram(gene_level_posteriors$total_expression_gam_posterior, variant_level_posteriors$total_expression_gam_posterior, "TE")

	legend <- get_legend(ase_histogram + theme(legend.position="bottom"))

	combined <- plot_grid(ase_histogram + theme(legend.position="none"), splicing_histogram+ theme(legend.position="none"),te_histogram+ theme(legend.position="none"), ncol=1)

	return(plot_grid(combined, legend, ncol=1, rel_heights=c(1,.05)))
}

make_outlier_type_specific_gam_posterior_comparison_stripchart <- function(gene_level_posteriors, variant_level_posteriors, title) {
	posteriors <- c((gene_level_posteriors), (variant_level_posteriors))
	typer <- c(rep("gene", length(gene_level_posteriors)), rep("variant", length(variant_level_posteriors)))

	df <- data.frame(posterior=posteriors, type=factor(typer))

	plotter <- ggplot(df, aes(x=type, y=posterior, color=type)) + 
		geom_jitter(position=position_jitter(0.2)) +
		scale_color_manual(values=c("#E69F00", "#56B4E9")) + 
		scale_y_continuous(breaks=c(0,.25, .5, .75,1)) +
		ylim(0,1.0) +
		labs(x="", y="GAM Posterior", colour="", title=title)

	return(plotter)
}

make_outlier_type_specific_gam_watershed_scatter <- function(gam_posterior, watershed_posterior,title, coloring) {
	df <- data.frame(gam=gam_posterior, watershed=watershed_posterior)

	plotter <- ggplot(df, aes(x=gam, y=watershed)) + geom_point(color=coloring) +
		scale_y_continuous(breaks=c(0,.25, .5, .75,1)) +
		ylim(0,1.0) +
		scale_x_continuous(breaks=c(0,.25, .5, .75,1)) +
		xlim(0,1.0) +
		labs(x="GAM Posterior", y="Watershed Posterior", colour="", title=title)
	return(plotter)
}

gam_posterior_comparison_stripchart <- function(gene_level_posteriors, variant_level_posteriors) {
	splicing_stripchart <- make_outlier_type_specific_gam_posterior_comparison_stripchart(gene_level_posteriors$splicing_gam_posterior, variant_level_posteriors$splicing_gam_posterior, "Splicing")
	ase_stripchart <- make_outlier_type_specific_gam_posterior_comparison_stripchart(gene_level_posteriors$ase_gam_posterior, variant_level_posteriors$ase_gam_posterior, "ASE")
	te_stripchart <- make_outlier_type_specific_gam_posterior_comparison_stripchart(gene_level_posteriors$total_expression_gam_posterior, variant_level_posteriors$total_expression_gam_posterior, "TE")

	legend <- get_legend(ase_stripchart + theme(legend.position="bottom"))

	combined <- plot_grid(ase_stripchart + theme(legend.position="none"), splicing_stripchart+ theme(legend.position="none"),te_stripchart+ theme(legend.position="none"), ncol=3)

	return(plot_grid(combined, legend, ncol=1, rel_heights=c(1,.08)))
}

gam_watershed_comparison_scatter <- function(posteriors) {
	splicing_scatter <- make_outlier_type_specific_gam_watershed_scatter(posteriors$splicing_gam_posterior, posteriors$splicing_watershed_posterior, "Splicing", "#0D324D")
	ase_scatter <- make_outlier_type_specific_gam_watershed_scatter(posteriors$ase_gam_posterior, posteriors$ase_watershed_posterior, "ASE", "#7F5A83")
	te_scatter <- make_outlier_type_specific_gam_watershed_scatter(posteriors$total_expression_gam_posterior, posteriors$total_expression_watershed_posterior, "TE", "#BFCDE0")
	combined <- plot_grid(ase_scatter + theme(legend.position="none"), splicing_scatter+ theme(legend.position="none"),te_scatter+ theme(legend.position="none"), ncol=1)


	return(combined)
}



unsupervised_learning_input_dir <- args[1]
watershed_3_class_run_dir <- args[2]


###################
# File names
###################
gene_level_genomic_annotation_file <- paste0(unsupervised_learning_input_dir, "all_availibile_samples_features_filter_partially_observed_expression_short.txt")
variant_level_genomic_annotation_file <- paste0(unsupervised_learning_input_dir, "all_availibile_samples_variant_level_features_filter_partially_observed_expression_ukb_short.txt")

gene_level_ukb_posterior_file <- paste0(watershed_3_class_run_dir, "fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_gene_level_posteriors_ukb_short.txt")
variant_level_ukb_posterior_file <- paste0(watershed_3_class_run_dir, "fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors_ukb_short.txt")

gene_level_posterior_file <- paste0(watershed_3_class_run_dir, "fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_gene_level_posteriors_short.txt")
variant_level_posterior_file <- paste0(watershed_3_class_run_dir, "fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors_short.txt")


###################
# Load in data
###################
gene_level_posteriors <- read.table(gene_level_posterior_file, header=TRUE)
variant_level_posteriors <- read.table(variant_level_posterior_file, header=TRUE)

gene_level_ukb_posteriors <- read.table(gene_level_ukb_posterior_file, header=TRUE)
variant_level_ukb_posteriors <- read.table(variant_level_ukb_posterior_file, header=TRUE)


###################
# Visualizations
###################
# Comparison of GAM Posteriors to watershed posteriors
scatter_plot <- gam_watershed_comparison_scatter(variant_level_ukb_posteriors)
output_file <- paste0(watershed_3_class_run_dir, "gam_watershed_posterior_scatter_ukb_variant_level.pdf")
ggsave(scatter_plot, file=output_file, width=7.2,height=9,units="in")

if (FALSE) {
# Compare gene level vs variant level GAM posteriors with stripchart
stripchart_plot <- gam_posterior_comparison_stripchart(gene_level_posteriors, variant_level_posteriors)
output_file <- paste0(watershed_3_class_run_dir, "gam_posterior_stripchart_comparison.pdf")
ggsave(stripchart_plot, file=output_file, width=7.2,height=3.5,units="in")

# Compare gene level vs variant level GAM posteriors with stripchart
stripchart_plot <- gam_posterior_comparison_stripchart(gene_level_ukb_posteriors, variant_level_ukb_posteriors)
output_file <- paste0(watershed_3_class_run_dir, "gam_posterior_ukb_stripchart_comparison.pdf")
ggsave(stripchart_plot, file=output_file, width=7.2,height=3.5,units="in")
}
