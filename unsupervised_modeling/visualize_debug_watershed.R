args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(PRROC)
library(cowplot)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}



compare_posteriors <- function(old_posterior, new_posterior, ase_change, titler) {
	df <- data.frame(old_posterior=old_posterior, new_posterior=new_posterior, ase_change=factor(ase_change))
	plotter <- ggplot(df, aes(x=old_posterior, y=new_posterior, colour=ase_change)) + geom_point() +
		labs(x="old posterior", y="new posterior",title=titler) + 
		figure_theme() + 
		geom_abline()
	return(plotter)
}






watershed_3_class_roc_run_dir <- args[1]
watershed_3_class_score_run_dir <- args[2]
watershed_debug_visualization_dir <- args[3]


############################
# Variant level comparison
############################

variant_level_score_comparison_file <- paste0(watershed_3_class_score_run_dir, "variant_level_comparison.txt")

variant_level_comparison <- read.table(variant_level_score_comparison_file, header=TRUE)


splicing_watershed_scatter = compare_posteriors(variant_level_comparison$splicing_watershed_old, variant_level_comparison$splicing_watershed, variant_level_comparison$ase_change, "Watershed: Splicing")
ggsave(splicing_watershed_scatter, file=paste0(watershed_debug_visualization_dir, "variant_level_splicing_watershed_scatter.pdf"), width=7.2, height=4.6, units="in")


te_watershed_scatter = compare_posteriors(variant_level_comparison$te_watershed_old, variant_level_comparison$te_watershed, variant_level_comparison$ase_change, "Watershed: Expression")
ggsave(te_watershed_scatter, file=paste0(watershed_debug_visualization_dir, "variant_level_te_watershed_scatter.pdf"), width=7.2, height=4.6, units="in")

ase_watershed_scatter = compare_posteriors(variant_level_comparison$ase_watershed_old, variant_level_comparison$ase_watershed, variant_level_comparison$ase_change, "Watershed: ASE")
ggsave(ase_watershed_scatter, file=paste0(watershed_debug_visualization_dir, "variant_level_ase_watershed_scatter.pdf"), width=7.2, height=4.6, units="in")



splicing_river_scatter = compare_posteriors(variant_level_comparison$splicing_river_old, variant_level_comparison$splicing_river, variant_level_comparison$ase_change, "RIVER: Splicing")
ggsave(splicing_river_scatter, file=paste0(watershed_debug_visualization_dir, "variant_level_splicing_river_scatter.pdf"), width=7.2, height=4.6, units="in")


te_river_scatter = compare_posteriors(variant_level_comparison$te_river_old, variant_level_comparison$te_river, variant_level_comparison$ase_change, "RIVER: Expression")
ggsave(te_river_scatter, file=paste0(watershed_debug_visualization_dir, "variant_level_te_river_scatter.pdf"), width=7.2, height=4.6, units="in")

ase_river_scatter = compare_posteriors(variant_level_comparison$ase_river_old, variant_level_comparison$ase_river, variant_level_comparison$ase_change, "RIVER: ASE")
ggsave(ase_river_scatter, file=paste0(watershed_debug_visualization_dir, "variant_level_ase_river_scatter.pdf"), width=7.2, height=4.6, units="in")


############################
# Gene level comparison
############################

gene_level_score_comparison_file <- paste0(watershed_3_class_score_run_dir, "gene_level_comparison.txt")

gene_level_comparison <- read.table(gene_level_score_comparison_file, header=TRUE)


splicing_watershed_scatter = compare_posteriors(gene_level_comparison$splicing_watershed_old, gene_level_comparison$splicing_watershed, gene_level_comparison$ase_change, "Watershed: Splicing")
ggsave(splicing_watershed_scatter, file=paste0(watershed_debug_visualization_dir, "gene_level_splicing_watershed_scatter.pdf"), width=7.2, height=4.6, units="in")


te_watershed_scatter = compare_posteriors(gene_level_comparison$te_watershed_old, gene_level_comparison$te_watershed, gene_level_comparison$ase_change, "Watershed: Expression")
ggsave(te_watershed_scatter, file=paste0(watershed_debug_visualization_dir, "gene_level_te_watershed_scatter.pdf"), width=7.2, height=4.6, units="in")

ase_watershed_scatter = compare_posteriors(gene_level_comparison$ase_watershed_old, gene_level_comparison$ase_watershed, gene_level_comparison$ase_change, "Watershed: ASE")
ggsave(ase_watershed_scatter, file=paste0(watershed_debug_visualization_dir, "gene_level_ase_watershed_scatter.pdf"), width=7.2, height=4.6, units="in")



splicing_river_scatter = compare_posteriors(gene_level_comparison$splicing_river_old, gene_level_comparison$splicing_river, gene_level_comparison$ase_change, "RIVER: Splicing")
ggsave(splicing_river_scatter, file=paste0(watershed_debug_visualization_dir, "gene_level_splicing_river_scatter.pdf"), width=7.2, height=4.6, units="in")


te_river_scatter = compare_posteriors(gene_level_comparison$te_river_old, gene_level_comparison$te_river, gene_level_comparison$ase_change, "RIVER: Expression")
ggsave(te_river_scatter, file=paste0(watershed_debug_visualization_dir, "gene_level_te_river_scatter.pdf"), width=7.2, height=4.6, units="in")

ase_river_scatter = compare_posteriors(gene_level_comparison$ase_river_old, gene_level_comparison$ase_river, gene_level_comparison$ase_change, "RIVER: ASE")
ggsave(ase_river_scatter, file=paste0(watershed_debug_visualization_dir, "gene_level_ase_river_scatter.pdf"), width=7.2, height=4.6, units="in")

