args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)
options(bitmapType = 'cairo', device = 'pdf')


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


single_power_simulation_plot <- function(cluster_data, num_replicates, pvalue_test_column, pvalue_test_name) {

	num_reads_ordered = c(30, 50, 100, 500)
	data_subset = cluster_data[cluster_data$num_replicates==num_replicates,]

	df <- data.frame(efficiency=factor(data_subset$efficiency), pvalue=-log10(data_subset[,pvalue_test_column] + 1e-5), num_reads=factor(data_subset$num_reads, levels=num_reads_ordered))

	p<-ggplot(df, aes(x=efficiency, y=pvalue, fill=num_reads)) +
 		 geom_boxplot() + gtex_v8_figure_theme() + 
 		 theme(plot.title = element_text(hjust = 0.5)) +
 		 theme(plot.title=element_text(face="plain",size=10)) +
 		 theme(legend.position="bottom") +
 		 labs(x="Efficiency", y="-log10(pvalue + 1e-5)", title=paste0(num_replicates, " replicates / ", pvalue_test_name))
 	return(p)


}






plot_cluster_power_simulation <- function(cluster_data) {
	num_replicate_arr <- c(2,3,4)

	plotter_2_chi <- single_power_simulation_plot(cluster_data, 2, 5, "chi_squared")	
	plotter_3_chi <- single_power_simulation_plot(cluster_data, 3, 5, "chi_squared")
	plotter_4_chi <- single_power_simulation_plot(cluster_data, 4, 5, "chi_squared")


	plotter_2_fisher <- single_power_simulation_plot(cluster_data, 2, 6, "fisher_exact")	
	plotter_3_fisher <- single_power_simulation_plot(cluster_data, 3, 6, "fisher_exact")
	plotter_4_fisher <- single_power_simulation_plot(cluster_data, 4, 6, "fisher_exact")


	#combined <- plot_grid(plotter_2_chi, plotter_3_chi, plotter_4_chi, plotter_2_fisher, plotter_3_fisher, plotter_4_fisher, ncol=3)
	combined <- plot_grid(plotter_2_chi, plotter_3_chi, plotter_4_chi,ncol=3)

	return(combined)
}








directory <- args[1]



simulation_results_file <- paste0(directory, "simulation_results.txt")


data <- read.table(simulation_results_file, header=TRUE)

clusters <- as.character(unique(data$cluster_id))
for (cluster_iter in 1:length(clusters)) {

	line_cluster_id <- clusters[cluster_iter]

	cluster_data <- data[as.character(data$cluster_id)==line_cluster_id, ]

	cluster_output_file <- paste0(directory, line_cluster_id, "_power_simulation.pdf")

	plotter <- plot_cluster_power_simulation(cluster_data)

	ggsave(plotter, file=cluster_output_file, width=15.2, height=5,units="in")

}
