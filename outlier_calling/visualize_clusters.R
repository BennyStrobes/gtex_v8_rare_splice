args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)


get_number_of_clusters <- function(jxn_names) {
	clusters <- c()
	for (itera in 1:length(jxn_names)) {
		jxn_name <- jxn_names[itera]
		cluster_id <- strsplit(jxn_name,split=":")[[1]][4]
		clusters <- c(clusters, cluster_id)
	}

	return(length(unique(clusters)))
}

get_number_of_genes <- function(jxn_names) {
	gene_list <- c()
	for (itera in 1:length(jxn_names)) {
		jxn_name <- jxn_names[itera]
		gene_string <- strsplit(jxn_name,split=":")[[1]][5]
		gene_info <- strsplit(gene_string,split=",")[[1]]
		for (itera in 1:length(gene_info)) {
			gene_list <- c(gene_list, gene_info[itera])
		}
	}

	return(length(unique(gene_list)))
}


get_number_of_jxns_per_cluster <- function(jxn_names) {
	clusters <- c()
	for (itera in 1:length(jxn_names)) {
		jxn_name <- jxn_names[itera]
		cluster_id <- strsplit(jxn_name,split=":")[[1]][4]
		clusters <- c(clusters, cluster_id)
	}
	return(count(clusters)$freq)
}

barplot_showing_number_of_clusters_per_tissue <- function(tissue_names, filtered_cluster_dir, output_file) {
	num_cluster_vec <- c()

	# Loop through tissues
	for (itera in 1:length(tissue_names)) {
		tissue_name <- tissue_names[itera]
		print(tissue_name)
		tissue_file_name <- paste0(filtered_cluster_dir, tissue_name, "_filtered_jxns_cross_tissue_clusters_gene_mapped.txt")
		# Load in data for this tissue 
		tissue_jxn_data <- read.table(tissue_file_name, header=TRUE)
		jxn_names <- as.character(tissue_jxn_data[,1])
		num_clusters <- get_number_of_clusters(jxn_names)
		num_cluster_vec <- c(num_cluster_vec, num_clusters)
	}
	tissue_names = factor(tissue_names, levels=rev(as.character(tissue_names)))

	df <- data.frame(tissue=tissue_names, number_of_clusters=num_cluster_vec)


	bar_plot <- ggplot(data=df, aes(x=tissue_names,y=number_of_clusters)) + geom_bar(stat="identity",fill="mediumpurple") + coord_flip() +
				labs(x = "Tissue", y = "Number of clusters") +
				gtex_v8_figure_theme()


	return(bar_plot)

}
gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

barplot_showing_number_of_genes_per_tissue <- function(tissue_names, filtered_cluster_dir, output_file) {
	num_cluster_vec <- c()

	# Loop through tissues
	for (itera in 1:length(tissue_names)) {
		tissue_name <- tissue_names[itera]
		print(tissue_name)
		tissue_file_name <- paste0(filtered_cluster_dir, tissue_name, "_filtered_jxns_cross_tissue_clusters_gene_mapped.txt")
		# Load in data for this tissue 
		tissue_jxn_data <- read.table(tissue_file_name, header=TRUE)
		jxn_names <- as.character(tissue_jxn_data[,1])
		num_clusters <- get_number_of_genes(jxn_names)
		num_cluster_vec <- c(num_cluster_vec, num_clusters)
	}

	tissue_names = factor(tissue_names, levels=rev(as.character(tissue_names)))
	df <- data.frame(tissue=tissue_names, number_of_clusters=num_cluster_vec)


    #melted_corr$X1 <- factor(melted_corr$X1, levels = rownames(correlation_matrix))


	bar_plot <- ggplot(data=df, aes(x=tissue_names,y=number_of_clusters)) + geom_bar(stat="identity",fill="darkcyan") + coord_flip() +
				labs(x = "Tissue", y = "Number of genes") +
				gtex_v8_figure_theme()
				#scale_y_discrete(limits = rev(levels(tissue_names))) +


	return(bar_plot)

}



barplot_showing_number_of_jxns_per_tissue <- function(tissue_names, filtered_cluster_dir, output_file) {
	num_jxns_vec <- c()

	# Loop through tissues
	for (itera in 1:length(tissue_names)) {
		tissue_name <- tissue_names[itera]
		print(tissue_name)
		tissue_file_name <- paste0(filtered_cluster_dir, tissue_name, "_filtered_jxns_cross_tissue_clusters_gene_mapped.txt")
		# Load in data for this tissue 
		tissue_jxn_data <- read.table(tissue_file_name, header=TRUE)
		num_jxns <- dim(tissue_jxn_data)[1]
		num_jxns_vec <- c(num_jxns_vec, num_jxns)
	}

	tissue_names = factor(tissue_names, levels=rev(as.character(tissue_names)))

	df <- data.frame(tissue=tissue_names, number_of_junctions=num_jxns_vec)


	bar_plot <- ggplot(data=df, aes(x=tissue_names,y=number_of_junctions)) + geom_bar(stat="identity",fill="steelblue3") + coord_flip() +
				labs(x = "Tissue", y = "Number of junctions") +
				gtex_v8_figure_theme()


	return(bar_plot)

}

boxplot_showing_number_of_jxns_per_cluster_per_tissue <- function(tissue_names, filtered_cluster_dir, output_file) {
	num_jxns_vec <- c()
	tissue_vec <- c()

	# Loop through tissues
	for (itera in 1:length(tissue_names)) {
		tissue_name <- tissue_names[itera]
		print(tissue_name)
		tissue_file_name <- paste0(filtered_cluster_dir, tissue_name, "_filtered_jxns_cross_tissue_clusters_gene_mapped.txt")
		# Load in data for this tissue 
		tissue_jxn_data <- read.table(tissue_file_name, header=TRUE)
		jxn_names <- as.character(tissue_jxn_data[,1])

		num_jxns_per_cluster <- get_number_of_jxns_per_cluster(jxn_names)

		num_jxns_vec <- c(num_jxns_vec, num_jxns_per_cluster)
		tissue_vec <- c(tissue_vec, rep(tissue_name, length(num_jxns_per_cluster)))

	}
	tissue_names = factor(tissue_vec, levels=rev(as.character(tissue_names)))

	df <- data.frame(tissue=tissue_vec, number_of_junctions_per_cluster=num_jxns_vec)
	box_plot <- ggplot(data=df, aes(x=tissue,y=number_of_junctions_per_cluster)) + geom_boxplot() + coord_flip() +
				labs(x = "Tissue", y = "Number of junctions/cluster") +
				theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 


	ggsave(box_plot, file=output_file,width = 15,height=12.5,units="cm")
}

histogram_showing_number_of_jxns_per_cluster_per_tissue <- function(tissue_name, filtered_cluster_dir, output_file) {
	num_jxns_vec <- c()
	tissue_vec <- c()

	tissue_file_name <- paste0(filtered_cluster_dir, tissue_name, "_filtered_jxns_cross_tissue_clusters_gene_mapped.txt")
	# Load in data for this tissue 
	tissue_jxn_data <- read.table(tissue_file_name, header=TRUE)
	jxn_names <- as.character(tissue_jxn_data[,1])

	num_jxns_per_cluster <- get_number_of_jxns_per_cluster(jxn_names)


	num_jxns_vec <- c(num_jxns_vec, num_jxns_per_cluster)
	tissue_vec <- c(tissue_vec, rep(tissue_name, length(num_jxns_per_cluster)))

	df <- data.frame(number_of_junctions_per_cluster=num_jxns_vec)
	histo <- ggplot(data=df, aes(x=number_of_junctions_per_cluster)) + geom_histogram(binwidth=2) +
				labs(x = "Number of junctions/cluster", title=tissue_name) +
				theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) 


	ggsave(histo, file=output_file,width = 15,height=12.5,units="cm")
}

histogram_showing_number_of_genes_per_cluster <- function(cluster_info_file, output_file) {
	num_genes_per_cluster <- c()
	data <- read.table(cluster_info_file, header=TRUE)
	genes <- as.character(data$genes)
	for (itera in 1:length(genes)) {
		gene <- genes[itera]
		num_genes <- length(strsplit(gene,split=",")[[1]])
		num_genes_per_cluster <- c(num_genes_per_cluster, num_genes)
	}
	df <- data.frame(num_genes_per_cluster=num_genes_per_cluster)

	histo <- ggplot(data=df, aes(x=num_genes_per_cluster)) + geom_histogram() +
		labs(x = "Number of genes/cluster") +
		theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) 
	ggsave(histo,file=output_file,width = 15,height=10.5,units="cm")

}

histogram_showing_fraction_of_expressed_samples <- function(tissue_name, min_reads, input_file, output_file) {
	data_table <- read.table(input_file,header=TRUE)[]
	row = dim(data_table)[1]
	col = dim(data_table)[2]
	#frac_expressed <- c()
	data_table <- data_table[,2:col]
	col = dim(data_table)[2]
	expressed <- data_table >= min_reads

	frac_expressed <- rowSums(expressed)/col

	df <- data.frame(fraction_expressed=frac_expressed)
	histo <- ggplot(data=df, aes(x=fraction_expressed)) + geom_histogram() +
		labs(x =paste0("Fraction of expressed (>= ", min_reads," read) samples / Cluster"), y="Number of clusters", title= paste0(tissue_name)) +
		theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), plot.title = element_text(hjust = 0.5), legend.title = element_text(size=10)) 
	ggsave(histo,file=output_file,width = 15,height=10.5,units="cm")

}	



tissue_names_file <- args[1]
filtered_cluster_dir <- args[2]
cluster_visualization_dir <- args[3]

options(bitmapType = 'cairo', device = 'pdf')


###################################################
# Histogram showing number of genes per cluster
###################################################
# Using exons for gene mapping
cluster_info_file <- paste0(filtered_cluster_dir, "cluster_info.txt")
output_file <- paste0(cluster_visualization_dir, "number_of_genes_per_cluster.pdf")
#histogram_showing_number_of_genes_per_cluster(cluster_info_file, output_file)

# Using whole gene for gene mapping
cluster_info_file <- paste0(filtered_cluster_dir, "cluster_info_full_gene_mapped.txt")
output_file <- paste0(cluster_visualization_dir, "number_of_genes_per_cluster_full_gene_mapped.pdf")
#histogram_showing_number_of_genes_per_cluster(cluster_info_file, output_file)


# Extract vector tissue names
tissue_names <- as.character(unlist(read.table(tissue_names_file,header=FALSE), use.names=FALSE))

if (FALSE) {
##################################################################
# Histogram showing fraction of samples w > $min_reads reads summed across all junctions
##################################################################

for (iter_num in 1:length(tissue_names)) {
	tissue_name <- tissue_names[iter_num]
	min_reads <- 3
	input_file <- paste0(filtered_cluster_dir, tissue_name, "_num_reads_per_cluster.txt")
	output_file <- paste0(cluster_visualization_dir, "fraction_of_", min_reads,"_expressed_samples_per_cluster_", tissue_name, "_histogram.pdf")
	histogram_showing_fraction_of_expressed_samples(tissue_name, min_reads, input_file, output_file)
}

}
###################################################
# Barplot showing number of genes per tissue
###################################################
print("A")

genes_per_tissue <- barplot_showing_number_of_genes_per_tissue(tissue_names, filtered_cluster_dir)
print("A")
###################################################
# Barplot showing number of clusters per tissue
###################################################
clusters_per_tissue <- barplot_showing_number_of_clusters_per_tissue(tissue_names, filtered_cluster_dir)
print("A")

###################################################
# Barplot showing number of jxns per tissue
###################################################
jxns_per_tissue <- barplot_showing_number_of_jxns_per_tissue(tissue_names, filtered_cluster_dir)
print("A")

combined_per_tissue <- plot_grid(jxns_per_tissue, clusters_per_tissue, genes_per_tissue, labels=c("A","B","C"), ncol=3)
ggsave(combined_per_tissue, file=paste0(cluster_visualization_dir, "combined_number_per_tissue_supplementary.pdf"),width=7.2, height=5, units="in")

combined_per_tissue <- plot_grid(jxns_per_tissue, clusters_per_tissue, genes_per_tissue, labels=c("A","B","C"), ncol=1)
ggsave(combined_per_tissue, file=paste0(cluster_visualization_dir, "combined_number_per_tissue_supplementary2.pdf"),width=7.2, height=7, units="in")


if (FALSE) {
###################################################
# Boxplot showing number of junctions per cluster per tissue
###################################################
output_file <- paste0(cluster_visualization_dir, "number_of_jxns_per_cluster_box_plot.pdf")
boxplot_showing_number_of_jxns_per_cluster_per_tissue(tissue_names, filtered_cluster_dir, output_file)


###################################################
# Histogram for specific tissue showing number of junctions per cluster
###################################################
for (iter_num in 1:length(tissue_names)) {
	tissue_name <- tissue_names[iter_num]
	output_file <- paste0(cluster_visualization_dir, "number_of_jxns_per_cluster_", tissue_name, "_histogram.pdf")
	histogram_showing_number_of_jxns_per_cluster_per_tissue(tissue_name, filtered_cluster_dir, output_file)
}
}

