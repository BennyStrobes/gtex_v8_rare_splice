args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)


get_number_of_clusters <- function(jxn_names) {
	clusters <- c()
	for (itera in 1:length(jxn_names)) {
		jxn_name <- jxn_names[itera]
		cluster_id <- strsplit(jxn_name,split=":")[[1]][4]
		clusters <- c(clusters, cluster_id)
	}

	return(length(unique(clusters)))
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
		tissue_file_name <- paste0(filtered_cluster_dir, tissue_name, "_filtered_jxns_cross_tissue_clusters.txt")
		# Load in data for this tissue 
		tissue_jxn_data <- read.table(tissue_file_name, header=TRUE)
		jxn_names <- as.character(tissue_jxn_data[,1])
		num_clusters <- get_number_of_clusters(jxn_names)
		num_cluster_vec <- c(num_cluster_vec, num_clusters)
	}
	df <- data.frame(tissue=tissue_names, number_of_clusters=num_cluster_vec)


	bar_plot <- ggplot(data=df, aes(x=tissue_names,y=number_of_clusters)) + geom_bar(stat="identity",fill="mediumpurple") + coord_flip() +
				labs(x = "Tissue", y = "Number of clusters") +
				theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) 


	ggsave(bar_plot, file=output_file,width = 15,height=10.5,units="cm")


}



barplot_showing_number_of_jxns_per_tissue <- function(tissue_names, filtered_cluster_dir, output_file) {
	num_jxns_vec <- c()

	# Loop through tissues
	for (itera in 1:length(tissue_names)) {
		tissue_name <- tissue_names[itera]
		print(tissue_name)
		tissue_file_name <- paste0(filtered_cluster_dir, tissue_name, "_filtered_jxns_cross_tissue_clusters.txt")
		# Load in data for this tissue 
		tissue_jxn_data <- read.table(tissue_file_name, header=TRUE)
		num_jxns <- dim(tissue_jxn_data)[1]
		num_jxns_vec <- c(num_jxns_vec, num_jxns)
	}
	df <- data.frame(tissue=tissue_names, number_of_junctions=num_jxns_vec)


	bar_plot <- ggplot(data=df, aes(x=tissue_names,y=number_of_junctions)) + geom_bar(stat="identity",fill="steelblue3") + coord_flip() +
				labs(x = "Tissue", y = "Number of junctions") +
				theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) 


	ggsave(bar_plot, file=output_file,width = 15,height=10.5,units="cm")


}

boxplot_showing_number_of_jxns_per_cluster_per_tissue <- function(tissue_names, filtered_cluster_dir, output_file) {
	num_jxns_vec <- c()
	tissue_vec <- c()

	# Loop through tissues
	for (itera in 1:length(tissue_names)) {
		tissue_name <- tissue_names[itera]
		print(tissue_name)
		tissue_file_name <- paste0(filtered_cluster_dir, tissue_name, "_filtered_jxns_cross_tissue_clusters.txt")
		# Load in data for this tissue 
		tissue_jxn_data <- read.table(tissue_file_name, header=TRUE)
		jxn_names <- as.character(tissue_jxn_data[,1])

		num_jxns_per_cluster <- get_number_of_jxns_per_cluster(jxn_names)

		num_jxns_vec <- c(num_jxns_vec, num_jxns_per_cluster)
		tissue_vec <- c(tissue_vec, rep(tissue_name, length(num_jxns_per_cluster)))

	}
	df <- data.frame(tissue=tissue_vec, number_of_junctions_per_cluster=num_jxns_vec)
	box_plot <- ggplot(data=df, aes(x=tissue,y=number_of_junctions_per_cluster)) + geom_boxplot() + coord_flip() +
				labs(x = "Tissue", y = "Number of junctions/cluster") +
				theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) 


	ggsave(box_plot, file=output_file,width = 15,height=10.5,units="cm")
}

histogram_showing_number_of_jxns_per_cluster_per_tissue <- function(tissue_name, filtered_cluster_dir, output_file) {
	num_jxns_vec <- c()
	tissue_vec <- c()

	tissue_file_name <- paste0(filtered_cluster_dir, tissue_name, "_filtered_jxns_cross_tissue_clusters.txt")
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


	ggsave(histo, file=output_file,width = 15,height=10.5,units="cm")
}



tissue_names_file <- args[1]
filtered_cluster_dir <- args[2]
cluster_visualization_dir <- args[3]




# Extract vector tissue names
tissue_names <- as.character(unlist(read.table(tissue_names_file,header=FALSE), use.names=FALSE))[1:2]



###################################################
# Barplot showing number of clusters per tissue
###################################################
output_file <- paste0(cluster_visualization_dir, "number_of_clusters_bar_plot.png")
barplot_showing_number_of_clusters_per_tissue(tissue_names, filtered_cluster_dir, output_file)


###################################################
# Barplot showing number of jxns per tissue
###################################################
output_file <- paste0(cluster_visualization_dir, "number_of_jxns_bar_plot.png")
barplot_showing_number_of_jxns_per_tissue(tissue_names, filtered_cluster_dir, output_file)



###################################################
# Boxplot showing number of junctions per cluster per tissue
###################################################
output_file <- paste0(cluster_visualization_dir, "number_of_jxns_per_cluster_box_plot.png")
boxplot_showing_number_of_jxns_per_cluster_per_tissue(tissue_names, filtered_cluster_dir, output_file)


###################################################
# Histogram for specific tissue showing number of junctions per cluster
###################################################
tissue_name <- "Heart_Left_Ventricle"
output_file <- paste0(cluster_visualization_dir, "number_of_jxns_per_cluster_", tissue_name, "_histogram.png")
histogram_showing_number_of_jxns_per_cluster_per_tissue(tissue_name, filtered_cluster_dir, output_file)

tissue_name <- "Muscle_Skeletal"
output_file <- paste0(cluster_visualization_dir, "number_of_jxns_per_cluster_", tissue_name, "_histogram.png")
histogram_showing_number_of_jxns_per_cluster_per_tissue(tissue_name, filtered_cluster_dir, output_file)

tissue_name <- "Adipose_Subcutaneous"
output_file <- paste0(cluster_visualization_dir, "number_of_jxns_per_cluster_", tissue_name, "_histogram.png")
histogram_showing_number_of_jxns_per_cluster_per_tissue(tissue_name, filtered_cluster_dir, output_file)

