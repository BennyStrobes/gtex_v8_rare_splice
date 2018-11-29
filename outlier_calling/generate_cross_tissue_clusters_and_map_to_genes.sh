#!/bin/bash -l

#SBATCH
#SBATCH --time=3:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1



tissue_names_file="$1"
filtered_cluster_dir="$2"
gencode_gene_annotation_file="$3"
cluster_visualization_dir="$4"

if false; then
python generate_cross_tissue_clusters.py $tissue_names_file $filtered_cluster_dir
fi



Rscript visualize_clusters.R $tissue_names_file $filtered_cluster_dir $cluster_visualization_dir