#!/bin/bash -l

#SBATCH
#SBATCH --time=17:00:00
#SBATCH --mem=10GB
#SBATCH --partition=shared
#SBATCH --nodes=1



tissue_names_file="$1"
filtered_cluster_dir="$2"
gencode_gene_annotation_file="$3"
cluster_visualization_dir="$4"
gene_list="$5"

if false; then
python generate_cross_tissue_clusters.py $tissue_names_file $filtered_cluster_dir
fi

python map_clusters_to_genes.py $tissue_names_file $filtered_cluster_dir $gencode_gene_annotation_file $gene_list

if false; then
Rscript visualize_clusters.R $tissue_names_file $filtered_cluster_dir $cluster_visualization_dir
fi