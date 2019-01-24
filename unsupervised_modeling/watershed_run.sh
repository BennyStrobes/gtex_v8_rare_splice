#!/bin/bash -l

#SBATCH
#SBATCH --time=16:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1




unsupervised_learning_input_dir="$1"
watershed_run_dir="$2"


gene_thresh="0.01"



pvalue_threshold=".01"
stem="fully_observed_merged_outliers_"$gene_thresh"_genes_intersection_between_te_splicing"
input_file=$unsupervised_learning_input_dir$stem"_features_filter_N2_pairs.txt"
Rscript watershed.R $pvalue_threshold $input_file $stem $watershed_run_dir