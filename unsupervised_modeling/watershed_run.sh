#!/bin/bash -l

#SBATCH
#SBATCH --time=22:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1




unsupervised_learning_input_dir="$1"
watershed_run_dir="$2"
pseudocount="$3"
pvalue_threshold="$4"
gene_thresh="$5"

# Parameters!
number_of_dimensions="3"
inference_method="exact"
input_stem="fully_observed_merged_outliers_"$gene_thresh"_genes_intersection_between_te_ase_splicing"
input_file=$unsupervised_learning_input_dir$input_stem"_features_filter_N2_pairs.txt"
output_stem=$watershed_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_"$gene_thresh"_outlier_thresh_"$pvalue_threshold"_pseudocount_"$pseudocount"_"$inference_method"_inference"

echo $pseudocount
echo $pvalue_threshold
echo $gene_thresh
Rscript watershed_roc.R $pvalue_threshold $input_file $output_stem $number_of_dimensions $inference_method $pseudocount
