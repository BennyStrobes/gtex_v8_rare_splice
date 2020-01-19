#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


unsupervised_learning_input_dir="$1"
watershed_run_dir="$2"
n2_pair_pvalue_fraction="$3"
binary_pvalue_threshold="$4"
gene_thresh="$5"



# Parameters!
number_of_dimensions="3"

input_stem="fully_observed_merged_outliers_"$gene_thresh"_genes_intersection_between_te_ase_splicing"

input_file=$unsupervised_learning_input_dir$input_stem"_features_filter_no_tissue_anno_same_N2_pairs_as_standard_3.txt"
standard_input_file=$unsupervised_learning_input_dir"fully_observed_merged_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_no_tissue_anno_N2_pairs_3.txt"


output_stem=$watershed_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_"$gene_thresh"_n2_pair_outlier_fraction_"$n2_pair_pvalue_fraction"_binary_pvalue_threshold_"$binary_pvalue_threshold"_gene_theshold_comparison"

Rscript watershed_roc_3_outlier_types_gene_threshold_comparison.R $input_file $output_stem $number_of_dimensions $n2_pair_pvalue_fraction $binary_pvalue_threshold $standard_input_file

