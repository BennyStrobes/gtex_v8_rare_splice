#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


unsupervised_learning_input_dir="$1"
watershed_run_dir="$2"
n2_pair_pvalue_fraction="$3"
binary_pvalue_threshold="$4"
input_file="$5"
output_stem="$6"




# Parameters!
number_of_dimensions="3"


standard_input_file=$unsupervised_learning_input_dir"fully_observed_merged_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_no_tissue_anno_N2_pairs_3.txt"


Rscript watershed_roc_3_outlier_types_comparison.R $input_file $output_stem $number_of_dimensions $n2_pair_pvalue_fraction $binary_pvalue_threshold $standard_input_file

