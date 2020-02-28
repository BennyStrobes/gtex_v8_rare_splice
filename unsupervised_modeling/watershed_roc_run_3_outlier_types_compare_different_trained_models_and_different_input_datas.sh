#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1

n2_pair_pvalue_fraction="$1"
binary_pvalue_threshold="$2"
unsupervised_learning_input_dir="$3"
watershed_3_class_roc_run_dir="$4"


# Parameters!
number_of_dimensions="3"
standard_input_file=$unsupervised_learning_input_dir"fully_observed_merged_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_no_tissue_anno_N2_pairs_3.txt"



##########################
# Version 1: Compare standard trained model to data with gene threshold of .05
###########################
model_stem=$watershed_3_class_roc_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_30_seed_3_"
n2_pairs_file=$unsupervised_learning_input_dir"fully_observed_merged_outliers_0.05_genes_intersection_between_te_ase_splicing_features_filter_no_tissue_anno_diff_N2_pairs_as_standard_450000_3.txt"
output_stem=$watershed_3_class_roc_run_dir"watershed_roc_3_class_model_data_comparison_standard_model_gene_0.05_intersection_data_"

watershed_model=$model_stem"inference_exact_independent_false_roc_object3.rds"
river_model=$model_stem"inference_exact_independent_true_roc_object3.rds"
Rscript watershed_roc_comparison_of_model_and_data.R $standard_input_file $n2_pairs_file $watershed_model $river_model $output_stem $number_of_dimensions $n2_pair_pvalue_fraction $binary_pvalue_threshold


##########################
# Version 2: Compare standard trained model to data with gene threshold of .1
###########################
model_stem=$watershed_3_class_roc_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_30_seed_3_"
n2_pairs_file=$unsupervised_learning_input_dir"fully_observed_merged_outliers_0.1_genes_intersection_between_te_ase_splicing_features_filter_no_tissue_anno_diff_N2_pairs_as_standard_450000_3.txt"
output_stem=$watershed_3_class_roc_run_dir"watershed_roc_3_class_model_data_comparison_standard_model_gene_0.1_intersection_data_"

watershed_model=$model_stem"inference_exact_independent_false_roc_object3.rds"
river_model=$model_stem"inference_exact_independent_true_roc_object3.rds"

Rscript watershed_roc_comparison_of_model_and_data.R $standard_input_file $n2_pairs_file $watershed_model $river_model $output_stem $number_of_dimensions $n2_pair_pvalue_fraction $binary_pvalue_threshold


##########################
# Version 3: Compare standard trained model to data with gene threshold of union of .01
###########################
model_stem=$watershed_3_class_roc_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_30_seed_3_"
n2_pairs_file=$unsupervised_learning_input_dir"fully_observed_merged_outliers_0.01_genes_union_between_te_ase_splicing_features_filter_no_tissue_anno_diff_N2_pairs_as_standard_450000_3.txt"
output_stem=$watershed_3_class_roc_run_dir"watershed_roc_3_class_model_data_comparison_standard_model_gene_0.01_union_data_"

watershed_model=$model_stem"inference_exact_independent_false_roc_object3.rds"
river_model=$model_stem"inference_exact_independent_true_roc_object3.rds"

Rscript watershed_roc_comparison_of_model_and_data.R $standard_input_file $n2_pairs_file $watershed_model $river_model $output_stem $number_of_dimensions $n2_pair_pvalue_fraction $binary_pvalue_threshold
