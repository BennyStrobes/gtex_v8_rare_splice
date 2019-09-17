#!/bin/bash -l

#SBATCH
#SBATCH --time=24:00:00
#SBATCH --mem=60GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1




unsupervised_learning_input_dir="$1"
watershed_run_dir="$2"
pseudocount="$3"
pvalue_fraction="$4"
binary_pvalue_threshold="$5"
watershed_3_class_roc_run_dir="$6"
# Parameters!
number_of_dimensions="3"
inference_method="exact"

gene_thresh="0.01"
input_stem="old_fully_observed_merged_outliers_"$gene_thresh"_genes_intersection_between_te_ase_splicing"
fully_observed_input_file=$unsupervised_learning_input_dir$input_stem"_features_filter_no_tissue_anno_N2_pairs.txt"


################
# Run analysis at the gene level
################
all_variants_input_file=$unsupervised_learning_input_dir"old_all_availibile_samples_features_filter_partially_observed_expression.txt"
output_stem=$watershed_run_dir"old_fully_observed_te_ase_splicing_outliers_gene_pvalue_"$gene_thresh"_outlier_fraction_"$pvalue_fraction"_pseudocount_"$pseudocount"_"$inference_method"_inference_apply_to_all_variants_gene_level"
Rscript watershed_score.R $pvalue_fraction $number_of_dimensions $inference_method $pseudocount $fully_observed_input_file $all_variants_input_file $output_stem $binary_pvalue_threshold $watershed_3_class_roc_run_dir

################
# Run analysis at the variant level
################
all_variants_variant_level_input_file=$unsupervised_learning_input_dir"old_all_availibile_samples_variant_level_features_filter_partially_observed_expression.txt"
output_stem=$watershed_run_dir"old_fully_observed_te_ase_splicing_outliers_gene_pvalue_"$gene_thresh"_outlier_fraction_"$pvalue_fraction"_pseudocount_"$pseudocount"_"$inference_method"_inference_apply_to_all_variants"
Rscript watershed_score.R $pvalue_fraction $number_of_dimensions $inference_method $pseudocount $fully_observed_input_file $all_variants_variant_level_input_file $output_stem $binary_pvalue_threshold $watershed_3_class_roc_run_dir


if false; then
Rscript debug_visualize_watershed_3_class_scores.R $unsupervised_learning_input_dir $watershed_run_dir
fi