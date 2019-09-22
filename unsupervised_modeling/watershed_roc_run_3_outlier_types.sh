#!/bin/bash -l

#SBATCH
#SBATCH --time=14:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


unsupervised_learning_input_dir="$1"
watershed_run_dir="$2"
pseudocount="$3"
n2_pair_pvalue_fraction="$4"
binary_pvalue_threshold="$5"



# Parameters!
number_of_dimensions="3"

gene_thresh="0.01"
input_stem="fully_observed_merged_outliers_"$gene_thresh"_genes_intersection_between_te_ase_splicing"
input_file=$unsupervised_learning_input_dir$input_stem"_features_filter_no_tissue_anno_N2_pairs.txt"
output_stem=$watershed_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_"$gene_thresh"_n2_pair_outlier_fraction_"$n2_pair_pvalue_fraction"_binary_pvalue_threshold_"$binary_pvalue_threshold"_pseudocount_"$pseudocount


echo $pseudocount
echo $n2_pair_pvalue_fraction
echo $gene_thresh
echo $binary_pvalue_threshold

Rscript watershed_roc_3_outlier_types.R $input_file $output_stem $number_of_dimensions $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold

if false; then
# Parameters!
number_of_dimensions="3"

gene_thresh="0.01"
input_stem="fully_observed_merged_outliers_"$gene_thresh"_genes_intersection_between_te_ase_splicing"
input_file=$unsupervised_learning_input_dir$input_stem"_features_filter_no_tissue_anno_N2_pairs.txt"
output_stem=$watershed_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_"$gene_thresh"_n2_pair_outlier_fraction_"$n2_pair_pvalue_fraction"_binary_pvalue_threshold_"$binary_pvalue_threshold"_pseudocount_"$pseudocount


echo $pseudocount
echo $n2_pair_pvalue_fraction
echo $gene_thresh
echo $binary_pvalue_threshold

Rscript watershed_roc_3_outlier_types.R $input_file $output_stem $number_of_dimensions $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold



# Parameters!
number_of_dimensions="3"

gene_thresh="0.01"
input_stem="old_fully_observed_merged_outliers_"$gene_thresh"_genes_intersection_between_te_ase_splicing"
input_file=$unsupervised_learning_input_dir$input_stem"_features_filter_no_tissue_anno_N2_pairs.txt"
output_stem=$watershed_run_dir"old_fully_observed_te_ase_splicing_outliers_gene_pvalue_"$gene_thresh"_n2_pair_outlier_fraction_"$n2_pair_pvalue_fraction"_binary_pvalue_threshold_"$binary_pvalue_threshold"_pseudocount_"$pseudocount


echo $pseudocount
echo $n2_pair_pvalue_fraction
echo $gene_thresh
echo $binary_pvalue_threshold

Rscript watershed_roc_3_outlier_types.R $input_file $output_stem $number_of_dimensions $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold
fi
if false; then
# Parameters!
number_of_dimensions="3"

gene_thresh="0.01"
input_stem="rep1_fully_observed_merged_outliers_"$gene_thresh"_genes_intersection_between_te_ase_splicing"
input_file=$unsupervised_learning_input_dir$input_stem"_features_filter_no_tissue_anno_N2_pairs.txt"
output_stem=$watershed_run_dir"rep1_fully_observed_te_ase_splicing_outliers_gene_pvalue_"$gene_thresh"_n2_pair_outlier_fraction_"$n2_pair_pvalue_fraction"_binary_pvalue_threshold_"$binary_pvalue_threshold"_pseudocount_"$pseudocount


echo $pseudocount
echo $n2_pair_pvalue_fraction
echo $gene_thresh
echo $binary_pvalue_threshold

Rscript watershed_roc_3_outlier_types.R $input_file $output_stem $number_of_dimensions $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold
fi
