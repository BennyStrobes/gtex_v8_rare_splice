#!/bin/bash -l

#SBATCH
#SBATCH --time=5:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1



unsupervised_learning_input_dir="$1"
watershed_run_dir="$2"
pseudocount="$3"
n2_pair_pvalue_fraction="$4"
binary_pvalue_threshold="$5"
phi_method="$6"
lambda_init="$7"
lambda_pair_init="$8"
independent_variables="$9"
inference_method="${10}"
outlier_type="${11}"
number_of_dimensions="${12}"


# Parameters!
gene_thresh="0.01"
input_stem=$outlier_type"_tbt_outliers_"$gene_thresh"_genes_intersection_between_te_ase_splicing"
input_file=$unsupervised_learning_input_dir$input_stem"_features_filter_N2_pairs.txt"
output_stem=$watershed_run_dir$outlier_type"_tbt_intersect_te_ase_splicing_out_gene_pvalue_"$gene_thresh"_n2_pair_outlier_fraction_"$n2_pair_pvalue_fraction"_binary_pvalue_threshold_"$binary_pvalue_threshold"_pseudocount_"$pseudocount"_"$phi_method"_"$lambda_init"_"$lambda_pair_init

echo $pseudocount
echo $n2_pair_pvalue_fraction
echo $gene_thresh
echo $binary_pvalue_threshold
echo $phi_method
echo $lambda_init
echo $lambda_pair_init
echo $independent_variables
echo $inference_method
echo $outlier_type


Rscript watershed_roc_tbt.R $input_file $output_stem $number_of_dimensions $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $phi_method $lambda_init $lambda_pair_init $independent_variables $inference_method $outlier_type
