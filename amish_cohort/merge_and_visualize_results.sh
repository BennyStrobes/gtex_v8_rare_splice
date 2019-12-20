#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


splicing_outlier_dir="$1"
processed_data_dir="$2"
merged_results_dir="$3"
visualize_results_dir="$4"


splicing_outlier_file=$splicing_outlier_dir"amish_cohort_covariate_method_none_merged_emperical_pvalue_gene_level.txt"
variant_bed_file=$processed_data_dir"variant_bed_file.txt"

merged_data_set_file=$merged_results_dir"merged_data_set.txt"
merged_compressed_data_set_file=$merged_results_dir"merged_compressed_data_set.txt"

if false; then
python merge_data_sets.py $splicing_outlier_file $variant_bed_file $merged_data_set_file $merged_compressed_data_set_file
fi

if false; then
Rscript visualize_results.R $merged_compressed_data_set_file $visualize_results_dir
fi