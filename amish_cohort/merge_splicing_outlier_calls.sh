#!/bin/bash -l

#SBATCH
#SBATCH --time=9:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


splicing_outlier_dir="$1"
covariate_method="$2"
cluster_info_file="$3"
total_jobs="$4"




# Merge splicing outlier calls across parrallelized runs
output_root=$splicing_outlier_dir"amish_cohort_covariate_method_"$covariate_method"_"
python merge_splicing_outlier_calls.py $output_root $total_jobs



####################
# Get gene level pvalues (accounting for the number of clusters we are taking the minimum over)
output_root=$splicing_outlier_dir"amish_cohort_covariate_method_"$covariate_method"_"
python convert_outlier_calls_from_cluster_to_gene_level.py $output_root $cluster_info_file




####################
# Remove global outliers
gene_level_outlier_file=$splicing_outlier_dir"amish_cohort_covariate_method_"$covariate_method"_merged_emperical_pvalue_gene_level.txt"
gene_level_outlier_file_no_global=$splicing_outlier_dir"amish_cohort_covariate_method_"$covariate_method"_merged_no_global_outliers_emperical_pvalue_gene_level.txt"
python remove_global_outliers.py $gene_level_outlier_file $gene_level_outlier_file_no_global

