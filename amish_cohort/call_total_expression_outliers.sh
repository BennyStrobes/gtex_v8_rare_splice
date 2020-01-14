#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1

raw_count_file="$1"
expression_outlier_dir="$2"
sample_mapping_file="$3"
splicing_individual_file="$4"
gtex_gene_list="$5"
gtex_lymphocyte_count_file="$6"
gtex_lymphocyte_standardized_expression_file="$7"
gtex_gene_length_file="$8"



##############################
# Filter gtex count file to only include samples used in rare analysis
##############################
filtered_gtex_lymphocyte_count_file=$expression_outlier_dir"gtex_lymphocyte_counts.txt"
python filter_samples_in_gtex_count_file.py $gtex_lymphocyte_count_file $gtex_lymphocyte_standardized_expression_file $filtered_gtex_lymphocyte_count_file



##############################
# Make feature count file have correct sample names (some manual processing here)
##############################
processed_feature_count_file=$expression_outlier_dir"processed_feature_counts.txt"
python process_feature_count_file_to_have_correct_sample_names.py $raw_count_file $sample_mapping_file $processed_feature_count_file $splicing_individual_file



##############################
# Merge (concatenate) gtex lymphocyte counts and amish counts into one file
##############################
concatenated_processed_feature_count_file=$expression_outlier_dir"concatenated_processed_feature_counts.txt"
python concatenate_gtex_lymphocyte_counts_and_amish_counts.py $filtered_gtex_lymphocyte_count_file $processed_feature_count_file $concatenated_processed_feature_count_file $gtex_gene_length_file


##############################
# Convert raw counts to tpm
##############################
count_file=$expression_outlier_dir"counts.txt"
tpm_file=$expression_outlier_dir"tpm.txt"
log_tpm_file=$expression_outlier_dir"log_tpm.txt"
Rscript raw_counts_to_tpm.R $concatenated_processed_feature_count_file $count_file $tpm_file $log_tpm_file



##############################
# filter to expressed genes
##############################
log_tpm_filtered_file=$expression_outlier_dir"log_tpm_filtered_genes.txt"
python filter_total_expression_genes.py $count_file $tpm_file $log_tpm_file $gtex_gene_list $log_tpm_filtered_file


##############################
# Standardize total expresssion
##############################
log_tpm_filtered_standardized_file=$expression_outlier_dir"log_tpm_filtered_genes_standardized.txt"
python standardize_total_expression.py $log_tpm_filtered_file $log_tpm_filtered_standardized_file



##############################
# Regress out covariates
##############################
log_tpm_filtered_standardized_cov_corrected_file=$expression_outlier_dir"log_tpm_filtered_genes_standardized_covariate_corrected.txt"
python regress_out_covariates_for_total_expression_outliers.py $log_tpm_filtered_standardized_file $log_tpm_filtered_standardized_cov_corrected_file


##############################
# Standardize total expresssion again
##############################
expression_outlier_file=$expression_outlier_dir"amish_expression_outlier_calls.txt"
python standardize_total_expression.py $log_tpm_filtered_standardized_cov_corrected_file $expression_outlier_file



##############################
# zscore to pvalue
##############################
pvalue_expression_outlier_file=$expression_outlier_dir"amish_pvalue_expression_outlier_calls.txt"
python zscore_to_pvalue_outliers.py $expression_outlier_file $pvalue_expression_outlier_file


##############################
# Remove global outliers
##############################
pvalue_expression_outlier_no_global_file=$expression_outlier_dir"amish_pvalue_expression_outlier_calls_no_global.txt"
python remove_global_outliers.py $pvalue_expression_outlier_file $pvalue_expression_outlier_no_global_file

