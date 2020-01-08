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



##############################
# Make feature count file have correct sample names (some manual processing here)
##############################
processed_feature_count_file=$expression_outlier_dir"processed_feature_counts.txt"
if false; then
python process_feature_count_file_to_have_correct_sample_names.py $raw_count_file $sample_mapping_file $processed_feature_count_file $splicing_individual_file
fi

##############################
# Convert raw counts to tpm
##############################
count_file=$expression_outlier_dir"counts.txt"
tpm_file=$expression_outlier_dir"tpm.txt"
log_tpm_file=$expression_outlier_dir"log_tpm.txt"
if false; then
Rscript raw_counts_to_tpm.R $processed_feature_count_file $count_file $tpm_file $log_tpm_file
fi


##############################
# filter to expressed genes
##############################
log_tpm_filtered_file=$expression_outlier_dir"log_tpm_filtered_genes.txt"
if false; then
python filter_total_expression_genes.py $count_file $tpm_file $log_tpm_file $gtex_gene_list $log_tpm_filtered_file
fi

##############################
# Standardize total expresssion
##############################
log_tpm_filtered_standardized_file=$expression_outlier_dir"log_tpm_filtered_genes_standardized.txt"
if false; then
python standardize_total_expression.py $log_tpm_filtered_file $log_tpm_filtered_standardized_file
fi


##############################
# Regress out covariates
##############################
log_tpm_filtered_standardized_cov_corrected_file=$expression_outlier_dir"log_tpm_filtered_genes_standardized_covariate_corrected.txt"
if false; then
python regress_out_covariates_for_total_expression_outliers.py $log_tpm_filtered_standardized_file $log_tpm_filtered_standardized_cov_corrected_file
fi

##############################
# Standardize total expresssion again
##############################
expression_outlier_file=$expression_outlier_dir"amish_expression_outlier_calls.txt"
if false; then
python standardize_total_expression.py $log_tpm_filtered_standardized_cov_corrected_file $expression_outlier_file
fi


##############################
# zscore to pvalue
##############################
pvalue_expression_outlier_file=$expression_outlier_dir"amish_pvalue_expression_outlier_calls.txt"
if false; then
python zscore_to_pvalue_outliers.py $expression_outlier_file $pvalue_expression_outlier_file
fi

##############################
# Remove global outliers
##############################
pvalue_expression_outlier_no_global_file=$expression_outlier_dir"amish_pvalue_expression_outlier_calls_no_global.txt"
if false; then
python remove_global_outliers.py $pvalue_expression_outlier_file $pvalue_expression_outlier_no_global_file
fi
