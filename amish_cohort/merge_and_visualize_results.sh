#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


splicing_outlier_dir="$1"
expression_outlier_dir="$2"
processed_data_dir="$3"
merged_results_dir="$4"
visualize_results_dir="$5"

maf_file=$processed_data_dir"snv.frq"

############################################
########## MERGE Splicing outlier data set
############################################
splicing_outlier_file=$splicing_outlier_dir"amish_cohort_covariate_method_none_merged_no_global_outliers_emperical_pvalue_gene_level.txt"
variant_bed_file=$processed_data_dir"variant_bed_splicing.txt"
merged_data_set_file=$merged_results_dir"merged_data_set_splicing.txt"
merged_compressed_data_set_file=$merged_results_dir"merged_compressed_data_set_splicing.txt"

python merge_data_sets.py $splicing_outlier_file $variant_bed_file $maf_file $merged_data_set_file $merged_compressed_data_set_file



############################################
########## MERGE ASE outlier data set
############################################
processed_ase_outlier_calls_file=$processed_data_dir"ase_outlier_calls.txt"

variant_bed_file=$processed_data_dir"variant_bed_ase.txt"
merged_data_set_file=$merged_results_dir"merged_data_set_ase.txt"
merged_compressed_data_set_file=$merged_results_dir"merged_compressed_data_set_ase.txt"
python merge_data_sets.py $processed_ase_outlier_calls_file $variant_bed_file $maf_file $merged_data_set_file $merged_compressed_data_set_file


############################################
########## MERGE expression outlier data set
############################################
expression_outlier_file=$expression_outlier_dir"amish_pvalue_expression_outlier_calls_no_global.txt"
variant_bed_file=$processed_data_dir"variant_bed_expression.txt"
merged_data_set_file=$merged_results_dir"merged_data_set_expression.txt"
merged_compressed_data_set_file=$merged_results_dir"merged_compressed_data_set_expression.txt"
python merge_data_sets.py $expression_outlier_file $variant_bed_file $maf_file $merged_data_set_file $merged_compressed_data_set_file


expression_outlier_file=$expression_outlier_dir"amish_pvalue_expression_outlier_calls_no_global.txt"
variant_bed_file=$processed_data_dir"variant_bed_rare_expression.txt"
merged_data_set_file=$merged_results_dir"merged_data_set_rare_expression.txt"
merged_compressed_data_set_file=$merged_results_dir"merged_compressed_data_set_rare_expression.txt"
python merge_data_sets.py $expression_outlier_file $variant_bed_file $maf_file $merged_data_set_file $merged_compressed_data_set_file


############################################
########## Visualize results
############################################


splicing_merged_results_file=$merged_results_dir"merged_compressed_data_set_splicing.txt"
ase_merged_results_file=$merged_results_dir"merged_compressed_data_set_ase.txt"
expression_merged_results_file=$merged_results_dir"merged_compressed_data_set_expression.txt"
rare_expression_merged_results_file=$merged_results_dir"merged_compressed_data_set_rare_expression.txt"
Rscript visualize_results.R $splicing_merged_results_file $ase_merged_results_file $expression_merged_results_file $rare_expression_merged_results_file $visualize_results_dir
