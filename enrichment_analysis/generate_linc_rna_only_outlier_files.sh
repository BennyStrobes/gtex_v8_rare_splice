#!/bin/bash -l

#SBATCH
#SBATCH --time=9:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1

tissue_names_file="$1"
gtex_v8_gene_list="$2"
linc_rna_gene_list="$3"
gtex_v8_outlier_calls_dir="$4"
cluster_info_file="$5"
outlier_calls_dir="$6"


python generate_linc_rna_only_outlier_files.py $tissue_names_file $gtex_v8_gene_list $linc_rna_gene_list $gtex_v8_outlier_calls_dir $cluster_info_file $outlier_calls_dir