#!/bin/bash -l

#SBATCH
#SBATCH --time=9:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


tissue_type="$1"
tissue_specific_junction_file="$2"
output_root="$3"
gencode_gene_annotation_file="$4"
gene_list="$5"


python call_heuristic_outliers.py $tissue_type $tissue_specific_junction_file $output_root $gencode_gene_annotation_file $gene_list