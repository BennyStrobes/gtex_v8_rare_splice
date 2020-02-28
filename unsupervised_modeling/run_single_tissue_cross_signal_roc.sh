#!/bin/bash -l

#SBATCH
#SBATCH --time=8:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


input_file="$1"
output_stem="$2"
number_of_dimensions="$3"
pseudocount="$4"
n2_pair_pvalue_fraction="$5"
binary_pvalue_threshold="$6"

Rscript watershed_roc_3_outlier_types.R $input_file $output_stem $number_of_dimensions $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold
