#!/bin/bash -l

#SBATCH
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


genomic_annotation_file="$1"
total_expression_outlier_file="$2"
ase_outlier_file="$3"
splicing_outlier_file="$4"
unsupervised_learning_input_dir="$5"


num_outlier_samples="2000"
python make_river_input_files_for_comparison_between_three_outlier_methods.py $genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $num_outlier_samples