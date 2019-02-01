#!/bin/bash -l

#SBATCH
#SBATCH --time=5:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


genomic_annotation_file="$1"
total_expression_outlier_file="$2"
ase_outlier_file="$3"
splicing_outlier_file="$4"
unsupervised_learning_input_dir="$5"
gene_individual_to_variant_mapping_file="$6"

if false; then
pvalue=".001"
python prepare_input_files_for_unsupervised_learning_union_te_splicing.py $genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $pvalue $gene_individual_to_variant_mapping_file
fi
if false; then
pvalue=".01"
python prepare_input_files_for_unsupervised_learning_intersection_te_splicing.py $genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $pvalue $gene_individual_to_variant_mapping_file

pvalue=".05"
python prepare_input_files_for_unsupervised_learning_intersection_te_splicing.py $genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $pvalue $gene_individual_to_variant_mapping_file
fi

pvalue=".04"
python prepare_input_files_for_unsupervised_learning_intersection_te_ase_splicing.py $genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $pvalue $gene_individual_to_variant_mapping_file

pvalue=".05"
python prepare_input_files_for_unsupervised_learning_intersection_te_ase_splicing.py $genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $pvalue $gene_individual_to_variant_mapping_file
