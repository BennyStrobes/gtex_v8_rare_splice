#!/bin/bash -l

#SBATCH
#SBATCH --time=5:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


genomic_annotation_file="$1"
variant_level_genomic_annotation_file="$2"
total_expression_outlier_file="$3"
ase_outlier_file="$4"
splicing_outlier_file="$5"
unsupervised_learning_input_dir="$6"
gene_individual_to_variant_mapping_file="$7"
splicing_outlier_dir="$8"
ase_outlier_dir="${9}"
te_outlier_dir="${10}"
tissue_names_file="${11}"
ase_old_outlier_file="${12}"

random_seed="3"


pvalue=".01"
if false; then
python prepare_input_files_for_unsupervised_learning_intersection_te_ase_splicing.py $genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $pvalue $gene_individual_to_variant_mapping_file $random_seed
fi


gene_pvalue_thresh_arr=(".05" ".1")
for pvalue in "${gene_pvalue_thresh_arr[@]}"; do
	echo $pvalue
	python prepare_input_files_for_unsupervised_learning_intersection_te_ase_splicing_at_other_gene_pvalue_thresholds.py $genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $pvalue $gene_individual_to_variant_mapping_file $random_seed
done



pvalue=".01"
python prepare_input_files_for_unsupervised_learning_union_te_ase_splicing_comparison.py $genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $pvalue $gene_individual_to_variant_mapping_file $random_seed




if false; then
##############################
# Generate TBT input files
##############################
# For splicing
python prepare_input_files_for_tbt_splicing.py $unsupervised_learning_input_dir $pvalue $splicing_outlier_dir $tissue_names_file
# For ASE
python prepare_input_files_for_tbt_ase.py $unsupervised_learning_input_dir $pvalue $ase_outlier_dir $tissue_names_file
# For TE
python prepare_input_files_for_tbt_te.py $unsupervised_learning_input_dir $pvalue $te_outlier_dir $tissue_names_file
fi

if false; then
##############################
# Generate Input files for all variants
##############################
python prepare_input_files_for_unsupervised_learning_all_variants.py $genomic_annotation_file $variant_level_genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir
fi

































#################################
# Old (no longer used)
#################################
if false; then
# When we were using old ase outlier file
python prepare_input_files_for_unsupervised_learning_intersection_te_ase_splicing.py $genomic_annotation_file $total_expression_outlier_file $ase_old_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir"old_" $pvalue $gene_individual_to_variant_mapping_file "1"
python prepare_input_files_for_unsupervised_learning_all_variants.py $genomic_annotation_file $variant_level_genomic_annotation_file $total_expression_outlier_file $ase_old_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir"old_"
fi
