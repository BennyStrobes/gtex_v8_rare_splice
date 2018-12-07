#!/bin/bash -l

#SBATCH
#SBATCH --time=2:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


rare_variant_dir="$1"
variant_enrichment_dir="$2"
splicing_outlier_dir="$3"
splicing_outlier_suffix="$4"
european_ancestry_individual_list="$5"
tissue_names_file="$6"


pvalue_threshold=".000001"
distance_window="8"
enrichment_version="all"
variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
output_root=$variant_enrichment_dir"multi_tissue_variant_outlier_enrichment_pvalue_"$pvalue_threshold"_distance_"$distance_window"_version_"$enrichment_version
python variant_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $tissue_names_file $enrichment_version $pvalue_threshold


pvalue_threshold=".00001"
distance_window="8"
enrichment_version="all"
variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
output_root=$variant_enrichment_dir"multi_tissue_variant_outlier_enrichment_pvalue_"$pvalue_threshold"_distance_"$distance_window"_version_"$enrichment_version
python variant_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $tissue_names_file $enrichment_version $pvalue_threshold

pvalue_threshold=".0001"
distance_window="8"
enrichment_version="all"
variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
output_root=$variant_enrichment_dir"multi_tissue_variant_outlier_enrichment_pvalue_"$pvalue_threshold"_distance_"$distance_window"_version_"$enrichment_version
python variant_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $tissue_names_file $enrichment_version $pvalue_threshold


pvalue_threshold=".001"
distance_window="8"
enrichment_version="all"
variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
output_root=$variant_enrichment_dir"multi_tissue_variant_outlier_enrichment_pvalue_"$pvalue_threshold"_distance_"$distance_window"_version_"$enrichment_version
python variant_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $tissue_names_file $enrichment_version $pvalue_threshold


