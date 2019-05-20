#!/bin/bash -l

#SBATCH
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1

rare_variant_dir="$1"
splicing_outlier_dir="$2"
splicing_outlier_suffix="$3"
european_ancestry_individual_list="$4"
gencode_gene_annotation_file="$5"
cluster_info_file="$6"
exon_file="$7"
jxn_usage_nearby_altered_ss_enrichment_dir="$8"
visualize_jxn_usage_nearby_altered_ss_enrichment_dir="$9"
tissue_names_file="${10}"
filtered_cluster_dir="${11}"



pvalue_outlier_threshold=".000001"
pvalue_inlier_threshold=".5"
variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_100.txt"

if false; then
python quantify_junction_usage_nearby_altered_splice_site.py $pvalue_outlier_threshold $pvalue_inlier_threshold $variant_bed_file $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $gencode_gene_annotation_file $cluster_info_file $exon_file $jxn_usage_nearby_altered_ss_enrichment_dir $tissue_names_file $filtered_cluster_dir
fi

Rscript visualize_junction_usage_nearby_altered_splice_site.R $jxn_usage_nearby_altered_ss_enrichment_dir $visualize_jxn_usage_nearby_altered_ss_enrichment_dir
