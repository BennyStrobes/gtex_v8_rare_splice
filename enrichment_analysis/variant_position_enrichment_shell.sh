#!/bin/bash -l

#SBATCH
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --mem=25GB
#SBATCH --nodes=1


rare_variant_dir="$1"
variant_position_enrichment_dir="$2"
visualize_variant_position_enrichment_dir="$3"
splicing_outlier_dir="$4"
splicing_outlier_suffix="$5"
european_ancestry_individual_list="$6"
gencode_gene_annotation_file="$7"
cluster_info_file="$8"
exon_file="$9"
jxn_usage_nearby_altered_ss_enrichment_dir="${10}"
tissue_names_file="${11}"
filtered_cluster_dir="${12}"


########################
# Compare distances between variants and splice sites for outliers vs non-outliers
########################
# Range of pvalues
pvalue_thresholds=( ".000001" ".00001" ".0001" ".001")
distance="1000"

# Loop through pvalue thresholds
if false; then
for pvalue_threshold in "${pvalue_thresholds[@]}"; do
	echo $distance"_"$pvalue_threshold
	python variant_position_enrichment_quantification.py $rare_variant_dir $variant_position_enrichment_dir $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $gencode_gene_annotation_file $cluster_info_file $pvalue_threshold $distance $exon_file
done
fi



########################
# Compare jxn usage nearby A: altered splice sites and B: altered PPT regions in outliers and non-outliers
########################
pvalue_outlier_threshold=".00001"
pvalue_inlier_threshold=".5"
variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_100.txt"
if false; then
python quantify_junction_usage_nearby_altered_splice_site.py $pvalue_outlier_threshold $pvalue_inlier_threshold $variant_bed_file $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $gencode_gene_annotation_file $cluster_info_file $exon_file $jxn_usage_nearby_altered_ss_enrichment_dir $tissue_names_file $filtered_cluster_dir
fi




########################
# Visualize results
########################
Rscript visualize_variant_position_enrichment.R $variant_position_enrichment_dir $jxn_usage_nearby_altered_ss_enrichment_dir $visualize_variant_position_enrichment_dir









