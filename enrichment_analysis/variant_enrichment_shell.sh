#!/bin/bash -l

#SBATCH
#SBATCH --time=14:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


rare_variant_dir="$1"
variant_enrichment_dir="$2"
splicing_outlier_dir="$3"
splicing_outlier_suffix="$4"
splicing_outlier_include_global_outliers_suffix="$5"
european_ancestry_individual_list="$6"
tissue_names_file="$7"
visualize_variant_enrichment_dir="$8"
tissue_colors_file="$9"

# Whether to take 'top_outlier' per cluster of 'all' variants per cluster
enrichment_version="all"


##################################
# Tissue-by Tissue variant outlier enrichment
##################################

# Range of Distances
distances=( "4" "6" "8" "10" "100" "1000")
# Range of pvalue thresholds
pvalue_thresholds=( ".000001" ".00001" ".0001" ".001" )

# Loop through distances
for distance_window in "${distances[@]}"; do
	# Loop through range of pvalue thresholds
	for pvalue_threshold in "${pvalue_thresholds[@]}"; do

		echo "Distance window="$distance_window" pvalue_threshold="$pvalue_threshold
		
		# Run Tissue by tissue enrichment of rare variants within outlier calls (before filtering individuals that were "global outliers")
		variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
		output_root=$variant_enrichment_dir"tbt_variant_outlier_include_global_outlier_enrichment_pvalue_"$pvalue_threshold"_distance_"$distance_window"_version_"$enrichment_version
		python variant_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_dir $splicing_outlier_include_global_outliers_suffix"_merged" $european_ancestry_individual_list $tissue_names_file $enrichment_version $pvalue_threshold

		
		# Run Tissue by tissue enrichment of rare variants within outlier calls (after filtering individuals that were "global outliers")
		variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
		output_root=$variant_enrichment_dir"tbt_variant_outlier_enrichment_pvalue_"$pvalue_threshold"_distance_"$distance_window"_version_"$enrichment_version
		python variant_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_dir $splicing_outlier_suffix"_merged" $european_ancestry_individual_list $tissue_names_file $enrichment_version $pvalue_threshold

	done
done





##################################
# Cross Tissue variant outlier enrichment
##################################

distances=( "4" "6" "8" "10" "100" "1000")
# Loop through distances
for distance_window in "${distances[@]}"; do

	echo "Multi-tissue distance="$distance_window

	# Run cross tissue enrichment using this sized window
	splicing_outlier_file=$splicing_outlier_dir"cross_tissue"$splicing_outlier_suffix"_emperical_pvalue.txt"
	variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
	output_root=$variant_enrichment_dir"cross_tissue_variant_outlier_enrichment_distance_"$distance_window"_version_"$enrichment_version
	python variant_cross_tissue_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_file $european_ancestry_individual_list $enrichment_version
done




##################################
# Visualize enrichments of rare variants within splicing outlier calls
##################################
if false; then
Rscript visualize_enrichment_of_rare_variants_within_outliers.R $visualize_variant_enrichment_dir $variant_enrichment_dir $tissue_names_file $tissue_colors_file
fi
