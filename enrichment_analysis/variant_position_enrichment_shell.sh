#!/bin/bash -l

#SBATCH
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --mem=25GB
#SBATCH --nodes=1


rare_variant_dir="$1"
variant_enrichment_dir="$2"
variant_position_enrichment_dir="$3"
visualize_variant_position_enrichment_dir="$4"
splicing_outlier_dir="$5"
splicing_outlier_suffix="$6"
european_ancestry_individual_list="$7"
gencode_gene_annotation_file="$8"
cluster_info_file="$9"
exon_file="${10}"
jxn_usage_nearby_altered_ss_enrichment_dir="${11}"
tissue_names_file="${12}"
filtered_cluster_dir="${13}"
tissue_colors_file="${14}"
heuristic_outlier_dir="${15}"
heuristic_outlier_suffix="${16}"
splice_site_cartoon="${17}"
figure_2_ab_data="${18}"


# Whether to take 'top_outlier' per cluster of 'all' variants per cluster
enrichment_version="all"


##################################
# Tissue-by Tissue variant outlier enrichment
##################################
# Range of Distances
distances=( "2" "4" "6" "8" "10" "100" "1000")
# Range of pvalue thresholds
pvalue_thresholds=( ".000001" ".00001" ".0001" )
# Loop through distances
if false; then
for distance_window in "${distances[@]}"; do
	# Loop through range of pvalue thresholds
	for pvalue_threshold in "${pvalue_thresholds[@]}"; do

		echo "Distance window="$distance_window" pvalue_threshold="$pvalue_threshold
		
		# Run Tissue by tissue enrichment of rare variants within outlier calls (after filtering individuals that were "global outliers")
		variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
		output_root=$variant_enrichment_dir"tbt_variant_outlier_enrichment_pvalue_"$pvalue_threshold"_distance_"$distance_window"_version_"$enrichment_version
		python variant_tbt_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_dir $splicing_outlier_suffix"_merged_emperical_pvalue.txt" $european_ancestry_individual_list $tissue_names_file $enrichment_version $pvalue_threshold

	done
done
fi




##################################
# Cross Tissue variant outlier enrichment
##################################
previous_distance="2"
distances=( "2" "4" "6" "8" "10" "100" "1000")
# Loop through distances
if false; then
for distance_window in "${distances[@]}"; do
	echo "Multi-tissue distance="$distance_window

	# Run cross tissue enrichment using this sized window
	splicing_outlier_file=$splicing_outlier_dir"cross_tissue"$splicing_outlier_suffix"_emperical_pvalue.txt"
	variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
	output_root=$variant_enrichment_dir"cross_tissue_variant_outlier_enrichment_distance_"$distance_window"_version_"$enrichment_version
	python variant_cross_tissue_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_file $european_ancestry_individual_list $enrichment_version

	# Run cross tissue enrichment using this sized window
	splicing_outlier_file=$splicing_outlier_dir"cross_tissue"$splicing_outlier_suffix"_emperical_pvalue.txt"
	variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
	previous_variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$previous_distance".txt"
	output_root=$variant_enrichment_dir"cross_tissue_variant_outlier_enrichment_mutually_exclusive_distance_"$distance_window"_version_"$enrichment_version
	python variant_cross_tissue_mutually_exclusive_distance_enrichment_quantification.py $variant_bed_file $previous_variant_bed_file $output_root $splicing_outlier_file $european_ancestry_individual_list $enrichment_version
	previous_distance=$distance_window
done
fi

##################################
# Heuristic approach variant outlier enrichment
##################################
distance_window="8"
pvalue_threshold=".05"
variant_bed_file=$rare_variant_dir"variant_heuristic_junction_bed_"$distance_window".txt"
output_root=$variant_enrichment_dir"tbt_variant_heuristic_outlier_enrichment_distance_"$distance_window"_version_"$enrichment_version
if false; then
python variant_tbt_enrichment_quantification.py $variant_bed_file $output_root $heuristic_outlier_dir $heuristic_outlier_suffix $european_ancestry_individual_list $tissue_names_file $enrichment_version $pvalue_threshold
fi


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
Rscript visualize_variant_position_enrichment.R $variant_position_enrichment_dir $jxn_usage_nearby_altered_ss_enrichment_dir $variant_enrichment_dir $tissue_names_file $tissue_colors_file $visualize_variant_position_enrichment_dir $splice_site_cartoon $figure_2_ab_data









