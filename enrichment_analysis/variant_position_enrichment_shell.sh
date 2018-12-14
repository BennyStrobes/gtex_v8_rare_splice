#!/bin/bash -l

#SBATCH
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1


rare_variant_dir="$1"
variant_position_enrichment_dir="$2"
visualize_variant_position_enrichment_dir="$3"
splicing_outlier_dir="$4"
splicing_outlier_suffix="$5"
european_ancestry_individual_list="$6"
gencode_gene_annotation_file="$7"
cluster_info_file="$8"
variant_position_enrichment_debug_dir="$9"


########################
# Compare distances between variants and splice sites for outliers vs non-outliers
########################

# Range of Distances
distances=( "10" "20" "100" )
# Range of pvalues
pvalue_thresholds=( ".000001" ".00001" ".0001" ".001")

# Loop through distances
if false; then
for distance in "${distances[@]}"; do
# Loop through pvalue thresholds
	for pvalue_threshold in "${pvalue_thresholds[@]}"; do
		echo $distance"_"$pvalue_threshold
		python variant_position_enrichment_quantification.py $rare_variant_dir $variant_position_enrichment_dir $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $gencode_gene_annotation_file $cluster_info_file $pvalue_threshold $distance
	done
done
fi

########################
# Visualize results
########################
if false; then
Rscript visualize_variant_position_enrichment.R $variant_position_enrichment_dir $visualize_variant_position_enrichment_dir
fi



# Debugging directory containing visualizations of variant position analysis (ie distance between RV and splice sites)
# With the goal of understanding difference between defining exons based on gencode vs calling exons via exon-exon junction data
pvalue_threshold=".0001"
distance="20"
python debug_variant_position_enrichment_quantification.py $rare_variant_dir $variant_position_enrichment_debug_dir $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $gencode_gene_annotation_file $cluster_info_file $pvalue_threshold $distance