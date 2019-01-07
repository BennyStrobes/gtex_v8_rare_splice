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
exon_file="$9"
sqtl_dir="${10}"



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
		python variant_position_enrichment_quantification.py $rare_variant_dir $variant_position_enrichment_dir $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $gencode_gene_annotation_file $cluster_info_file $pvalue_threshold $distance $exon_file
	done
done
fi

########################
# Compare distances between variants and splice sites for sQTLs vs background
########################
tissue_type="Muscle_Skeletal"
sqtl_file=$sqtl_dir$tissue_type".v8.sqtl_allpairs.txt"
distance_window="10"
pvalue_thresh=".00000000001"
if false; then
python variant_position_enrichment_quantification_in_sqtls.py $sqtl_file $gencode_gene_annotation_file $variant_position_enrichment_dir $distance_window $pvalue_thresh
fi



########################
# Visualize results
########################
Rscript visualize_variant_position_enrichment.R $variant_position_enrichment_dir $visualize_variant_position_enrichment_dir




