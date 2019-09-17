#!/bin/bash -l

#SBATCH
#SBATCH --time=3:00:00
#SBATCH --mem=10GB
#SBATCH --partition=shared
#SBATCH --nodes=1


variant_bed_file="$1"
cluster_info_file="$2"
rare_variant_dir="$3"
ctcf_interaction_file="$4"
liftover_directory="$5"




if false; then
# Range of Distances
distances=( "2" "4" "6" "8" "10" "20" "100" "1000")


# Loop through distances
for distance in "${distances[@]}"
do
	echo $distance
	variant_cluster_bed_file=$rare_variant_dir"variant_cluster_bed_"$distance".txt"
	variant_cluster_only_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance".txt"
	python map_variants_to_clusters.py $variant_bed_file $variant_cluster_bed_file $variant_cluster_only_bed_file $cluster_info_file $distance
done
fi


# Two parts:
# 1. Map genomic coordinates to hg38
# 2. Map genomic coordinates to hg38 ctcf interactions
python map_variants_to_cluters_based_on_ctcf_interactions.py $variant_bed_file $cluster_info_file $ctcf_interaction_file $liftover_directory $rare_variant_dir



