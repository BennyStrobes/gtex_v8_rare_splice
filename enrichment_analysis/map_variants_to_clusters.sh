#!/bin/bash -l

#SBATCH
#SBATCH --time=8:00:00
#SBATCH --mem=10GB
#SBATCH --partition=shared
#SBATCH --nodes=1


variant_bed_file="$1"
cluster_info_file="$2"
rare_variant_dir="$3"




# Range of Distances
distances=( "4" "6" "8" "10" "20" "100" "1000")


# Loop through distances
if false; then
for distance in "${distances[@]}"
do
	echo $distance
	variant_cluster_bed_file=$rare_variant_dir"variant_cluster_bed_"$distance".txt"
	variant_cluster_only_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance".txt"
	python map_variants_to_clusters.py $variant_bed_file $variant_cluster_bed_file $variant_cluster_only_bed_file $cluster_info_file $distance
done
fi






distances=( "100" )

for distance in "${distances[@]}"
do
	variant_cluster_bed_file=$rare_variant_dir"variant_cluster_intron_body_"$distance"_bed.txt"
	variant_cluster_no_consensus_bed_file=$rare_variant_dir"variant_cluster_intron_body_"$distance"_no_consensus_bed.txt"
	python map_variants_to_cluster_intron_bodies.py $variant_bed_file $variant_cluster_bed_file $variant_cluster_no_consensus_bed_file $cluster_info_file $distance
done
