#!/bin/bash -l

#SBATCH
#SBATCH --time=3:00:00
#SBATCH --mem=10GB
#SBATCH --partition=shared
#SBATCH --nodes=1


variant_bed_file="$1"
heuristic_outlier_dir="$2"
heuristic_outlier_suffix="$3"
tissue_names_file="$4"
rare_variant_dir="$5"



# Range of Distances
distances=( "2" "4" "6" "8" "10" "20" "100" "1000")
distances=( "8" )

# Loop through distances
for distance in "${distances[@]}"
do
	echo $distance
	variant_junction_bed_file=$rare_variant_dir"variant_heuristic_junction_bed_"$distance".txt"
	python map_variants_to_junctions.py $variant_bed_file $variant_junction_bed_file $heuristic_outlier_dir $heuristic_outlier_suffix $tissue_names_file $distance
done




