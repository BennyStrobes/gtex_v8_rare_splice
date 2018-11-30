#!/bin/bash -l

#SBATCH
#SBATCH --time=5:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1

raw_leafcutter_cluster_file="$1"
filtered_leafcutter_cluster_file="$2"
min_reads="$3"
min_reads_per_sample_in_cluster="$4"
min_samples_per_cluster="$5"
individual_list="$6"
tissue_name="$7"

python filter_clusters.py $raw_leafcutter_cluster_file $filtered_leafcutter_cluster_file $min_reads $min_reads_per_sample_in_cluster $min_samples_per_cluster $individual_list $tissue_name