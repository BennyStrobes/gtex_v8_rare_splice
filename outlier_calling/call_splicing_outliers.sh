#!/bin/bash -l

#SBATCH
#SBATCH --time=6:00:00
#SBATCH --mem=10GB
#SBATCH --partition=shared
#SBATCH --nodes=1


tissue_type="$1"
tissue_specific_junction_file="$2"
covariate_method="$3"
min_reads_per_sample_in_cluster="$4"
max_number_of_junctions_per_cluster="$5"
output_root="$6"
job_number="$7"
total_jobs="$8"




python call_splicing_outliers.py $tissue_type $tissue_specific_junction_file $covariate_method $min_reads_per_sample_in_cluster $max_number_of_junctions_per_cluster $output_root $job_number $total_jobs