#!/bin/bash -l

#SBATCH
#SBATCH --time=6:00:00
#SBATCH --mem=10GB
#SBATCH --partition=shared
#SBATCH --nodes=1


tissue_type="$1"
tissue_specific_junction_file="$2"
covariate_method="$3"
max_number_of_junctions_per_cluster="$4"
output_root="$5"
job_number="$6"
total_jobs="$7"




python call_splicing_outliers.py $tissue_type $tissue_specific_junction_file $covariate_method $max_number_of_junctions_per_cluster $output_root $job_number $total_jobs