#!/bin/bash -l

#SBATCH
#SBATCH --time=9:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


test_junction_file="$1"
train_junction_file="$2"
covariate_method="$3"
max_number_of_junctions_per_cluster="$4"
output_root="$5"
num_reads="$6"
model_version="$7"
job_number="$8"
total_jobs="$9"



date
python call_splicing_outliers_train_test.py $test_junction_file $train_junction_file $covariate_method $max_number_of_junctions_per_cluster $output_root $num_reads $model_version $job_number $total_jobs
date