#!/bin/bash -l

#SBATCH
#SBATCH --time=9:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


jxn_file="$1"
outlier_pvalue_file="$2"
output_dir="$3"



filtered_jxn_file=$output_dir"exon_exon_junction_file.txt"
if false; then
python filter_junction_file_for_github_repo.py $jxn_file $outlier_pvalue_file $filtered_jxn_file
fi
covariate_method="none"
max_number_of_junctions_per_cluster="20"
output_root=$output_dir"practice_muscle"
job_number="0"
total_jobs="1"
sh call_splicing_outliers.sh "Muscle_Skeletal" $filtered_jxn_file $covariate_method $max_number_of_junctions_per_cluster $output_root $job_number $total_jobs
