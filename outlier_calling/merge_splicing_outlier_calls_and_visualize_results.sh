#!/bin/bash -l

#SBATCH
#SBATCH --time=6:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


tissue_names_file="$1"
min_reads_per_sample_in_cluster="$2"
covariate_method="$3"
total_jobs="$4"
splicing_outlier_dir="$5"
splicing_outlier_visualization_dir="$6"
european_ancestry_individual_list="$7"






#################
# For each tissue merge outlier calls
if false; then
while read tissue_type; do
	echo $tissue_type
	output_root=$splicing_outlier_dir$tissue_type"_min_reads_per_sample_"$min_reads_per_sample_in_cluster"_covariate_method_"$covariate_method"_"
	python merge_splicing_outlier_calls.py $output_root $total_jobs

done<$tissue_names_file
fi
	if false; then
	tissue_type="Muscle_Skeletal"
	output_root=$splicing_outlier_dir$tissue_type"_min_reads_per_sample_"$min_reads_per_sample_in_cluster"_covariate_method_"$covariate_method"_"
	python merge_splicing_outlier_calls.py $output_root $total_jobs
	fi


#################
# Code identifying "global outliers" as well as multi-tissue outliers. Also some code limiting to European Ancestry individuals