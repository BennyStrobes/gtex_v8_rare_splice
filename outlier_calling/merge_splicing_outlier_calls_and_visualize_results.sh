#!/bin/bash -l

#SBATCH
#SBATCH --time=6:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


tissue_names_file="$1"
covariate_method="$2"
total_jobs="$3"
splicing_outlier_dir="$4"
splicing_outlier_visualization_dir="$5"
european_ancestry_individual_list="$6"






#################
# For each tissue merge outlier calls
if false; then
while read tissue_type; do
	echo $tissue_type
	output_root=$splicing_outlier_dir$tissue_type"_covariate_method_"$covariate_method"_"
	python merge_splicing_outlier_calls.py $output_root $total_jobs

done<$tissue_names_file
fi
	tissue_type="Muscle_Skeletal"
	output_root=$splicing_outlier_dir$tissue_type"_covariate_method_"$covariate_method"_"
	python merge_splicing_outlier_calls.py $output_root $total_jobs


#################
# Code identifying "global outliers" as well as multi-tissue outliers. Also some code limiting to European Ancestry individuals