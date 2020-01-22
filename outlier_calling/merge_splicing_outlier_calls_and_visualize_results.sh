#!/bin/bash -l

#SBATCH
#SBATCH --time=6:00:00
#SBATCH --mem=25GB
#SBATCH --partition=shared
#SBATCH --nodes=1


tissue_names_file="$1"
covariate_method="$2"
total_jobs="$3"
splicing_outlier_dir="$4"
splicing_outlier_visualization_dir="$5"
european_ancestry_individual_list="$6"
filtered_cluster_dir="$7"



# File containing mapping from junction to cluster to gene
cluster_info_file=$filtered_cluster_dir"cluster_info.txt"




#################
# For each tissue merge outlier calls
if false; then
while read tissue_type; do
	echo $tissue_type
	output_root=$splicing_outlier_dir$tissue_type"_covariate_method_"$covariate_method"_"
	python merge_splicing_outlier_calls.py $output_root $total_jobs
done<$tissue_names_file
fi

#################
# Merge splicing outlier calls in Muscle while varying hyperparameters (num reads in emperical distribution)
tissue_type="Muscle_Skeletal"
# Number of reads to simulate for emperical distribution
num_read_array=("20000" "30000" "10000" "5000" "1000" "500" "100")
num_read_array=("100000")
if false; then
for num_reads in "${num_read_array[@]}"; do
	echo $num_reads
	output_root=$splicing_outlier_dir$tissue_type"_compare_num_read_hyperparam_"$num_reads"_covariate_method_"$covariate_method"_"
	python merge_splicing_outlier_calls.py $output_root $total_jobs
done
fi 

#################
# Merge splicing outlier calls in Muscle while varying hyperparameters (num reads in emperical distribution)
tissue_type="Muscle_Skeletal"
# Number of reads to simulate for emperical distribution
num_read_array=("20000" "10000" "1000" "100")
if false; then
for num_reads in "${num_read_array[@]}"; do
	echo $num_reads
	output_root=$splicing_outlier_dir$tissue_type"_compare_num_read_hyperparam_"$num_reads"_standard_pseudocount_covariate_method_"$covariate_method"_"
	python merge_splicing_outlier_calls.py $output_root $total_jobs
done
fi

#################
# Merge splicing outlier calls in Muscle while varying hyperparameters (# Type of Dirichlet multinomial model to use (either 'standard' or 'no_prior')
tissue_type="Muscle_Skeletal"
# Type of Dirichlet multinomial model to use (either 'standard' or 'no_prior')
model_version_array=("no_prior_multiple_initializations" "standard")
if false; then
for model_version in "${model_version_array[@]}"; do
	echo $model_version
	output_root=$splicing_outlier_dir$tissue_type"_compare_model_version_hyperparam_"$model_version"_covariate_method_"$covariate_method"_"
	python merge_splicing_outlier_calls.py $output_root $total_jobs
done
fi

##################
# Compute fraction of reads coming from each junction
tissue_type="Muscle_Skeletal"
output_file=$splicing_outlier_dir$tissue_type"_max_fraction_of_reads_from_a_junction.txt"
splicing_outlier_file=$splicing_outlier_dir$tissue_type"_covariate_method_none_merged_emperical_pvalue.txt"
tissue_specific_junction_file=$filtered_cluster_dir$tissue_type"_filtered_jxns_cross_tissue_clusters_gene_mapped.txt"
if false; then
python compute_fraction_of_reads_coming_from_each_junction.py $splicing_outlier_file $tissue_specific_junction_file $output_file
fi
####################
# Get gene level pvalues (accounting for the number of clusters we are taking the minimum over)
if false; then
while read tissue_type; do
	echo $tissue_type
	output_root=$splicing_outlier_dir$tissue_type"_covariate_method_"$covariate_method"_"
	python convert_outlier_calls_from_cluster_to_gene_level.py $output_root $cluster_info_file
done<$tissue_names_file
fi



####################
# Get gene level pvalues!!!
# change global outlier identification based on gene level file!
############################

#################
# Identify and remove "global outliers" (ie. donors that are median(pvalue) outliers in a large number of clusters)
# To be eligible to be a "global outlier", (individual, cluster) pair must be expressed in at least $min_number_of_expressed_tissues
min_number_of_expressed_tissues="5"
# Anything less than p=$pvalue_threshold is considered an outlier
pvalue_threshold=".01"
if false; then
python identify_and_remove_global_outliers.py $tissue_names_file $splicing_outlier_dir $covariate_method $min_number_of_expressed_tissues $pvalue_threshold $european_ancestry_individual_list
fi


#################
# Create file containing cross-tissue outliers (ie median(pvalue) across observed tissues)
# To be eligible to be a "cross tissue outlier", (individual, cluster) pair must be expressed in at least $min_number_of_expressed_tissues
if false; then
min_number_of_expressed_tissues="5"
python identify_cross_tissue_outliers.py $tissue_names_file $splicing_outlier_dir $covariate_method $min_number_of_expressed_tissues $european_ancestry_individual_list "emperical_pvalue"
# Do the same for gene level outlier calls
min_number_of_expressed_tissues="5"
python identify_cross_tissue_outliers.py $tissue_names_file $splicing_outlier_dir $covariate_method $min_number_of_expressed_tissues $european_ancestry_individual_list "emperical_pvalue_gene_level"
fi


#################
# Visualize distribution of outier calls
Rscript visualize_outlier_calls.R $tissue_names_file $covariate_method $splicing_outlier_dir $splicing_outlier_visualization_dir

