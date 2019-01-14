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
# A donor is called a "global outlier" if it is a median(pvalue) outlier in at least $num_global_outlier_clusters clusters
num_global_outlier_clusters="500"
# Anything less than p=$pvalue_threshold is considered an outlier
pvalue_threshold=".0456"

if false; then
python identify_and_remove_global_outliers.py $tissue_names_file $splicing_outlier_dir $covariate_method $min_number_of_expressed_tissues $num_global_outlier_clusters $pvalue_threshold $european_ancestry_individual_list
fi


#################
# Create file containing cross-tissue outliers (ie median(pvalue) across observed tissues)
# To be eligible to be a "cross tissue outlier", (individual, cluster) pair must be expressed in at least $min_number_of_expressed_tissues
min_number_of_expressed_tissues="5"
# A donor is called a "global outlier" if it is a median(pvalue) outlier in at least $num_global_outlier_clusters clusters
num_global_outlier_clusters="500"
if false; then
python identify_cross_tissue_outliers.py $tissue_names_file $splicing_outlier_dir $covariate_method $min_number_of_expressed_tissues $european_ancestry_individual_list $num_global_outlier_clusters
fi


#################
# Visualize distribution of outier calls
if false; then
Rscript visualize_outlier_calls.R $tissue_names_file $covariate_method $splicing_outlier_dir $splicing_outlier_visualization_dir
fi
