#!/bin/bash -l

#SBATCH
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1

rare_variant_dir="$1"
splicing_outlier_dir="$2"
splicing_outlier_suffix="$3"
european_ancestry_individual_list="$4"
cluster_info_file="$5"
exon_file="$6"
visualize_cluster_distribution_dir="$7"
tissue_names_file="$8"
filtered_cluster_dir="$9"


############################################
# Could theoretically loop through all tissues. But not necessary now!
#############################################
tissue_name="Muscle_Skeletal"

distance_window="6"
pvalue_threshold=".000001"
enrichment_version="all"
##############
# Part 1: Extract file containing cluster_id\tvariant_position\toutlier_indi\tinlier_individuals
variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
clusters_to_plot_file=$visualize_cluster_distribution_dir$tissue_name"_outliers_with_rv_nearby_ss_"$distance_window"_pvalue_"$pvalue_threshold".txt"
if false; then
python extract_list_of_outliers_with_rv_nearby_ss.py $tissue_name $pvalue_threshold $variant_bed_file $clusters_to_plot_file $splicing_outlier_dir $splicing_outlier_suffix"_merged" $european_ancestry_individual_list $enrichment_version $cluster_info_file $exon_file
fi

##############
# Part 2: For each cluster in $clusters_to_plot_file, save matrix of cluster counts for outlier indi and inlier indi
tissue_specific_junction_file=$filtered_cluster_dir$tissue_name"_filtered_jxns_cross_tissue_clusters_gene_mapped.txt"
python extract_cluster_counts.py $clusters_to_plot_file $tissue_name $tissue_specific_junction_file $visualize_cluster_distribution_dir


##############
# Part 3: Visualize clusters in $clusters_to_plot
# This code is strongly based off of code from Leafcutter (https://github.com/davidaknowles/leafcutter/blob/master/leafcutter/R/junction_plot.R)
# (NEEDS TO BE RUN LOCALLY!)
clusters_to_plot_file="Muscle_Skeletal_outliers_with_rv_nearby_ss_6_pvalue_.000001.txt"
tissue_name="Muscle_Skeletal"
exon_file="gencode_v26_exons.txt"
if false; then
Rscript visualize_cluster_distribution.R $clusters_to_plot_file $tissue_name $exon_file
fi