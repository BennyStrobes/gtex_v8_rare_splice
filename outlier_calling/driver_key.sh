######## Ben Strober
######## 11/27/18




#######################################
# Input Data
#######################################

# Ordered list of GTEx v8 tissue names
tissue_names_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtex_v8_tissues.txt"

# Directory containing one file for each gtex tissue. Each file has leafcutter cluster results for that tissue
# Each file is of the format $tissue_name"_perind.counts.gz"
# Files were created by Francois Aguet with Leafcutter (https://github.com/broadinstitute/gtex-v8/tree/master/pipelines/leafcutter)
# Leafcutter settings min_clu_ratio=0 and min_clu_reads was 30
leafcutter_cluster_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/leafcutter_clusters/"

# File containing list of genes (ensamble IDs) to be used in this analysis
# File provided by Nicoloe Ferraro
gene_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gencode.v26.GRCh38.genes_genetypes_autosomal_PCandlinc_only.txt"

# Gencode Gene annotation file
# Downloaded from scp.nygenome.org:/data/delivery/gtex-rare/data/gencode.v26.GRCh38.genes.gtf on 1/10/19
# Xin put the file there. He said he got in from the GTEx Exchange
gencode_gene_annotation_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gencode.v26.GRCh38.genes.gtf"

# List of individuals used in each tissue generated by Nicole Ferraro
# Downloaded from scp.nygenome.org:/data/delivery/gtex-rare/data/gtexV8_inds_by_tissue.txt on 11/30/18
individual_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtexV8_inds_by_tissue.txt"

# List of european ancestry individuals (extracted from Nicole Ferraro's variant list: "all_rare_variants_SNPs_10kb_genebody.txt")
european_ancestry_individual_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/european_ancestry_individuals.txt"




#############################################################
#Parameters
#############################################################
#We only consider junctions that have at least one sample with greater than or equal to $min_reads
min_reads="15"
# We only consider clusters where $min_fraction_samples_per_cluster fraction of samples have >= $min_reads_per_sample_in_cluster reads summed across all jxns in that cluster
min_reads_per_sample_in_cluster="3"
min_fraction_expressed_samples_per_cluster=".9"
# Max number of junctions per cluster
max_number_of_junctions_per_cluster="20"






#############################################################
#Used Directories (directories need to be created and empty before starting)
#############################################################

# Root directories that all other directories will be placed in 
output_root="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/"

# Directory containing filtered leafcutter clusters
filtered_cluster_dir=$output_root"clusters_filtered/"

# Directory containing visualizations describing clusters/junctions data
cluster_visualization_dir=$output_root"visualize_clusters_filtered/"

# Directory containing splicing outlier calls
splicing_outlier_dir=$output_root"splicing_outlier_calls/"

# Directory containing visualizations of outlier calls
splicing_outlier_visualization_dir=$output_root"visualize_splicing_outlier_calls/"

# Directory containing outlier calls from heuristic approach
heuristic_outlier_dir=$output_root"heuristic_outlier_calls/"






#############################################################
# Scripts: Run the following X parts in Series
#############################################################


#################
# Part 1: Filter Clusters
##	This script pre-processes the junction/cluster data. It works directly on the leafcutter output ($tissue_name"_perind.counts.gz")
##	It first filters junctions if:
#####1. junction contains no individuals with >= $min_reads mapped to the junction
#####2. Junction is on a non-autosomal chromosome
##After these filters, it will re-make leafcutter clusters (after junctions are removed)
# Lastly, only include clusters that:
##### 1. Contain more than 1 junction
##### 2. We only consider clusters where $min_fraction_samples_per_cluster fraction of samples have >= $min_reads_per_sample_in_cluster reads summed across all jxns in that cluster
if false; then
while read tissue_name; do
	# Input file
	raw_leafcutter_cluster_file=$leafcutter_cluster_dir$tissue_name"_perind.counts.gz"
	# Output file
	filtered_leafcutter_cluster_file=$filtered_cluster_dir$tissue_name"_filtered_jxns.txt"
	# Submit batch job
	sbatch filter_clusters.sh $raw_leafcutter_cluster_file $filtered_leafcutter_cluster_file $min_reads $min_reads_per_sample_in_cluster $min_fraction_expressed_samples_per_cluster $individual_list $tissue_name
done<$tissue_names_file
fi




#################
# Part 2: More cluster/junction filtering
######## 1. Generate clusters that are consistent across tissues (ie Cluster 1 in tissue 1 corresponds to the same set of junctions in all other tissues)
######## 2. Map Clusters to genes
######## 3. Visualize clusters
if false; then
sbatch generate_cross_tissue_clusters_and_map_to_genes.sh $tissue_names_file $filtered_cluster_dir $gencode_gene_annotation_file $cluster_visualization_dir $gene_list
fi


#################
# Part 3: Call splicing outliers in each tissue
# How many nodes to run in parallel
total_jobs="10"
# Whether to include covariates in GLM
covariate_method="none"
if false; then
while read tissue_type; do
	for job_number in $(seq 0 `expr $total_jobs - "1"`); do
		tissue_specific_junction_file=$filtered_cluster_dir$tissue_type"_filtered_jxns_cross_tissue_clusters_gene_mapped.txt"
		output_root=$splicing_outlier_dir$tissue_type"_covariate_method_"$covariate_method"_"$job_number"_"$total_jobs
		sbatch call_splicing_outliers.sh $tissue_type $tissue_specific_junction_file $covariate_method $max_number_of_junctions_per_cluster $output_root $job_number $total_jobs
	done
done<$tissue_names_file
fi






#################
# Part 4: Merge outlier calls (across parallelization runs)
# get gene level pvalues (accounting for the number of clusters we are taking the minimum over)
# visualize outlier calls (in each tissue seperately)
if false; then
sh merge_splicing_outlier_calls_and_visualize_results.sh $tissue_names_file $covariate_method $total_jobs $splicing_outlier_dir $splicing_outlier_visualization_dir $european_ancestry_individual_list $filtered_cluster_dir
fi


#################
# Part 5: Make outlier calls in each tissue using heuristic approach from Cummings et. al
if false; then
while read tissue_type; do
	echo $tissue_type
	tissue_specific_junction_file=$filtered_cluster_dir$tissue_type"_filtered_jxns_cross_tissue_clusters_gene_mapped.txt"
	output_root=$heuristic_outlier_dir$tissue_type"_heuristic_approach_"
	sh call_heuristic_outliers.sh $tissue_type $tissue_specific_junction_file $output_root $gencode_gene_annotation_file $gene_list
done<$tissue_names_file
fi
