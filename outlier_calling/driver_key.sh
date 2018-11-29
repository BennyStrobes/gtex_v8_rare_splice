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
# Downloaded from https://www.gencodegenes.org/human/release_26.html on 11/26/18
gencode_gene_annotation_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gencode.v26.annotation.gff3.gz"



#############################################################
#Parameters
#############################################################
#We only consider junctions that have at least one sample with greater than or equal to $min_reads
min_reads="15"
# We only consider clusters that have at least $min_samples_per_cluster samples that have >= $min_reads_per_sample_in_cluster reads summed across all jxns in that cluster
min_reads_per_sample_in_cluster="5"
min_samples_per_cluster="50"






#############################################################
#Used Directories (directories need to be created and empty before starting)
#############################################################

# Root directories that all other directories will be placed in 
output_root="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/"

# Directory containing filtered leafcutter clusters
filtered_cluster_dir=$output_root"clusters_filtered/"

# Directory containing visualizations describing clusters/junctions data
cluster_visualization_dir=$output_root"visualize_clusters_filtered/"























#################
# Part 1: Filter Clusters
##	This script does a number of different cleaning/organization related activities. It works directly on the leafcutter output ($tissue_name"_perind.counts.gz")
##	It first filters junctions if:
#####1. junction contains no individuals with >= $min_reads mapped to the junction
#####2. Junction is on a non-autosomal chromosome
##After these filters, it will re-make leafcutter clusters (after junctions are removed)
# Lastly, only include clusters that:
##### 1. Contain more than 1 junction
##### 2. Has at least $min_samples_per_cluster with $min_reads_per_sample_in_cluster summed across all valid junctions
if false; then
while read tissue_name; do
	# Input file
	raw_leafcutter_cluster_file=$leafcutter_cluster_dir$tissue_name"_perind.counts.gz"
	# Output file
	filtered_leafcutter_cluster_file=$filtered_cluster_dir$tissue_name"_filtered_jxns.txt"
	# Submit batch job
	sbatch filter_clusters.sh $raw_leafcutter_cluster_file $filtered_leafcutter_cluster_file $min_reads $min_reads_per_sample_in_cluster $min_samples_per_cluster
done<$tissue_names_file
fi



#################
# Part 2: Generate clusters across tissues and map clusters to genes

if false; then
sh generate_cross_tissue_clusters_and_map_to_genes.sh $tissue_names_file $filtered_cluster_dir $gencode_gene_annotation_file $cluster_visualization_dir
fi


