######## Ben Strober
######## 12/19/19




#######################################
# Input Data (Assumes everything else has previously been run)
#######################################

# File containing mapping from clusterID to to genes and junctions
cluster_info_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/cluster_info.txt"

# Directory containing filtered junction read counts
whole_blood_junction_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/Whole_Blood_filtered_jxns_cross_tissue_clusters_gene_mapped.txt"

# Manually generated file by me containing variants we will test
test_variant_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/select_crispr_splice_variants/input_data/test_variants.txt"



#######################################
# Output Directories
#######################################
power_analysis_data_output_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/select_crispr_splice_variants/power_analysis/"




python run_power_analysis.py $cluster_info_file $whole_blood_junction_file $test_variant_file $power_analysis_data_output_dir