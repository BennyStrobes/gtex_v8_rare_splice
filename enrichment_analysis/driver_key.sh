######## Ben Strober
######## 11/27/18




#######################################
# Input Data (Assumes "outlier_calling" (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/outlier_calling) has been already run)
#######################################


# Ordered list of GTEx v8 tissue names
tissue_names_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtex_v8_tissues.txt"

# List of individuals used in each tissue generated by Nicole Ferraro
# Downloaded from scp.nygenome.org:/data/delivery/gtex-rare/data/gtexV8_inds_by_tissue.txt on 11/30/18
individual_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtexV8_inds_by_tissue.txt"

# List of european ancestry individuals (extracted from Nicole Ferraro's variant list: "all_rare_variants_SNPs_10kb_genebody.txt")
european_ancestry_individual_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/european_ancestry_individuals.txt"

# List of rare variants across all European ancestry individuals (Created by Nicole Ferraro. Downloaded from stroberb-420@scp.nygenome.org:/data/delivery/gtex-rare/data/all_rare_variants_noWindow_SNPs.txt on 12/5/18)
variant_bed_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/all_rare_variants_noWindow_SNPs.txt"

# File containing all cluster_ids and their corresponding junction positions
# Generated with "outlier_calling" (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/outlier_calling)
cluster_info_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/cluster_info.txt"

# Directory containing splicing outlier calling results
# By running "outlier_calling" (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/outlier_calling)
splicing_outlier_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/splicing_outlier_calls/"

# Suffix for tissue-specific file_names in $splicing_outlier_dir
# Filenames are of form: $splicing_outlier_calling_dir$tissue_name$splicing_outlier_suffix"_emperical_pvalue.txt"
splicing_outlier_suffix="_covariate_method_none_merged"





#############################################################
#Used Directories (directories need to be created and empty before starting)
#############################################################

# Root directories that all other directories will be placed in 
output_root="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/enrichment_analysis/"

# Directory containing processed/filtered rare variant calls
rare_variant_dir=$output_root"processed_rare_variants/"

# Directory containing results of enrichments of rare variants within outlier calls
variant_enrichment_dir=$output_root"variant_enrichment/"





#############################################################
# Scripts: Run the following X parts in Series
#############################################################




#################
# Part 1: Map variants to clusters
# We have variant calls ($variant_bed_file) for all European Ancestry individuals
# We now map these variant calls to clusters if that variant is in a window (of a range of sizes) around a splice site in that cluster
if false; then
sbatch map_variants_to_clusters.sh $variant_bed_file $cluster_info_file $rare_variant_dir
fi

if false; then
sh variant_enrichment_shell.sh $rare_variant_dir $variant_enrichment_dir $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $tissue_names_file
fi


