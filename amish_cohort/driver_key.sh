######## Ben Strober
######## 12/17/19


#######################################
# Input Data
#######################################

# VCF File for amish cohort
# Sent from Casey Brown
vcf_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/amish_cohort/chrAll.imputed.poly.vcf.gz"

# Directory containing junction count files for each sample
junction_count_input_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/amish_cohort/"

# Gtex watershsed file
gtex_watershed_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/watershed_three_class_scores/fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors.txt.gz"

# Processed junction count file for lymphocytes
gtex_lymphocyte_jxn_count_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/Cells_EBV-transformed_lymphocytes_filtered_jxns_cross_tissue_clusters_gene_mapped.txt"

#######################################
# output Directories
#######################################
output_root="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/amish_cohort/"

# Directory containing preprocessed data
processed_data_dir=$output_root"processed_data/"


sh preprocess_data.sh $vcf_file $junction_count_input_dir $gtex_watershed_file $gtex_lymphocyte_jxn_count_file $processed_data_dir