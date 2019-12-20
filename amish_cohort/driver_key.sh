######## Ben Strober
######## 12/17/19


#######################################
# Input Data
#######################################

# Sample mapping file
# sent from Casey brown
sample_mapping_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/amish_cohort/GM_Blood_IDs.txt"

# VCF File for amish cohort
# Sent from Casey Brown
vcf_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/amish_cohort/chrAll.imputed.poly.vcf.gz"

# Directory containing junction count files for each sample
junction_count_input_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/amish_cohort/"

# Gtex watershsed file
gtex_watershed_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/watershed_three_class_scores/fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors.txt.gz"

# Processed junction count file for lymphocytes
gtex_lymphocyte_jxn_count_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/Cells_EBV-transformed_lymphocytes_filtered_jxns_cross_tissue_clusters_gene_mapped.txt"

# File containing all cluster_ids and their corresponding junction positions
# Generated with "outlier_calling" (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/outlier_calling)
cluster_info_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/cluster_info.txt"

#######################################
# output Directories
#######################################
output_root="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/amish_cohort/"

# Directory containing preprocessed data
processed_data_dir=$output_root"processed_data/"

# Directory containing outlier cals
splicing_outlier_dir=$output_root"splicing_outlier_calls/"

# Directory containing outlier cals
merged_results_dir=$output_root"merged_results/"

# Directory containing outlier cals
visualize_results_dir=$output_root"visualize_results/"


########################
# Part 1 preprocess the data
#########################
if false; then
sh preprocess_data.sh $vcf_file $junction_count_input_dir $gtex_watershed_file $gtex_lymphocyte_jxn_count_file $sample_mapping_file $processed_data_dir
fi


#########################
# Part 2: Call splicing 
#########################
train_junction_file=$gtex_lymphocyte_jxn_count_file
test_junction_file=$processed_data_dir"amish_junction_count.txt"

# How many nodes to run in parallel
total_jobs="20"
# Whether to include covariates in GLM
covariate_method="none"
# Number of reads to simulate for emperical distribution
num_reads="20000"
# Type of Dirichlet multinomial to use (either 'standard' or 'no_prior')
model_version="standard"
# Max number of junctions per cluster
max_number_of_junctions_per_cluster="20"
if false; then
for job_number in $(seq 0 `expr $total_jobs - "1"`); do
	output_root=$splicing_outlier_dir"amish_cohort_covariate_method_"$covariate_method"_"$job_number"_"$total_jobs
	sbatch call_splicing_outliers_train_test.sh $test_junction_file $train_junction_file $covariate_method $max_number_of_junctions_per_cluster $output_root $num_reads $model_version $job_number $total_jobs
done
fi

if false; then
sh merge_splicing_outlier_calls.sh $splicing_outlier_dir $covariate_method $cluster_info_file $total_jobs
fi

#########################
# Part 3: Merge results
#########################
sh merge_and_visualize_results.sh $splicing_outlier_dir $processed_data_dir $merged_results_dir $visualize_results_dir

