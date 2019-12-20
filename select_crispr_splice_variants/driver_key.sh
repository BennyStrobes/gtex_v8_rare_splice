######## Ben Strober
######## 12/14/19




#######################################
# Input Data (Assumes everything else has previously been run)
#######################################

# Outlier calls at cluster level in Whole Blood (global outliers and EA outliers removed)
whole_blood_outlier_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/splicing_outlier_calls/Whole_Blood_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue.txt"

# Outlier calls at cluster level in Median tissue sense (global outliers and EA outliers removed)
xt_outlier_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/splicing_outlier_calls/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue.txt"

# File containing mapping from clusterID to to genes and junctions
cluster_info_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/cluster_info.txt"

# TBT junction usage file
tbt_junction_usage_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/enrichment_analysis/jxn_usage_nearby_altered_ss_enrichment/tissue_by_tissue_outliers_with_rv_in_concensus_sites_outlier_individuals_1e-05_inlier_individuals_0.5_with_read_counts.txt"

# XT Junction usage file
xt_junction_usage_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/enrichment_analysis/jxn_usage_nearby_altered_ss_enrichment/cross_tissue_outliers_with_rv_in_concensus_sites_outlier_individuals_1e-05_inlier_individuals_1e-05_with_median_read_counts.txt"

# Watershed posterior file
watershed_posterior_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/watershed_three_class_scores/fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors.txt.gz"

# File containing exons of genes relevent to our analysis
# File created in outlier calling
exon_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/gencode_v26_exons.txt"

# Inlier variant position file
inlier_variant_position_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/enrichment_analysis/variant_position_enrichment/inlier_distance_to_observed_splice_site_distance_1000_pvalue_thresh_0.001.txt"


# Inlier variant position file
outlier_variant_position_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/enrichment_analysis/variant_position_enrichment/outlier_distance_to_observed_splice_site_distance_1000_pvalue_thresh_1e-05.txt"


# Directory containing filtered junction read counts
filtered_cluster_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/"

# File containing genomic annotations describing rare variants
# Downloaded from stroberb-420@scp.nygenome.org:/data/delivery/gtex-rare/data/gtex_v8_vep88_loftee_cadd_gnomad.tsv.gz on 1/9/19
raw_genomic_annotation_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/gtex_v8_rare_indivgeno_vep88_loftee_cadd_gnomad.tsv"

#######################################
# Output Directories
#######################################

processed_data_output_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/select_crispr_splice_variants/processed_data/"





#######################################
# Scripts
#######################################


####################
# Select test variants on consensus exon splice sites
#####################
if false; then
python select_test_variants.py $cluster_info_file $whole_blood_outlier_file $xt_outlier_file $tbt_junction_usage_file $xt_junction_usage_file $watershed_posterior_file $exon_file $raw_genomic_annotation_file $processed_data_output_dir
fi

clusters_to_plot_file=$processed_data_output_dir"clusters_to_plot.txt"
tissue_name="Whole_Blood"
tissue_specific_junction_file=$filtered_cluster_dir$tissue_name"_filtered_jxns_cross_tissue_clusters_gene_mapped.txt"
if false; then
python extract_cluster_counts.py $clusters_to_plot_file $tissue_name $tissue_specific_junction_file $processed_data_output_dir
fi

if false; then
Rscript visualize_cluster_distribution.R $clusters_to_plot_file $tissue_name $exon_file $processed_data_output_dir
fi

####################
# Select test variants nearby consensus exonic splice sites
#####################
if false; then
python select_test_variants_near_splice_sites.py $cluster_info_file $whole_blood_outlier_file $xt_outlier_file $watershed_posterior_file $exon_file $outlier_variant_position_file $processed_data_output_dir
fi
clusters_to_plot_file=$processed_data_output_dir"test_variants_near_splice_sites_clusters_to_plot.txt"
tissue_name="Whole_Blood"
tissue_specific_junction_file=$filtered_cluster_dir$tissue_name"_filtered_jxns_cross_tissue_clusters_gene_mapped.txt"
if false; then
python extract_cluster_counts.py $clusters_to_plot_file $tissue_name $tissue_specific_junction_file $processed_data_output_dir"test_near_splice_sites_"
fi
if false; then
Rscript visualize_cluster_distribution.R $clusters_to_plot_file $tissue_name $exon_file $processed_data_output_dir"test_near_splice_sites_"
fi
####################
# Select control variants
#####################
if false; then
python select_background_variants.py $cluster_info_file $whole_blood_outlier_file $xt_outlier_file $watershed_posterior_file $exon_file $inlier_variant_position_file $processed_data_output_dir
fi
clusters_to_plot_file=$processed_data_output_dir"background_clusters_to_plot.txt"
tissue_name="Whole_Blood"
tissue_specific_junction_file=$filtered_cluster_dir$tissue_name"_filtered_jxns_cross_tissue_clusters_gene_mapped.txt"
if false; then
python extract_cluster_counts.py $clusters_to_plot_file $tissue_name $tissue_specific_junction_file $processed_data_output_dir"background_"
fi

if false; then
Rscript visualize_cluster_distribution.R $clusters_to_plot_file $tissue_name $exon_file $processed_data_output_dir"background_"
fi