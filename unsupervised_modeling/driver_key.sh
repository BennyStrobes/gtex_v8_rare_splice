######## Ben Strober
######## 11/27/18


#######################################
# Input Data (Assumes "outlier_calling" (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/outlier_calling) and enrichment analysis (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/enrichment_analysis) has been already run)
#######################################


# List of european ancestry individuals (extracted from Nicole Ferraro's variant list: "all_rare_variants_SNPs_10kb_genebody.txt")
european_ancestry_individual_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/european_ancestry_individuals.txt"

# File containing genomic annotations describing rare variants
# Downloaded from stroberb-420@scp.nygenome.org:/data/delivery/gtex-rare/data/gtex_v8_vep88_loftee_cadd_gnomad.tsv.gz on 1/9/19
raw_genomic_annotation_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/gtex_v8_rare_indivgeno_vep88_loftee_cadd_gnomad.tsv"

# File containing mapping of rare variants to genes
rare_variant_to_gene_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/all_rare_variants_SNPs_10kb_genebody.txt"

# File containing all rare variants
variant_bed_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/all_rare_variants_noWindow_SNPs.txt"

# Total expression outlier file
total_expression_outlier_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/outlier_calls/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt"

# Allele specific expression outlier file
ase_outlier_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/unsupervised_learning_input/ase_v8_data_from_jonah/median.ad.scores.uncorrected.no.global.outliers.tsv.gz"
ase_old_outlier_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/outlier_calls/median_uncorrected_DOT_scores.tsv"

# Splicing outlier file
splicing_outlier_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/splicing_outlier_calls/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt"

# Directory containing processed genomic annotations
genomic_annotation_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/genomic_annotations/"

# File containing genomic annotations
genomic_annotation_file=$genomic_annotation_dir"filtered_real_valued_gene_level_variant_compressed_genomic_annotations.txt"

# File containing genomic annotations for each variant
variant_level_genomic_annotation_file=$genomic_annotation_dir"filtered_real_valued_gene_level_clean_genomic_annotations.txt"

# File containing mapping from (gene, individual) to rare variants mapped within 10KB of the gene
gene_individual_to_variant_mapping_file=$genomic_annotation_dir"variants_compressed_onto_genes.txt"

# Ordered list of GTEx v8 tissue names
tissue_names_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtex_v8_tissues.txt"

# Directory containing tbt splicing outlier calls
splicing_outlier_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/splicing_outlier_calls/"

# Directory containing tbt ase outlier calls
ase_outlier_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/unsupervised_learning_input/ase_v8_data_from_jonah/"

# Directory containing tbt te outlier calls
te_outlier_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/unsupervised_learning_input/total_expression_data_from_nicole/"

# File containing names of genomic annotations
genomic_annotations_names_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/unsupervised_learning_input/feature_names.txt"

# File containing mapping from GTEx tissue name to chromHMM cell type (provided by Nicole Ferraro)
chrom_hmm_to_tissue_mapping_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/tissue_specific_chromHMM_annotations/roadmap_gtex_map_castel.txt"

# File containing mapping from GTEx v8 tissue names to tissue colors
tissue_colors_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtex_colors.txt"

# File in RIVER package containing simulated data
river_sim_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/sim_river.txt"

# File containing MFASS results
mfass_file="/work-zfs/abattle4/lab_data/genomic_annotation_data/MFASS/processed_data/snv/snv_data_clean.txt"

#############################################################
#Used Directories (directories need to be created and empty before starting)
#############################################################

# Root directories that all other directories will be placed in 
output_root="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/"


# Directory containing unsupervised learning input files
unsupervised_learning_input_dir=$output_root"unsupervised_learning_input/"

# Directory containing results from single outlier type RIVER runs
river_run_dir=$output_root"river/"

# Directory containing results from watershed ROC analysis
watershed_3_class_roc_run_dir=$output_root"watershed_three_class_roc/"

# Directory containing results from watershed tbt ROC analysis
watershed_tbt_roc_run_dir=$output_root"watershed_tbt_roc/"

# Directory containing results from watershed analysis applied to all variants
watershed_3_class_score_run_dir=$output_root"watershed_three_class_scores/"

# Directory containing visualizations
watershed_visualization_dir=$output_root"visualize_watershed/"

# Directory comparing to MFASS results
mfass_comparison_dir=$output_root"mfass_comparison/"

# Directory containing visualization to debug watershed
watershed_debug_visualization_dir=$output_root"visualize_debug_watershed/"

# Directory containing processed data for github repo 
github_repo_dir=$output_root"github_repo_data/"


###############################################
# Scripts 
###############################################
if false; then
sh prepare_input_files_for_unsupervised_learning_methods.sh $genomic_annotation_file $variant_level_genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $gene_individual_to_variant_mapping_file $splicing_outlier_dir $ase_outlier_dir $te_outlier_dir $tissue_names_file $ase_old_outlier_file
fi



pseudocount="30"
n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
gene_thresh="0.01"
if false; then
sh watershed_roc_run_3_outlier_types.sh $unsupervised_learning_input_dir $watershed_3_class_roc_run_dir $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $gene_thresh
fi

n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
input_file=$unsupervised_learning_input_dir"fully_observed_merged_outliers_0.05_genes_intersection_between_te_ase_splicing_features_filter_no_tissue_anno_same_N2_pairs_as_standard_3.txt"
output_stem=$watershed_3_class_roc_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_0.05_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_gene_theshold_comparison_seperate_pseudoc"
if false; then
sbatch watershed_roc_run_3_outlier_types_comparison.sh $unsupervised_learning_input_dir $watershed_3_class_roc_run_dir $n2_pair_pvalue_fraction $binary_pvalue_threshold $input_file $output_stem
fi

n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
input_file=$unsupervised_learning_input_dir"fully_observed_merged_outliers_0.1_genes_intersection_between_te_ase_splicing_features_filter_no_tissue_anno_same_N2_pairs_as_standard_3.txt"
output_stem=$watershed_3_class_roc_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_0.1_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_gene_theshold_comparison_seperate_pseudoc"
if false; then
sbatch watershed_roc_run_3_outlier_types_comparison.sh $unsupervised_learning_input_dir $watershed_3_class_roc_run_dir $n2_pair_pvalue_fraction $binary_pvalue_threshold $input_file $output_stem
fi

n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
input_file=$unsupervised_learning_input_dir"fully_observed_merged_outliers_0.01_genes_union_between_te_ase_splicing_features_filter_no_tissue_anno_same_N2_pairs_as_standard_3.txt"
output_stem=$watershed_3_class_roc_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_union_intersection_comparison_seperate_pseudoc"
if false; then
sbatch watershed_roc_run_3_outlier_types_comparison.sh $unsupervised_learning_input_dir $watershed_3_class_roc_run_dir $n2_pair_pvalue_fraction $binary_pvalue_threshold $input_file $output_stem
fi

pseudocount="30"
pvalue_fraction=".01"
binary_pvalue_threshold=".01"
if false; then
sbatch watershed_score_run.sh $unsupervised_learning_input_dir $watershed_3_class_score_run_dir $pseudocount $pvalue_fraction $binary_pvalue_threshold $watershed_3_class_roc_run_dir

fi


#####################
# TBT Model
#####################
pseudocount="10"
n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
phi_method="fixed"  # fixed, sample_size, marginal
lambda_init=".001"
lambda_pair_init=".001"
independent_variables="false"  # false or true
inference_method="pseudolikelihood" # pseudolikelihood or exact
outlier_type="total_expression"  # splicing, total_expression, ase
number_of_dimensions="49"
if false; then
sh watershed_roc_run_tbt.sh $unsupervised_learning_input_dir $watershed_tbt_roc_run_dir $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $phi_method $lambda_init $lambda_pair_init $independent_variables $inference_method $outlier_type $number_of_dimensions
fi


pseudocount="10"
n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
phi_method="fixed"  # fixed, sample_size
lambda_init=".001"
lambda_pair_init=".001"
independent_variables="true"  # false or true
inference_method="exact" # pseudolikelihood or exact
outlier_type="total_expression"  # splicing, total_expression, ase
number_of_dimensions="49"
if false; then
sh watershed_roc_run_tbt.sh $unsupervised_learning_input_dir $watershed_tbt_roc_run_dir $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $phi_method $lambda_init $lambda_pair_init $independent_variables $inference_method $outlier_type $number_of_dimensions
fi

Rscript visualize_watershed_results.R $watershed_3_class_roc_run_dir $watershed_tbt_roc_run_dir $genomic_annotations_names_file $tissue_names_file $chrom_hmm_to_tissue_mapping_file $watershed_visualization_dir $tissue_colors_file


































#############
# Old anylsis
# No longer used
###############
if false; then
sh organize_data_for_github_repo.sh $river_sim_file $github_repo_dir 
sh compare_to_mfass.sh $watershed_3_class_score_run_dir"fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors.txt.gz" $mfass_file $mfass_comparison_dir
sh river_copy_run.sh $unsupervised_learning_input_dir $river_run_dir
fi




