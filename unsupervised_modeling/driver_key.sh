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
ase_outlier_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/outlier_calls/median_uncorrected_DOT_scores.tsv"

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
ase_outlier_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/unsupervised_learning_input/ase_data_from_jonah/"

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



###############################################
# Scripts 
###############################################
if false; then
sh prepare_input_files_for_unsupervised_learning_methods.sh $genomic_annotation_file $variant_level_genomic_annotation_file $total_expression_outlier_file $ase_outlier_file $splicing_outlier_file $unsupervised_learning_input_dir $gene_individual_to_variant_mapping_file $splicing_outlier_dir $ase_outlier_dir $tissue_names_file
fi

if false; then
sh river_copy_run.sh $unsupervised_learning_input_dir $river_run_dir
fi

if false; then
pseudocount="30"
n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
sh watershed_roc_run_3_outlier_types.sh $unsupervised_learning_input_dir $watershed_3_class_roc_run_dir $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold
fi


if false; then
pseudocount="30"
pvalue_fraction=".01"
sbatch watershed_score_run.sh $unsupervised_learning_input_dir $watershed_3_class_score_run_dir $pseudocount $pvalue_fraction
fi


#####################
# TBT Model
#####################
if false; then
pseudocount="30"
n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
phi_method="fixed"  # fixed, sample_size
lambda_init=".001"
lambda_pair_init=".001"
independent_variables="false"  # false or true
inference_method="pseudolikelihood" # pseudolikelihood or exact
outlier_type="total_expression"  # splicing, total_expression, ase
sh watershed_roc_run_tbt.sh $unsupervised_learning_input_dir $watershed_tbt_roc_run_dir $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $phi_method $lambda_init $lambda_pair_init $independent_variables $inference_method $outlier_type

pseudocount="10"
n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
phi_method="fixed"  # fixed, sample_size
lambda_init=".001"
lambda_pair_init=".001"
independent_variables="false"  # false or true
inference_method="pseudolikelihood" # pseudolikelihood or exact
outlier_type="total_expression"  # splicing, total_expression, ase
sh watershed_roc_run_tbt.sh $unsupervised_learning_input_dir $watershed_tbt_roc_run_dir $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $phi_method $lambda_init $lambda_pair_init $independent_variables $inference_method $outlier_type

pseudocount="30"
n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
phi_method="sample_size"  # fixed, sample_size
lambda_init=".001"
lambda_pair_init=".001"
independent_variables="false"  # false or true
inference_method="pseudolikelihood" # pseudolikelihood or exact
outlier_type="total_expression"  # splicing, total_expression, ase
sh watershed_roc_run_tbt.sh $unsupervised_learning_input_dir $watershed_tbt_roc_run_dir $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $phi_method $lambda_init $lambda_pair_init $independent_variables $inference_method $outlier_type



pseudocount="30"
n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
phi_method="fixed"  # fixed, sample_size
lambda_init=".001"
lambda_pair_init=".001"
independent_variables="true"  # false or true
inference_method="exact" # pseudolikelihood or exact
outlier_type="total_expression"  # splicing, total_expression, ase
sh watershed_roc_run_tbt.sh $unsupervised_learning_input_dir $watershed_tbt_roc_run_dir $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $phi_method $lambda_init $lambda_pair_init $independent_variables $inference_method $outlier_type


pseudocount="10"
n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
phi_method="fixed"  # fixed, sample_size
lambda_init=".001"
lambda_pair_init=".001"
independent_variables="true"  # false or true
inference_method="exact" # pseudolikelihood or exact
outlier_type="total_expression"  # splicing, total_expression, ase
sh watershed_roc_run_tbt.sh $unsupervised_learning_input_dir $watershed_tbt_roc_run_dir $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $phi_method $lambda_init $lambda_pair_init $independent_variables $inference_method $outlier_type

pseudocount="30"
n2_pair_pvalue_fraction=".01"
binary_pvalue_threshold=".01"
phi_method="sample_size"  # fixed, sample_size
lambda_init=".001"
lambda_pair_init=".001"
independent_variables="true"  # false or true
inference_method="exact" # pseudolikelihood or exact
outlier_type="total_expression"  # splicing, total_expression, ase
sh watershed_roc_run_tbt.sh $unsupervised_learning_input_dir $watershed_tbt_roc_run_dir $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $phi_method $lambda_init $lambda_pair_init $independent_variables $inference_method $outlier_type
fi


Rscript visualize_watershed_tbt_results.R $watershed_tbt_roc_run_dir

