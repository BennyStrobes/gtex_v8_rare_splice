######## Ben Strober
######## 11/27/18




#######################################
# Input Data (Assumes "outlier_calling" (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/outlier_calling) has been already run)
#######################################


# Ordered list of GTEx v8 tissue names
tissue_names_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtex_v8_tissues.txt"

# File containing mapping from GTEx v8 tissue names to tissue colors
tissue_colors_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtex_colors.txt"

# File containing genomic annotations describing rare variants
# Downloaded from stroberb-420@scp.nygenome.org:/data/delivery/gtex-rare/data/gtex_v8_vep88_loftee_cadd_gnomad.tsv.gz on 1/9/19
raw_genomic_annotation_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/gtex_v8_rare_indivgeno_vep88_loftee_cadd_gnomad.tsv"

# Gencode Gene annotation file
# Downloaded from https://www.gencodegenes.org/human/release_26.html on 11/26/18
gencode_gene_annotation_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gencode.v26.annotation.gff3.gz"

# List of individuals used in each tissue generated by Nicole Ferraro
# Downloaded from scp.nygenome.org:/data/delivery/gtex-rare/data/gtexV8_inds_by_tissue.txt on 11/30/18
individual_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtexV8_inds_by_tissue.txt"

# List of european ancestry individuals (extracted from Nicole Ferraro's variant list: "all_rare_variants_SNPs_10kb_genebody.txt")
european_ancestry_individual_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/european_ancestry_individuals.txt"

# List of rare variants across all European ancestry individuals (Created by Nicole Ferraro. Downloaded from stroberb-420@scp.nygenome.org:/data/delivery/gtex-rare/data/all_rare_variants_noWindow_SNPs.txt on 12/5/18)
variant_bed_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/all_rare_variants_noWindow_SNPs_refalt.txt"

# List of rare variants across all European ancestry individuals mapped to gene bodies (Created by Nicole Ferraro. Downloaded from stroberb-420@scp.nygenome.org:/data/delivery/gtex-rare/data/all_rare_variants_noWindow_SNPs.txt on 12/5/18)
variant_bed_gene_mapped_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/all_rare_variants_SNPs_10kb_genebody.txt"

# File containing all cluster_ids and their corresponding junction positions
# Generated with "outlier_calling" (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/outlier_calling)
cluster_info_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/cluster_info.txt"

# Directory containing splicing outlier calling results
# By running "outlier_calling" (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/outlier_calling)
splicing_outlier_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/splicing_outlier_calls/"

# Suffix for tissue-specific file_names in $splicing_outlier_dir
# Filenames are of form: $splicing_outlier_calling_dir$tissue_name$splicing_outlier_suffix"_emperical_pvalue.txt"
splicing_outlier_suffix="_covariate_method_none_no_global_outliers_ea_only"

# Suffix for tissue-specific file_names in $splicing_outlier_dir
# Filenames are of form: $splicing_outlier_calling_dir$tissue_name$splicing_outlier_suffix"_emperical_pvalue.txt"
splicing_outlier_include_global_outliers_suffix="_covariate_method_none"

# Directory containing heuristic outlier calls
heuristic_outlier_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/heuristic_outlier_calls/"

# Suffix for tissue-specific file names in heuristic_outlier_dir
heuristic_outlier_suffix="_heuristic_approach_outlier_calls.txt"

# File containing exons of genes relevent to our analysis
# File created in outlier calling
exon_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/gencode_v26_exons.txt"

# Directory containing filtered junction read counts
filtered_cluster_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/"

# File containing high-confidence human branchpoints according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4315302/
# Their analysis was done in K562 cells
# Downloaded on 1/7/19
# It is in hg19 reference. Need to liftover!
branchpoint_bed_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/supp_gr.182899.114_Supplemental_DataS2.bed"

# Directory containing liftover executable
liftover_directory="/work-zfs/abattle4/bstrober/tools/liftOver_x86/"

# Polypyrimidine tract binding protein (PTBP) bed file in hg18 coordinates
# Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19323 on 1/8/19
ptbp_bed_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/GSE19323_ptb_peak_cluster.bed"




#############################################################
#Used Directories (directories need to be created and empty before starting)
#############################################################

# Root directories that all other directories will be placed in 
output_root="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/enrichment_analysis/"

# Directory containing processed/filtered rare variant calls
rare_variant_dir=$output_root"processed_rare_variants/"

# Directory containing results of enrichments of rare variants within outlier calls
variant_enrichment_dir=$output_root"variant_enrichment/"

# Directory containing visualizations of results of enrichments of rare variants within outlier calls
visualize_variant_enrichment_dir=$output_root"visualize_variant_enrichment/"

# Directory containing results of variant position analysis (ie distance between RV and splice sites)
variant_position_enrichment_dir=$output_root"variant_position_enrichment/"

# Directory containing visualizations of variant position analysis (ie distance between RV and splice sites)
visualize_variant_position_enrichment_dir=$output_root"visualize_variant_position_enrichment/"

# Directory containing comparison of jxn usage nearby altered splice sites in outliers and non-outliers
jxn_usage_nearby_altered_ss_enrichment_dir=$output_root"jxn_usage_nearby_altered_ss_enrichment/"

# Directory containing visualizations of comparison of jxn usage nearby altered splice sites in outliers and non-outliers
visualize_jxn_usage_nearby_altered_ss_enrichment_dir=$output_root"visualize_jxn_usage_nearby_altered_ss_enrichment/"

# Directory containing visualizations of cluster distributions for outliers compared to non-outliers
visualize_cluster_distribution_dir=$output_root"visualize_cluster_distribution/"

# Directory containing data for branchpoint enrichments
branchpoint_enrichment_dir=$output_root"branchpoint_enrichment/"

# Directory containing visualizations for branchpoint enrichments
visualize_branchpoint_enrichment_dir=$output_root"visualize_branchpoint_enrichment/"






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

# Part 1b: Map variants to jxns
# We have variant calls ($variant_bed_file) for all European Ancestry individuals
# We now map these variant calls to jxns if that variant is in a window (of a range of sizes) around a junction (this is used for heuristic outlier analysis)
if false; then
sh map_variants_to_junctions.sh $variant_bed_file $heuristic_outlier_dir $heuristic_outlier_suffix $tissue_names_file $rare_variant_dir
fi


#################
# Part 2: Compute enrichment of rare variants within spliding outlier calls
# Do this enrichments:
#    1. For each of the tissues, independently
#    2. For cross tissue outliers (median pvalue)
# Compute enrichments of rare variants within heuristic outilers
# Then visualize enrichments
if false; then
sh variant_enrichment_shell.sh $rare_variant_dir $variant_enrichment_dir $splicing_outlier_dir $splicing_outlier_suffix $splicing_outlier_include_global_outliers_suffix $european_ancestry_individual_list $tissue_names_file $visualize_variant_enrichment_dir $tissue_colors_file $heuristic_outlier_dir $heuristic_outlier_suffix
fi


#################
# Part 3: Compare distances between variants and splice sites for outliers vs non-outliers
# Then visualize results
if false; then
sh variant_position_enrichment_shell.sh $rare_variant_dir $variant_position_enrichment_dir $visualize_variant_position_enrichment_dir $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $gencode_gene_annotation_file $cluster_info_file $exon_file
fi

#################
# Part 4: Compare jxn usage nearby altered splice sites in outliers and non-outliers
# Then visualize results
if false; then
sh junction_usage_nearby_altered_splice_site_enrichment_shell.sh $rare_variant_dir $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $gencode_gene_annotation_file $cluster_info_file $exon_file $jxn_usage_nearby_altered_ss_enrichment_dir $visualize_jxn_usage_nearby_altered_ss_enrichment_dir $tissue_names_file $filtered_cluster_dir
fi

#################
# Part 5: Visualize cluster distributions for outliers compared to non-outliers
if false; then
sh visualize_cluster_distribution_shell.sh $rare_variant_dir $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $cluster_info_file $exon_file $visualize_cluster_distribution_dir $tissue_names_file $filtered_cluster_dir
fi

#################
# Part 6: Compute enrichments in branchpoints for outliers vs non-outliers
# Then visualize results
variant_cluster_intron_mapped_file=$rare_variant_dir"variant_cluster_only_bed_100.txt"
if false; then
sh branchpoint_enrichment_shell.sh $variant_cluster_intron_mapped_file $branchpoint_bed_file $splicing_outlier_dir $splicing_outlier_suffix $european_ancestry_individual_list $branchpoint_enrichment_dir $visualize_branchpoint_enrichment_dir $liftover_directory $cluster_info_file $tissue_names_file
fi




