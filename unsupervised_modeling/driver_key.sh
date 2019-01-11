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





#############################################################
#Used Directories (directories need to be created and empty before starting)
#############################################################

# Root directories that all other directories will be placed in 
output_root="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/"

# Directory containing processed genomic annotations
genomic_annotation_dir=$output_root"genomic_annotation/"







sh process_genomic_annotations.sh $raw_genomic_annotation_file $variant_bed_file $rare_variant_to_gene_file $genomic_annotation_dir