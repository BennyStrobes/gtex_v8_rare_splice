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

# File containing exons of genes relevent to our analysis
# File created in outlier calling
exon_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/gencode_v26_exons.txt"

# File containing all cluster_ids and their corresponding junction positions
# Generated with "outlier_calling" (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/outlier_calling)
cluster_info_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/cluster_info.txt"

# File containing number of rare variants per individual (stratitfied by MAF bins)
# Created by Xin (5/5/19)
num_variants_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtex_v8_variant_count.txt"

# File containing CADD scores (also containing $cadd_file".tbi")
#cadd_file="/work-zfs/abattle4/lab_data/genomic_annotation_data/hg38/whole_genome_SNVs.tsv.gz"

# File containing CADD annotations (also containing $cadd_anno_file".tbi")
#cadd_anno_file="/work-zfs/abattle4/lab_data/genomic_annotation_data/hg38/whole_genome_SNVs_inclAnno.tsv.gz"


#############################################################
#Used Directories (directories need to be created and empty before starting)
#############################################################



# Directory containing genomic annotations
genomic_annotation_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/genomic_annotations/"

# Directory containing visualizations
visualization_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/genomic_annotations/visualize/"





##########################
# Part 1: Preprocess genomic annotations from Xin
# Essentially compress genomic annotation lines from transcript level to the gene level
if false; then
sh preprocess_genomic_annotations.sh $raw_genomic_annotation_file $variant_bed_file $rare_variant_to_gene_file $genomic_annotation_dir $exon_file $cluster_info_file
fi


Rscript visualize_genomic_variants.R $num_variants_file $visualization_dir

