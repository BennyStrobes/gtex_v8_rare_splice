#!/bin/bash -l

#SBATCH
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


variant_bed_file="$1"
variant_bed_no_concensus_file="$2"
ptbp_bed_file="$3"
splicing_outlier_dir="$4"
splicing_outlier_suffix="$5"
european_ancestry_individual_list="$6"
ptbp_enrichment_dir="$7"
visualize_ptbp_enrichment_dir="$8"
liftover_directory="$9"
cluster_info_file="${10}"
tissue_names_file="${11}"


####################################################
# Liftover $branchpoint_bed_file from hg18 to hg38
####################################################
ptbp_hg38_bed_file=$ptbp_enrichment_dir"ptbp_hg38.bed"
# python liftover_shell.py $ptbp_bed_file $ptbp_hg38_bed_file $liftover_directory "hg18_to_hg38"


cluster_level_outlier_file=$splicing_outlier_dir"cross_tissue"$splicing_outlier_suffix"_emperical_pvalue.txt"



####################################################
# Map variants to branchpoints and then compute enrichment of outliers within ptbp 
####################################################

# Map variants to branchpoints
variant_bed_ptbp_mapped=$ptbp_enrichment_dir"all_rare_variants_noWindow_SNPs_ptbp_mapped.txt"
# python map_variants_to_bed_file.py $variant_bed_file $ptbp_hg38_bed_file $variant_bed_ptbp_mapped 

pvalue="1e-05"
tbt_enrichment_file=$ptbp_enrichment_dir"tissue_by_tissue_"$pvalue"_ptbp_enrichment.txt"
# python bed_file_enrichment_tbt.py $tissue_names_file $splicing_outlier_dir $splicing_outlier_suffix $variant_bed_ptbp_mapped $tbt_enrichment_file $european_ancestry_individual_list $pvalue

# Compute enrichment of outliers within branchpoints
cross_tissue_enrichment_file=$ptbp_enrichment_dir"cross_tissue_ptbp_enrichment.txt"
# python bed_file_enrichment.py $cluster_level_outlier_file $variant_bed_ptbp_mapped $cross_tissue_enrichment_file $european_ancestry_individual_list

####################################################
# Map variants to branchpoints and then compute enrichment of outliers within ptbp
####################################################

# Map variants to branchpoints
variant_bed_ptbp_mapped=$ptbp_enrichment_dir"all_rare_variants_no_concensus_noWindow_SNPs_ptbp_mapped.txt"
python map_variants_to_bed_file.py $variant_bed_no_concensus_file $ptbp_hg38_bed_file $variant_bed_ptbp_mapped 

pvalue="1e-05"
tbt_enrichment_file=$ptbp_enrichment_dir"tissue_by_tissue_no_concensus_"$pvalue"_ptbp_enrichment.txt"
python bed_file_enrichment_tbt.py $tissue_names_file $splicing_outlier_dir $splicing_outlier_suffix $variant_bed_ptbp_mapped $tbt_enrichment_file $european_ancestry_individual_list $pvalue

# Compute enrichment of outliers within branchpoints
cross_tissue_enrichment_file=$ptbp_enrichment_dir"cross_tissue_no_concensus_ptbp_enrichment.txt"
python bed_file_enrichment.py $cluster_level_outlier_file $variant_bed_ptbp_mapped $cross_tissue_enrichment_file $european_ancestry_individual_list

