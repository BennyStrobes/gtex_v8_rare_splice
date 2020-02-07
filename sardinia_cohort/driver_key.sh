
#################
# Input data
#################
sardinia_variant_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/sardinia/sardinia_variant_gene_pairs.txt"
gtex_watershed_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/watershed_three_class_scores/fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors.txt.gz"




#################
# Output directory
#################
sardinia_output_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/sardinia/"





##################
# Overlap sardinia variants with gtex variants
###################
overlap_file=$sardinia_output_dir"sardinia_gtex_overlap.txt"


python check_overlap.py $sardinia_variant_file $gtex_watershed_file $overlap_file