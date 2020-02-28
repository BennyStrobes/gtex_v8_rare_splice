
unsupervised_learning_input_dir="$1"
watershed_run_dir="$2"
pseudocount="$3"
n2_pair_pvalue_fraction="$4"
binary_pvalue_threshold="$5"
gene_thresh="$6"
tissue_names_file="$7"






while read tissue_name; do
	input_stem="single_tissue_"$tissue_name"_pvalue_cross_signal_outliers_"$gene_thresh"_genes_intersection_between_te_ase_splicing"
 	input_file=$unsupervised_learning_input_dir$input_stem"_features_filter_N2_pairs.txt"
 	number_of_dimensions="3"
 	output_stem=$watershed_run_dir"single_tissue_cross_signal_somewhat_observed_"$tissue_name"_fully_observed_te_ase_splicing_outliers_gene_pvalue_"$gene_thresh"_n2_pair_outlier_fraction_"$n2_pair_pvalue_fraction"_binary_pvalue_threshold_"$binary_pvalue_threshold"_pseudocount_"$pseudocount"_seed_3"
 	sh run_single_tissue_cross_signal_roc.sh $input_file $output_stem $number_of_dimensions $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold
done <$tissue_names_file

