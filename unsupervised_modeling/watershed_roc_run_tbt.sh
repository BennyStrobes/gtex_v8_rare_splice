#!/bin/bash -l

#SBATCH
#SBATCH --time=24:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1




unsupervised_learning_input_dir="$1"
watershed_run_dir="$2"
pseudocount="$3"
pvalue_fraction="$4"
gradient_descent_threshold="$5"
theta_pair_init="$6"
lambda="$7"
lambda_pair="$8"
gradient_descent_stepsize="$9"
vi_step_size="${10}"
vi_thresh="${11}"

# Parameters!
number_of_dimensions="49"

gene_thresh="0.01"
input_stem="splicing_tbt_outliers_"$gene_thresh"_genes_intersection_between_te_ase_splicing"
input_file=$unsupervised_learning_input_dir$input_stem"_features_filter_N2_pairs.txt"
output_stem=$watershed_run_dir"splicing_tbt_intersect_te_ase_splicing_out_gene_pvalue_"$gene_thresh"_out_fraction_"$pvalue_fraction"_pseudocount_"$pseudocount"_grad_thresh_"$gradient_descent_threshold"_theta_pair_i_"$theta_pair_init"_lambda_"$lambda"_lambda_pair_"$lambda_pair"__stepsize_"$gradient_descent_stepsize$vi_step_size$vi_thresh"_log_reg_init2"

echo $pseudocount
echo $pvalue_fraction
echo $gene_thresh
echo $gradient_descent_threshold
echo $theta_pair_init
echo $lambda
echo $lambda_pair
echo $gradient_descent_stepsize
echo $vi_step_size
echo $vi_thresh

Rscript watershed_roc_tbt.R $pvalue_fraction $input_file $output_stem $number_of_dimensions $pseudocount $gradient_descent_threshold $theta_pair_init $lambda_pair $gradient_descent_stepsize $lambda $vi_step_size $vi_thresh
