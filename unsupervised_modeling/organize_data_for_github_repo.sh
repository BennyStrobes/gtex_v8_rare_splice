#!/bin/bash -l

#SBATCH
#SBATCH --time=3:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


river_sim_file="$1"
github_repo_dir="$2"


watershed_sim_file=$github_repo_dir"watershed_simulated_data.txt"
if false; then
Rscript generate_watershed_simulation_data.R $river_sim_file $watershed_sim_file
fi
watershed_sim_file=$github_repo_dir"watershed_simulated_data_v3.txt"
river_sim_file_stem=$github_repo_dir"river_simulated_"
if false; then
python generate_watershed_simulation_data_v2.py $river_sim_file $watershed_sim_file $river_sim_file_stem
fi

n2_pair_pvalue_fraction=".1"
binary_pvalue_threshold=".1"
pseudocount="10"
number_of_dimensions="3"
lambda=".01"
seed="1"
output_stem=$github_repo_dir"sim_3_class_n2_pair_outlier_fraction_"$n2_pair_pvalue_fraction"_binary_pvalue_threshold_"$binary_pvalue_threshold"_pseudocount_"$pseudocount"_lambda_"$lambda"_seed_"$seed
output_stem="sim_3_class_n2_pair_outlier_fraction_"$n2_pair_pvalue_fraction"_binary_pvalue_threshold_"$binary_pvalue_threshold"_pseudocount_"$pseudocount"_lambda_"$lambda"_seed_"$seed

river_file=$river_sim_file_stem"pheno_1.txt"
Rscript watershed_roc_3_outlier_types_git.R $watershed_sim_file $output_stem $number_of_dimensions $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $lambda $seed


echo $lambda
echo $seed

n2_pair_pvalue_fraction=".1"
binary_pvalue_threshold=".1"
pseudocount="10"
number_of_dimensions="3"
output_stem=$github_repo_dir"sim_3_class_n2_pair_outlier_fraction_"$n2_pair_pvalue_fraction"_binary_pvalue_threshold_"$binary_pvalue_threshold"_pseudocount_"$pseudocount"_lambda_"$lambda"_seed_"$seed

if false; then
Rscript watershed_roc_3_outlier_types_git.R $watershed_sim_file $output_stem $number_of_dimensions $pseudocount $n2_pair_pvalue_fraction $binary_pvalue_threshold $lambda $seed



Rscript git_repo_viz.R $output_stem
fi