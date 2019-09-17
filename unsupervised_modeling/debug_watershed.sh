#!/bin/bash -l

#SBATCH
#SBATCH --time=38:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1



watershed_3_class_roc_run_dir="$1"
watershed_3_class_score_run_dir="$2"
watershed_debug_visualization_dir="$3"


if false; then
python shorten_score_files.py $watershed_3_class_score_run_dir
fi

if false; then
python shorten_score_files_variant_level.py $watershed_3_class_score_run_dir
fi




if false; then
Rscript visualize_debug_watershed.R $watershed_3_class_roc_run_dir $watershed_3_class_score_run_dir $watershed_debug_visualization_dir
fi