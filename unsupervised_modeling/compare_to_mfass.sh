#!/bin/bash -l

#SBATCH
#SBATCH --time=14:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


watershed_scores="$1"
mfass_file="$2"
mfass_comparison_dir="$3"


python process_mfass_comparison.py $watershed_scores $mfass_file $mfass_comparison_dir