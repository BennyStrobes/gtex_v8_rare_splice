#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1

vcf_file="$1"
junction_count_input_dir="$2"
gtex_watershed_file="$3"
gtex_lymphocyte_jxn_count_file="$4"
processed_data_dir="$5"



################################
# Create list of individuals that have both RNA-seq and WGS
################################
individual_list=$processed_data_dir"individuals_with_WGS_and_RNA.txt"
individual_to_junction_file_list=$processed_data_dir"individual_id_to_junction_file.txt"
if false; then
python extract_individuals_with_WGS_and_RNA.py $vcf_file $junction_count_input_dir $individual_list $individual_to_junction_file_list
fi

################################
# Extract list of positions that we have gtex variants
################################
gtex_variant_position_file=$processed_data_dir"gtex_variant_positions.txt"
if false; then
python extract_gtex_variant_positions.py $gtex_watershed_file $gtex_variant_position_file
fi

################################
# Process VCF FILE
################################
variant_prefix=$processed_data_dir"snv"
if false; then
vcftools --gzvcf $vcf_file --out $variant_prefix --remove-filtered-all --keep $individual_list --positions $gtex_variant_position_file --remove-indels --max-missing-count 10 --extract-FORMAT-info DS
vcftools --gzvcf $vcf_file --out $variant_prefix --remove-filtered-all --keep $individual_list --positions $gtex_variant_position_file --remove-indels --max-missing-count 10 --freq
fi



################################
# Create junction count matrix file for amish cohort based on gtex lcl junctions
################################
amish_junction_count_file=$processed_data_dir"amish_junction_count.txt"

python generate_amish_junction_count_file.py $gtex_lymphocyte_jxn_count_file $individual_to_junction_file_list $amish_junction_count_file