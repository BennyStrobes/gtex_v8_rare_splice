#!/bin/bash -l

#SBATCH
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1

raw_genomic_annotation_file="$1"
variant_bed_file="$2"
rare_variant_to_gene_file="$3"
genomic_annotation_dir="$4"



python process_genomic_annotations.py $raw_genomic_annotation_file $variant_bed_file $rare_variant_to_gene_file $genomic_annotation_dir