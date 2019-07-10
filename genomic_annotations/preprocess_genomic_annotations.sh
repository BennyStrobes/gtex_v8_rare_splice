#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1

raw_genomic_annotation_file="$1"
variant_bed_file="$2"
rare_variant_to_gene_file="$3"
genomic_annotation_dir="$4"
exon_file="$5"
cluster_info_file="$6"
chrom_hmm_file="$7"




python preprocess_genomic_annotations.py $raw_genomic_annotation_file $variant_bed_file $rare_variant_to_gene_file $genomic_annotation_dir $exon_file $cluster_info_file $chrom_hmm_file