#!/bin/bash -l

#SBATCH
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


variant_bed_file="$1"
variant_bed_no_concensus_file="$2"
branchpoint_bed_file="$3"
splicing_outlier_dir="$4"
splicing_outlier_suffix="$5"
european_ancestry_individual_list="$6"
branchpoint_enrichment_dir="$7"
visualize_branchpoint_enrichment_dir="$8"
liftover_directory="$9"
cluster_info_file="${10}"
tissue_names_file="${11}"




####################################################
# Liftover $branchpoint_bed_file from hg19 to hg38
####################################################
branchpoint_hg38_bed_file=$branchpoint_enrichment_dir"branchpoint_hg38.bed"
# python liftover_shell.py $branchpoint_bed_file $branchpoint_hg38_bed_file $liftover_directory "hg19_to_hg38"




cluster_level_outlier_file=$splicing_outlier_dir"cross_tissue"$splicing_outlier_suffix"_emperical_pvalue.txt"

####################################################
# Map variants to branchpoints and then compute enrichment of outliers within branchpoints
####################################################
# Range of Distances
distances=( "0" "1" "2" "3" "4" "5" "10" "20")
distances=( "0" )


# Loop through distances
for distance in "${distances[@]}"
do
	# Map variants to branchpoints
	variant_bed_branchpoint_mapped=$branchpoint_enrichment_dir"all_rare_variants_noWindow_SNPs_branchpoint_"$distance"_mapped.txt"
	#python map_variants_to_branchpoints.py $variant_bed_file $branchpoint_hg38_bed_file $variant_bed_branchpoint_mapped $distance

	pvalue="1e-05"
	tbt_enrichment_file=$branchpoint_enrichment_dir"tissue_by_tissue_"$pvalue"_branchpoint_"$distance"_enrichment.txt"
	python branchpoint_enrichment_tbt.py $tissue_names_file $splicing_outlier_dir $splicing_outlier_suffix $variant_bed_branchpoint_mapped $tbt_enrichment_file $european_ancestry_individual_list $pvalue

	# Compute enrichment of outliers within branchpoints
	cross_tissue_enrichment_file=$branchpoint_enrichment_dir"cross_tissue_branchpoint_"$distance"_enrichment.txt"
	python branchpoint_enrichment.py $cluster_level_outlier_file $variant_bed_branchpoint_mapped $cross_tissue_enrichment_file $european_ancestry_individual_list
done


####################################################
# Map variants to branchpoints and then compute enrichment of outliers within branchpoints
####################################################
# Range of Distances
distances=( "0" "1" "2" "3" "4" "5" "10" "20")
distances=( "0" "1" "2" "3" "4" "5" "10" )


# Loop through distances
if false; then
for distance in "${distances[@]}"
do
	echo $distance
	# Map variants to branchpoints
	variant_bed_branchpoint_mapped=$branchpoint_enrichment_dir"all_rare_variants_no_concensus_noWindow_SNPs_branchpoint_"$distance"_mapped.txt"
	python map_variants_to_branchpoints.py $variant_bed_no_concensus_file $branchpoint_hg38_bed_file $variant_bed_branchpoint_mapped $distance

	pvalue="1e-05"
	tbt_enrichment_file=$branchpoint_enrichment_dir"tissue_by_tissue_"$pvalue"_no_concensus_branchpoint_"$distance"_enrichment.txt"
	python branchpoint_enrichment_tbt.py $tissue_names_file $splicing_outlier_dir $splicing_outlier_suffix $variant_bed_branchpoint_mapped $tbt_enrichment_file $european_ancestry_individual_list $pvalue

	# Compute enrichment of outliers within branchpoints
	cross_tissue_enrichment_file=$branchpoint_enrichment_dir"cross_tissue_no_concensus_branchpoint_"$distance"_enrichment.txt"
	python branchpoint_enrichment.py $cluster_level_outlier_file $variant_bed_branchpoint_mapped $cross_tissue_enrichment_file $european_ancestry_individual_list
done
fi

####################################################
# Visualize results
####################################################
if false; then
Rscript visualize_branchpoint_enrichment.R $branchpoint_enrichment_dir $visualize_branchpoint_enrichment_dir
fi






