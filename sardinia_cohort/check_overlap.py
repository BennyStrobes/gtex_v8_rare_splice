import numpy as np 
import os
import sys
import pdb
import gzip




def get_sardinia_variant_gene_pairs(sardinia_variant_file):
	f = open(sardinia_variant_file)
	variants = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_info = data[0].split(':')
		new_variant_id = variant_info[0] + '_' + variant_info[1].split('-')[0]
		variant_gene_pair = data[1] + ':' + new_variant_id
		variants[variant_gene_pair] = []
	f.close()
	return variants



def get_gtex_sardinia_overlap(sardinia_variants, gtex_watershed_file):
	f = gzip.open(gtex_watershed_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			print(data)
			continue
		# Extract relevent fields
		line_info = data[0].split(':')
		ensamble_id = line_info[1].split('.')[0]
		variant_info = line_info[2].split('_')
		variant_id = variant_info[0] + '_' + variant_info[1]
		test_name = ensamble_id + ':' + variant_id
		# Ignore lines not in sardinia
		if test_name not in sardinia_variants:
			continue
		sardinia_variants[test_name].append(np.asarray(data))
	f.close()
	return sardinia_variants

#####################
# Command line args
#####################
sardinia_variant_file = sys.argv[1]
gtex_watershed_file = sys.argv[2]
overlap_file = sys.argv[3]



#####################
# Extract sardinia variants
sardinia_variants = get_sardinia_variant_gene_pairs(sardinia_variant_file)

#####################
# Extract GTEx variants
variants = get_gtex_sardinia_overlap(sardinia_variants, gtex_watershed_file)



t = open(overlap_file, 'w')

t.write('variant_id\tmedian_expression_pvalue\tmedian_ase_pvalue\tmedian_splicing_pvalue\tmedian_expression_watershed_posterior\tmedian_ase_watershed_posterior\tmedian_splicing_watershed_posterior\tnum_expression\tnum_ase\tnum_splicing\n')

for variant in variants.keys():
	if len(variants[variant]) == 0:
		continue
	eOutliers = []
	sOutliers = []
	aseOutliers = []
	e_watersheds = []
	s_watersheds = []
	ase_watersheds = []
	for info in variants[variant]:
		if info[2] != 'NaN':
			eOutliers.append(np.abs(float(info[2])))
			e_watersheds.append(float(info[11]))
		if info[1] != 'NaN':
			sOutliers.append(np.abs(float(info[1])))
			s_watersheds.append(float(info[10]))
		if info[3] != 'NaN':	
			aseOutliers.append(np.abs(float(info[3])))
			ase_watersheds.append(float(info[12]))
	t.write(variant + '\t' + str(np.median(eOutliers)) + '\t' + str(np.median(aseOutliers)) + '\t' + str(np.median(sOutliers)) + '\t')
	t.write(str(np.median(e_watersheds)) + '\t' + str(np.median(ase_watersheds)) + '\t' + str(np.median(s_watersheds)) + '\t')
	t.write(str(np.median(sOutliers)) + '\t' + str(len(eOutliers)) + '\t' + str(len(aseOutliers)) + '\t' + str(len(sOutliers)) + '\n')

t.close()

