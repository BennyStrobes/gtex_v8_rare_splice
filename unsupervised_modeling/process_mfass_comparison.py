import numpy as np 
import os
import sys
import pdb
import gzip




def create_mfass_mapping(mfass_file):
	f = open(mfass_file)
	counter = 0
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 54:
			continue
		#delta_psi = data[51]
		ensamble_id = data[1]
		ref_allele = data[10]
		alt_allele = data[11]
		snp_pos = data[18]
		chromosome = data[3]
		strong_lof = data[52]
		test_name =  chromosome + '_' + ref_allele + '_' + alt_allele + '_' + snp_pos
		if ensamble_id != 'NA' and ref_allele != 'NA' and alt_allele != 'NA' and snp_pos != 'NA' and chromosome != 'NA' and strong_lof != 'NA':
			if test_name in dicti:
				pdb.set_trace()
				print('assumption error')
			dicti[test_name] = strong_lof
	f.close()
	return dicti





######################
# command line args
#####################
watershed_score_file = sys.argv[1]
mfass_file = sys.argv[2]
output_dir = sys.argv[3]





# Create mapping from variant to mfass score
mfass_mapping = create_mfass_mapping(mfass_file)
hit1 = 0
hit2 = 0
f = gzip.open(watershed_score_file)
overlapping_variants = {}
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	test_id = data[0]
	test_info = test_id.split(':')
	ensamble_id = test_info[1].split('.')[0]
	chromosome = test_info[2].split('_')[0]
	snp_pos = test_info[2].split('_')[1]
	ref = test_info[2].split('_')[2]
	alt = test_info[2].split('_')[3]
	
	test_name1 = chromosome + '_' + ref + '_' + alt + '_' + snp_pos

	if test_name1 in mfass_mapping:
		if test_name1 not in overlapping_variants:
			overlapping_variants[test_name1] = (mfass_mapping[test_name1], float(data[10]))
		else:
			old_watershed = overlapping_variants[test_name1][1]
			overlapping_variants[test_name1] = (mfass_mapping[test_name1], max(old_watershed, float(data[10])))

