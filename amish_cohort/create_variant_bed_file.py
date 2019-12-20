import numpy as np 
import os
import sys
import pdb
import gzip

def extract_watershed_variants(gtex_watershed_file, valid_positions):
	f = gzip.open(gtex_watershed_file)
	head_count = 0
	var_dict = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		indi = data[0].split(':')[0]
		gene = data[0].split(':')[1]
		variant_id = data[0].split(':')[2]
		variant_id_info = variant_id.split('_')
		short_variant_id = '_'.join(variant_id_info[:4])
		really_short_variant_id = '_'.join(variant_id_info[:2])
		if really_short_variant_id not in valid_positions:
			continue
		num_variants = variant_id_info[4]
		if short_variant_id not in var_dict:
			var_dict[short_variant_id] = []
		var_dict[short_variant_id].append(line)
	return var_dict

def extract_watershed_variants_v2(gtex_watershed_file, valid_positions, outlier_file):
	f = open(outlier_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			samples = np.asarray(data[1:])
			continue
		pvals = np.asarray(data[1:])
		ensamble_id = data[0]
		for i, pval in enumerate(pvals):
			sample = samples[i]
			dicti[ensamble_id + '_' + sample] = pval
	f.close()
	f = gzip.open(gtex_watershed_file)
	head_count = 0
	var_dict = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		indi = data[0].split(':')[0]
		gene = data[0].split(':')[1]
		if gene + '_' + sample not in dicti:
			continue
		variant_id = data[0].split(':')[2]
		variant_id_info = variant_id.split('_')
		short_variant_id = '_'.join(variant_id_info[:4])
		really_short_variant_id = '_'.join(variant_id_info[:2])
		if really_short_variant_id not in valid_positions:
			continue
		num_variants = variant_id_info[4]
		data[1] = dicti[gene + '_' + sample]
		data[10] = dicti[gene + '_' + sample]
		liner = '\t'.join(data)
		if short_variant_id not in var_dict:
			var_dict[short_variant_id] = []
		var_dict[short_variant_id].append(liner)
	return var_dict


def get_valid_positions(variant_frequency_file):
	f = open(variant_frequency_file)
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		pos = data[0] + '_' + data[1]
		dicti[pos] = 1
	return dicti

def compress_watershed_across_individuals(watershed_arr):
	genes = {}
	for sample in watershed_arr:
		data = sample.split()
		if data[1] == 'NaN':
			continue
		gene_id = data[0].split(':')[1]
		watershed_posterior = float(data[10])
		num_rv = data[0].split('_')[-1]
		if num_rv == '0':
			pdb.set_trace()
		gene_id = gene_id + '_' + num_rv
		if gene_id not in genes:
			genes[gene_id] = [watershed_posterior]
		else:
			old_arr = genes[gene_id]
			old_arr.append(watershed_posterior)
			genes[gene_id] = old_arr
	arr = []
	for gene in genes:
		arr.append((gene, np.median(genes[gene])))
	return arr

def get_alternate_count(standard_order, num_rv):
	if standard_order == True:
		alternate_count = num_rv
	else:
		alternate_count = 2 - num_rv
	return alternate_count

variant_dosage_file = sys.argv[1]
variant_frequency_file = sys.argv[2]
gtex_watershed_file = sys.argv[3]
variant_bed_file = sys.argv[4]



valid_positions = get_valid_positions(variant_frequency_file)
# Extract dictionary of watershed variants
# Will be mapping from variant id to lines involving that variant
#watershed_variants = extract_watershed_variants(gtex_watershed_file, valid_positions)
watershed_variants = extract_watershed_variants(gtex_watershed_file, valid_positions)

f = open(variant_dosage_file)
g = open(variant_frequency_file)
t = open(variant_bed_file, 'w')
t.write('variant_id\tensamble_id\tamish_sample_id\tmedian_watershed_splice_posterior\n')
aa=0
bb=0
head_count = 0
for line1 in f:
	line1 = line1.rstrip()
	data1 = line1.split()
	line2 = g.next().rstrip()
	data2 = line2.split()
	if head_count == 0:
		head_count = head_count + 1
		amish_samples = np.asarray(data1[2:])
		continue
	allele1 = data2[4].split(':')[0]
	allele2 = data2[5].split(':')[0]
	var1 = data2[0] + '_' + data2[1] + '_' + allele1 + '_' + allele2
	var2 = data2[0] + '_' + data2[1] + '_' + allele2 + '_' + allele1
	# Simple error checking
	if var1 in watershed_variants and var2 in watershed_variants:
		print('assumption error')
		pdb.set_trace()
	if var1 in watershed_variants:
		var = var1
		major_allele = allele1
		variant_allele = allele2
		standard_order = True
	elif var2 in watershed_variants:
		var = var2
		major_allele = allele2
		variant_allele = allele1
		standard_order = False
	else:
		continue
	watershed_array = compress_watershed_across_individuals(watershed_variants[var])
	for ele in watershed_array:
		average_watershed_splice_posterior = ele[1]
		ensamble_id = ele[0].split('_')[0]
		num_rv = int(ele[0].split('_')[1])
		if num_rv < 1 or num_rv > 2:
			print('assumption error')
			pdb.set_trace()
		alternate_count = get_alternate_count(standard_order, num_rv)
		genotypes = np.asarray(data1[2:]).astype(float)
		hits = 0
		for index, genotype in enumerate(genotypes):
			amish_sample = amish_samples[index]
			if int(np.round(genotype)) == alternate_count:
				rv = var + '_' + str(num_rv)
				t.write(rv + '\t' + ensamble_id + '\t' + amish_sample + '\t' + str(average_watershed_splice_posterior) + '\n')
				hits = hits + 1
		if hits == len(genotypes):
			print('all rv')
			pdb.set_trace()
t.close()
g.close()
f.close()

