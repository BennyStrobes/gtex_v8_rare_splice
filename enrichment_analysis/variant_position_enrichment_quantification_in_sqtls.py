import numpy as np 
import os
import sys
import pdb
import gzip
import random




# Get mapping from ensamble_id to strand
def get_mapping_from_ensamble_id_to_strand(gencode_gene_annotation_file):
	mapping = {}
	f = gzip.open(gencode_gene_annotation_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if line.startswith('#'):  # ignore header lines
			continue
		# Parse line
		gene_info_arr = data[8].split(';')
		hits = 0
		for ele in gene_info_arr:
			info = ele.split('=')
			if info[0] == 'gene_id':
				hits = hits + 1
				gene_name = info[1]
		strand = data[6]
		if strand != '+' and strand != '-':
			print('assumtpion erroro!')
			pdb.set_trace()
		if gene_name not in mapping:
			mapping[gene_name] = strand
		else:
			if mapping[gene_name] != strand:
				print('mapping consistency error')
				pdb.set_trace()
	return mapping

# Make dictionary list of valid chromosomes
def make_dictionary_list_of_autosomal_chromosomes():
	dicti = {}
	for chrom_int in range(1,23):
		dicti['chr' + str(chrom_int)] = 1
	return dicti


# Filter sQTL file to variant-jxn papers where variant is within a $distance_window BP around a splice site of the junction
def filter_sqtl_file_to_variants_near_ss(sqtl_file, filtered_sqtl_file, distance_window):
	# Make dictionary list of valid chromosomes
	valid_chromosomes = make_dictionary_list_of_autosomal_chromosomes()
	# Used to skip header
	head_count = 0
	# Open input and output file handles
	f = open(sqtl_file)
	t = open(filtered_sqtl_file, 'w')
	# Stream sqtl file
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		# Parse line
		jxn_info = data[0].split(':')
		chromosome = jxn_info[0]
		ss_start_pos = int(jxn_info[1])
		ss_end_pos = int(jxn_info[2])
		variant_info = data[1].split('_')
		variant_pos = int(variant_info[1])
		ref_allele = variant_info[2]
		alt_allele = variant_info[3]
		# Error checking
		if ss_end_pos < ss_start_pos:
			print('ASSUMPTIONER EROROOR')
			pdb.set_trace()
		# Filter out non autosomal tests
		if chromosome not in valid_chromosomes:
			continue
		# Filter out variants farther from $distance_window BP from SS
		if abs(variant_pos - ss_start_pos) > distance_window and abs(variant_pos - ss_end_pos) > distance_window:
			continue
		# Filter out non-SNVs
		if len(ref_allele) > 1 or len(alt_allele) > 1:
			continue
		# PASSED ALL FILTERS!
		t.write(line + '\n')
	t.close()
	f.close()

def get_distance_to_nearest_ss(ss_start_pos, ss_end_pos, variant_pos, strand):
	min_distance = 100000000000000000
	ss_type = 'null'
	disty = ss_start_pos - variant_pos
	if abs(disty) < abs(min_distance):
		min_distance = disty
		if strand == '+':
			ss_type = 'donor'
		elif strand == '-':
			ss_type = 'acceptor'
		else:
			print('assumptionerror')
			pdb.set_trace()
	disty = variant_pos - ss_end_pos
	if abs(disty) < abs(min_distance):
		min_distance = disty
		if strand == '+':
			ss_type = 'acceptor'
		elif strand == '-':
			ss_type = 'donor'
		else:
			print('asssumptioner eoror')
			pdb.set_trace()
	return min_distance, ss_type

def extract_significant_sqtl_distances(filtered_sqtl_file, pvalue_thresh, significant_output_file, ensamble_to_strand_mapping):
	# Data structure to keep track of clusters
	clusters = {}
	# To skip header
	head_count = 0
	# Stream input file
	f = open(filtered_sqtl_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1 
			continue
		# Parse line
		jxn_info = data[0].split(':')
		chromosome = jxn_info[0]
		ss_start_pos = int(jxn_info[1])
		ss_end_pos = int(jxn_info[2])
		ensamble_id = jxn_info[4]
		cluster_id = jxn_info[3]
		variant_info = data[1].split('_')
		variant_pos = int(variant_info[1])
		ref_allele = variant_info[2]
		alt_allele = variant_info[3]
		pvalue = float(data[6])
		# Get strand
		strand = ensamble_to_strand_mapping[ensamble_id]
		# Get distance to nearest ss
		distance, ss_type = get_distance_to_nearest_ss(ss_start_pos, ss_end_pos, variant_pos, strand)
		if pvalue <= pvalue_thresh:
			if cluster_id not in clusters:
				clusters[cluster_id] = []
			clusters[cluster_id].append((pvalue, distance, ss_type))
	f.close()
	t = open(significant_output_file, 'w')
	t.write('distance\tsplice_site_type\tref\talt\n')
	for cluster_id in clusters.keys():
		sig_array = clusters[cluster_id]
		sig_array.sort(key=lambda tup: tup[0])
		t.write(str(sig_array[0][1]) + '\t' + sig_array[0][2] + '\t' + ref_allele + '\t' + alt_allele + '\n')
	t.close()

def extract_background_sqtl_distances(filtered_sqtl_file, pvalue_thresh, significant_output_file, ensamble_to_strand_mapping):
	# Data structure to keep track of clusters
	clusters = {}
	# To skip header
	head_count = 0
	# Stream input file
	f = open(filtered_sqtl_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1 
			continue
		# Parse line
		jxn_info = data[0].split(':')
		chromosome = jxn_info[0]
		ss_start_pos = int(jxn_info[1])
		ss_end_pos = int(jxn_info[2])
		ensamble_id = jxn_info[4]
		cluster_id = jxn_info[3]
		variant_info = data[1].split('_')
		variant_pos = int(variant_info[1])
		ref_allele = variant_info[2]
		alt_allele = variant_info[3]
		pvalue = float(data[6])
		# Get strand
		strand = ensamble_to_strand_mapping[ensamble_id]
		# Get distance to nearest ss
		distance, ss_type = get_distance_to_nearest_ss(ss_start_pos, ss_end_pos, variant_pos, strand)
		if pvalue > pvalue_thresh:
			if cluster_id not in clusters:
				clusters[cluster_id] = []
			clusters[cluster_id].append((pvalue, distance, ss_type))
	f.close()
	t = open(significant_output_file, 'w')
	t.write('distance\tsplice_site_type\tref\talt\n')
	for cluster_id in clusters.keys():
		sig_array = clusters[cluster_id]
		rand_index = random.randint(0,len(sig_array) -1)
		t.write(str(sig_array[rand_index][1]) + '\t' + sig_array[rand_index][2] + '\t' + ref_allele + '\t' + alt_allele + '\n')
	t.close()

sqtl_file = sys.argv[1]
gencode_gene_annotation_file = sys.argv[2]
variant_position_enrichment_dir = sys.argv[3]
distance_window = sys.argv[4]
pvalue_thresh = float(sys.argv[5])

# Filter sQTL file to variant-jxn papers where variant is within a $distance_window BP around a splice site of the junction
# Also require variant to be a SNV and autosomal
tissue_name = sqtl_file.split('/')[-1].split('.')[0]
filtered_sqtl_file = variant_position_enrichment_dir + tissue_name + '_filtered_sqtls_' + distance_window + '_around_ss_snvs_only.txt'
# filter_sqtl_file_to_variants_near_ss(sqtl_file, filtered_sqtl_file, int(distance_window))

# Get mapping from ensamble_id to strand
ensamble_to_strand_mapping = get_mapping_from_ensamble_id_to_strand(gencode_gene_annotation_file)



significant_output_file = variant_position_enrichment_dir + tissue_name + '_sqtl_significant_distance_to_observed_splice_site_distance_' + distance_window + '_pvalue_thresh_' + str(pvalue_thresh) + '.txt'
extract_significant_sqtl_distances(filtered_sqtl_file, pvalue_thresh, significant_output_file, ensamble_to_strand_mapping)

background_output_file = variant_position_enrichment_dir + tissue_name + '_sqtl_background_distance_to_observed_splice_site_distance_' + distance_window + '_pvalue_thresh_' + str(pvalue_thresh) + '.txt'
extract_background_sqtl_distances(filtered_sqtl_file, pvalue_thresh, background_output_file, ensamble_to_strand_mapping)

