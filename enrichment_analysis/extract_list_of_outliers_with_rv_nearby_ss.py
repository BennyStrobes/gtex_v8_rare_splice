import numpy as np
import os
import sys
import pdb




# In each tissue, extract list of individuals that we have RNA-seq for AND Have WGS and are european ancestry
def extract_individuals_that_have_rna_and_are_european_ancestry(outlier_file, european_ancestry_individual_list):
	# Extract dictionary list of european ancestry individuals
	european_indiz = {}
	f = open(european_ancestry_individual_list)
	for line in f:
		line = line.rstrip()
		european_indiz[line] = 1
	f.close()
	# Get array of individauls for which we have RNA
	f = open(outlier_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			rna_seq_indiz = data[1:]
			continue
	f.close()
	# Extract dictionary list of the union of the above two
	union_indiz = {}
	for indi in rna_seq_indiz:
		if indi in european_indiz:
			union_indiz[indi] = 1
	return union_indiz

# Extract outliers in each tissue and save results in cluster_struct
def extract_outliers(outlier_file, individuals, pvalue_threshold, enrichment_version):
	# Initialize object keeping track of outlier calls
	cluster_struct = {}
	# Loop through outlier file
	f = open(outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if line.startswith('CLUSTER_ID'):
			#create mapping from index to individual id
			position_to_indi = {}
			for i, val in enumerate(data[1:]):
				position_to_indi[i] = val 
			continue
		# Standard Line
		cluster_id = data[0]
		# Add key (cluster_id) to cluster_struct object
		cluster_struct[cluster_id] = {}
		cluster_struct[cluster_id]['outlier_individuals'] = {}
		cluster_struct[cluster_id]['inlier_individuals'] = {}
		cluster_struct[cluster_id]['rv_individuals'] = {}
		# Don't limit to most extreme outlier / cluster. Just take everyone that passes a threshold
		if enrichment_version == 'all':
			pvalz = np.asarray(data[1:]).astype(float)
			for position, pvalue in enumerate(pvalz):
				indi = position_to_indi[position]
				if indi in individuals and pvalue < pvalue_threshold and np.isnan(pvalue) == False:
					cluster_struct[cluster_id]['outlier_individuals'][indi] = 1
				if indi in individuals and pvalue > .75 and np.isnan(pvalue) == False:
					cluster_struct[cluster_id]['inlier_individuals'][indi] = 1
		if len(cluster_struct[cluster_id]['inlier_individuals']) < 10:
			print('low')
			pdb.set_trace()
	f.close()
	return cluster_struct

def correct_alleles_for_strand(strand, major_allele, variant_allele):
	if strand == '-':
		converter = {'A':'T', 'T':'A','G':'C','C':'G'}
		return converter[major_allele], converter[variant_allele]
	else:
		return major_allele, variant_allele

# Add RV calls to cluster_struct object
def extract_rare_variants(variant_bed_file, cluster_struct, individuals, cluster_to_strand_mapping):
	# Stream variant file
	f = open(variant_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Extract relevent fields
		indi = data[3]
		cluster_id = data[10]
		var_pos = data[1]
		chrom_num = data[0]
		# Add RV in each tissue that the individual has RNA-seq
		# Don't have RNA-seq for this individual in this tissue
		if indi not in individuals:
			continue
		# Cluster_id not in this tissue
		if cluster_id not in cluster_struct:
			continue
		ref_allele = data[6]
		alt_allele = data[7]
		ref_allele_af_string = data[8]
		alt_allele_af_string = data[9]
		# Error checking
		if ref_allele_af_string.split(':')[0] != ref_allele:
			print('assumptionerror')
			pdb.set_trace()
		if alt_allele_af_string.split(':')[0] != alt_allele:
			print('assumption error')
			pdb.set_trace()
		ref_allele_af = float(ref_allele_af_string.split(':')[1])
		alt_allele_af = float(alt_allele_af_string.split(':')[1])
		if ref_allele_af <= .01:
			major_allele = alt_allele
			variant_allele = ref_allele
		elif alt_allele_af <= .01:
			major_allele = ref_allele
			variant_allele = alt_allele
		else:
			print('assumption eroror!')
		major_allele, variant_allele = correct_alleles_for_strand(cluster_to_strand_mapping[cluster_id], major_allele, variant_allele)
		cluster_struct[cluster_id]['rv_individuals'][indi] = (chrom_num, var_pos, major_allele + '->' + variant_allele)
	f.close()
	return cluster_struct

# Create mapping from cluster id to junctions
def extract_mapping_from_cluster_id_to_junctions(cluster_info_file):
	# Initialize mapping
	mapping = {}
	# Used to skip header
	head_count = 0
	# Stream cluster_info_file
	f = open(cluster_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent info
		cluster_id = data[0]
		jxn_string = data[1]
		chromosome = jxn_string.split(',')[0].split(':')[0]
		jxn_array = jxn_string.split(',')
		# Add info to mapping
		if cluster_id in mapping:
			print('assumption error')
		mapping[cluster_id] = {}
		mapping[cluster_id]['chromosome'] = chromosome
		mapping[cluster_id]['start_splice_sites'] = []
		mapping[cluster_id]['end_splice_sites'] = []
		# Add splice sites to mapping
		for jxn in jxn_array:
			jxn_info = jxn.split(':')
			start_pos = int(jxn_info[1])
			end_pos = int(jxn_info[2])
			if end_pos <= start_pos:
				print('assumptioner eroror')
				pdb.set_trace()
			mapping[cluster_id]['start_splice_sites'].append(start_pos)
			mapping[cluster_id]['end_splice_sites'].append(end_pos)
	f.close()
	return mapping

def get_cluster_to_strand_mapping(cluster_info_file, exon_file):
	gene_to_strand = {}
	cluster_to_strand = {}
	f = open(exon_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[4]
		strand = data[3]
		gene_to_strand[ensamble_id] = strand
	f.close()
	head_count = 0
	f = open(cluster_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_id = data[0]
		ensamble_ids = data[2].split(',')
		for ensamble_id in ensamble_ids:
			cluster_to_strand[cluster_id] = gene_to_strand[ensamble_id]
	f.close()
	return cluster_to_strand

# Create dictionary of annotated splice sites
# Where keys are "chr"$chrom_num"_"#splice_sites
def get_dictionary_of_annotated_splice_sites(exon_file):
	# initialize dictionary
	ss = {}
	# used to skip header
	head_count = 0
	# stream exon file
	f = open(exon_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ss[data[0] + '_' + data[1]] = 1
		ss[data[0] + '_' + data[2]] = 1
	return ss


# Compute distance from variant to nearest junction (positive values correspond to exons)
def get_distance_to_nearest_splice_site(cluster_mapping, var_pos, chromosome_string, strand, annotated_splice_sites):
	# Quick error checking
	if cluster_mapping['chromosome'] != chromosome_string:
		print('chromosome mismatch error!')
		pdb.set_trace()
	# Initialize min distance to really high number (bigger than any distance possible in genomics)
	min_distance = 10000000000000000
	ss_type = 'null'  # Corresponding to donor or acceptor
	ss_annotated = 'null'  # corresponding to annotated or novel
	# Loop through start_splice_sites
	for start_splice_site in cluster_mapping['start_splice_sites']:
		distance = start_splice_site - var_pos
		if abs(distance) <= abs(min_distance):
			min_distance = distance
			if strand == '+':
				ss_type = 'donor'
			elif strand == '-':
				ss_type = 'acceptor'
			else:
				print('strand assumption error')
				pdb.set_trace()
			if chromosome_string + '_' + str(start_splice_site) in annotated_splice_sites:
				ss_annotated = 'annotated'
			else:
				ss_annotated = 'novel'
	# Loop through end splice sites
	for end_splice_site in cluster_mapping['end_splice_sites']:
		distance = var_pos - end_splice_site
		if abs(distance) <= abs(min_distance):
			min_distance = distance
			if strand == '+':
				ss_type = 'acceptor'
			elif strand == '-':
				ss_type = 'donor'
			else:
				print('strand assumption error')
				pdb.set_trace()
			if chromosome_string + '_' + str(end_splice_site) in annotated_splice_sites:
				ss_annotated = 'annotated'
			else:
				ss_annotated = 'novel'
	if ss_type == 'null':
		print('strand assumption error!')
		pdb.set_trace()
	if ss_annotated == 'null':
		print('ss annotation assumption error')
		pdb.set_trace()
	return min_distance, ss_type, ss_annotated

def get_variant_position_around_ss_name(distance, ss_type):
	if ss_type == 'donor':
		if distance < 0:
			namer = 'D+' + str(abs(distance))
		else:
			namer = 'D-' + str(distance + 1)
	elif ss_type == 'acceptor':
		if distance < 0:
			namer = 'A-' + str(abs(distance))
		else:
			namer = 'A+' + str(distance+1)
	else:
		print('assumption error!')
		pdb.set_trace()
	return namer


def extract_readible_exon_exon_junctions(dicti):
	arr = []
	start = dicti['start_splice_sites']
	end = dicti['end_splice_sites']
	for i, ele in enumerate(start):
		arr.append(str(ele) + ':' + str(end[i]))
	return arr

#######################
# Command line args
#######################
tissue = sys.argv[1]
pvalue_threshold = float(sys.argv[2])
variant_bed_file = sys.argv[3]
clusters_to_plot_file = sys.argv[4]
splicing_outlier_dir = sys.argv[5]
splicing_outlier_suffix = sys.argv[6]
european_ancestry_individual_list = sys.argv[7]
enrichment_version = sys.argv[8]
cluster_info_file = sys.argv[9]
exon_file = sys.argv[10]






cluster_to_ss_mapping = extract_mapping_from_cluster_id_to_junctions(cluster_info_file)

cluster_to_strand_mapping = get_cluster_to_strand_mapping(cluster_info_file, exon_file)
# Create dictionary of annotated splice sites
# Where keys are "chr"$chrom_num"_"#splice_sites
annotated_splice_sites = get_dictionary_of_annotated_splice_sites(exon_file)

# In each tissue, extract list of individuals that we have RNA-seq for AND Have WGS and are european ancestry
# Extract list of individuals for this tissue
# Use outlier file to get list of individuals we have RNA-seq for
outlier_file = splicing_outlier_dir + tissue + splicing_outlier_suffix + '_emperical_pvalue.txt'
individuals = extract_individuals_that_have_rna_and_are_european_ancestry(outlier_file, european_ancestry_individual_list)


# Extract outliers in this tissue and save results in cluster_struct
outlier_file = splicing_outlier_dir + tissue + splicing_outlier_suffix + '_emperical_pvalue.txt'
cluster_struct = extract_outliers(outlier_file, individuals, pvalue_threshold, enrichment_version)

# Add RV calls to cluster_struct object
cluster_struct = extract_rare_variants(variant_bed_file, cluster_struct, individuals, cluster_to_strand_mapping)

#################
# print to output
####################
# Open output file handle and print header
t = open(clusters_to_plot_file, 'w')
t.write('cluster_id\tchrom_num\tvariant_position\toutlier_individual\tinlier_individuals\tdistance_to_ss\tstrand\tannotated_ss\tvariant_allele\n')

# Loop through cluster ids
for cluster_id in cluster_struct.keys():
	# dictionary containing all individuals that have a RV for this cluster
	rv_indi = cluster_struct[cluster_id]['rv_individuals']
	# dictionary containing all individuals that have a RV for this cluster
	outlier_indi = cluster_struct[cluster_id]['outlier_individuals']

	# Skip clusters with no outliers or no rv
	if len(outlier_indi) == 0 or len(rv_indi) == 0:
		continue

	# Loop through outliers
	for outlier_individual in outlier_indi.keys():
		# Check if outlier is also a RV
		if outlier_individual in rv_indi:
			chrom_num = rv_indi[outlier_individual][0]
			var_pos = rv_indi[outlier_individual][1]
			variant_allele = rv_indi[outlier_individual][2]

			distance, ss_type, ss_annotated = get_distance_to_nearest_splice_site(cluster_to_ss_mapping[cluster_id], int(var_pos), chrom_num, cluster_to_strand_mapping[cluster_id], annotated_splice_sites)
			variant_position_around_ss_name = get_variant_position_around_ss_name(distance, ss_type)
			strand = cluster_to_strand_mapping[cluster_id]
			inlier_string = ','.join(cluster_struct[cluster_id]['inlier_individuals'].keys())

			# Classify cluster
			#exon_exon_junctions = extract_readible_exon_exon_junctions(cluster_to_ss_mapping[cluster_id])
			#exon_skipping_bool, alternate_5_bool, alternate_3_bool = classify_cluster(exon_exon_junctions, strand)

			# Print to output file
			t.write(cluster_id + '\t' + chrom_num + '\t' + var_pos + '\t' + outlier_individual + '\t' + inlier_string + '\t' + variant_position_around_ss_name + '\t' + strand + '\t' + ss_annotated + '\t' + variant_allele + '\n')
t.close()
