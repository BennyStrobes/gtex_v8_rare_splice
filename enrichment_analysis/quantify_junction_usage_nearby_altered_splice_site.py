import numpy as np 
import os
import sys
import pdb

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

def get_tissue_names(tissue_names_file):
	tissue_names = []
	f = open(tissue_names_file)
	for line in f:
		line = line.rstrip()
		tissue_names.append(line)
	f.close()
	return tissue_names

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
def extract_outliers(outlier_file, individuals, pvalue_outlier_threshold, pvalue_inlier_threshold, enrichment_version):
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
				if indi in individuals and pvalue < pvalue_outlier_threshold and np.isnan(pvalue) == False:
					cluster_struct[cluster_id]['outlier_individuals'][indi] = 1
				if indi in individuals and pvalue > pvalue_inlier_threshold and np.isnan(pvalue) == False:
					cluster_struct[cluster_id]['inlier_individuals'][indi] = 1
		if len(cluster_struct[cluster_id]['inlier_individuals']) < 1:
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
	ss_name = 'null'  # Corresponding to chromNum:Pos
	# Loop through start_splice_sites
	for start_splice_site in cluster_mapping['start_splice_sites']:
		distance = start_splice_site - var_pos
		if abs(distance) <= abs(min_distance):
			min_distance = distance
			ss_name = chromosome_string + ':' + str(start_splice_site)
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
			ss_name = chromosome_string + ':' + str(end_splice_site)
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
	return min_distance, ss_type, ss_annotated, ss_name


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

def get_to_or_from_concensus(concensus_allele, variant_allele):
	old = variant_allele.split('->')[0]
	new = variant_allele.split('->')[1]
	if old == concensus_allele:
		return 'from_concensus'
	elif new == concensus_allele:
		return 'to_concensus'
	else:
		return 'neither'


# Part 1: Extract file containing cases where outlier (individual, cluster) has rare variant in nearby concensus site around splice site for that cluster
# Concensus sites include A-2, A-1, A+1, D-1, D+1, D+2,D+3,D+4,D+5,D+6
def extract_cases_where_outlier_individuals_has_concensus_variant(tissue_names, cluster_to_ss_mapping, cluster_to_strand_mapping, annotated_splice_sites, pvalue_outlier_threshold, pvalue_inlier_threshold, variant_bed_file, splicing_outlier_dir, splicing_outlier_suffix, european_ancestry_individual_list, output_file):
	valid_positions = {'A-2':'A','A-1':'G', 'A+1':'G', 'D-1':'G', 'D+1':'G', 'D+2':'T', 'D+3':'A', 'D+4':'A', 'D+5':'G', 'D+6':'T'}
	# Open output file handle
	t = open(output_file, 'w')
	# Print header
	t.write('tissue\tcluster_id\tchrom_num\tvariant_position\tss_name\toutlier_individual\tinlier_individuals\tdistance_to_ss\tstrand\tannotated_ss\tvariant_allele\tto_or_from_concensus\n')

	# Loop through tissues
	for tissue in tissue_names:
		print(tissue)
		# In each tissue, extract list of individuals that we have RNA-seq for AND Have WGS and are european ancestry
		# Extract list of individuals for this tissue
		# Use outlier file to get list of individuals we have RNA-seq for
		outlier_file = splicing_outlier_dir + tissue + splicing_outlier_suffix + '_merged_emperical_pvalue.txt'
		individuals = extract_individuals_that_have_rna_and_are_european_ancestry(outlier_file, european_ancestry_individual_list)

		# Extract outliers in this tissue and save results in cluster_struct
		outlier_file = splicing_outlier_dir + tissue + splicing_outlier_suffix + '_merged_emperical_pvalue.txt'
		cluster_struct = extract_outliers(outlier_file, individuals, pvalue_outlier_threshold, pvalue_inlier_threshold, "all")

		# Add RV calls to cluster_struct object
		cluster_struct = extract_rare_variants(variant_bed_file, cluster_struct, individuals, cluster_to_strand_mapping)
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

					distance, ss_type, ss_annotated, ss_name = get_distance_to_nearest_splice_site(cluster_to_ss_mapping[cluster_id], int(var_pos), chrom_num, cluster_to_strand_mapping[cluster_id], annotated_splice_sites)
					variant_position_around_ss_name = get_variant_position_around_ss_name(distance, ss_type)
					if variant_position_around_ss_name not in valid_positions:
						continue
					strand = cluster_to_strand_mapping[cluster_id]
					
					to_or_from_concensus = get_to_or_from_concensus(valid_positions[variant_position_around_ss_name], variant_allele)
					inlier_string = ','.join(cluster_struct[cluster_id]['inlier_individuals'].keys())
					# Print to output file
					t.write(tissue + '\t' + cluster_id + '\t' + chrom_num + '\t' + var_pos + '\t' + ss_name + '\t' + outlier_individual + '\t' + inlier_string + '\t' + variant_position_around_ss_name + '\t' + strand + '\t' + ss_annotated + '\t' + variant_allele + '\t' + to_or_from_concensus + '\n')
	t.close()

# Extract object containing junction counts
def extract_junction_counts_object(jxn_file):
	jxn_object = {}
	f = open(jxn_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			indi = np.asarray(data[1:])
			continue
		jxn_name = data[0]
		counts = np.asarray(data[1:]).astype(float)
		jxn_id = jxn_name.split(':')[0] + ':' + jxn_name.split(':')[1] + ':' + jxn_name.split(':')[2]
		cluster_id = jxn_name.split(':')[3]
		if cluster_id not in jxn_object:
			jxn_object[cluster_id] = {}
		# Error checking
		if jxn_id in jxn_object[cluster_id]:
			print('assumption error')
			pdb.set_trace()
		jxn_object[cluster_id][jxn_id] = counts
	f.close()
	return jxn_object, indi

def get_dictionary_from_comma_seperated_string(stringy):
	arr = stringy.split(',')
	dicti = {}
	for ele in arr:
		dicti[ele] = 1
	return dicti

def ss_in_rna_seq_data(ss_position, cluster_counts_object):
	ss_in_data = False
	for jxn in cluster_counts_object.keys():
		jxn_start = jxn.split(':')[1]
		jxn_end = jxn.split(':')[2]
		if jxn_start == ss_position:
			ss_in_data = True
		if jxn_end == ss_position:
			ss_in_data = True
	return ss_in_data

def get_read_counts(ss_position, individuals, cluster_counts_object, indi_arr):
	ss_count = 0
	cluster_count = 0
	for jxn in cluster_counts_object.keys():
			desired_jxn = False
			jxn_start = jxn.split(':')[1]
			jxn_end = jxn.split(':')[2]
			if jxn_end == ss_position or jxn_start == ss_position:
				desired_jxn = True
			count_arr = cluster_counts_object[jxn]
			for index, count in enumerate(count_arr):
				index_individual = indi_arr[index]
				if index_individual in individuals:
					cluster_count = cluster_count + count
					if desired_jxn == True:
						ss_count = ss_count + count
	return ss_count, cluster_count

# Part 2: Add read count contingency table for each variant-jxn pair
def add_read_count_contingency_table(input_file, output_file, tissue_names, filtered_cluster_dir):
	t = open(output_file, 'w')
	# Loop through tissues
	for i, tissue_name in enumerate(tissue_names):
		print(tissue_name)
		# Extract object containing junction counts
		jxn_file = filtered_cluster_dir + tissue_name + '_filtered_jxns_cross_tissue_clusters_gene_mapped.txt'
		junction_counts_object, indi_arr = extract_junction_counts_object(jxn_file)
		# Stream input file, but limit to cases that are in this tissue
		f = open(input_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				if i == 0:
					t.write(line + '\toutlier_ss_counts\toutlier_cluster_counts\tinlier_ss_counts\tinlier_cluster_counts\todds_ratio\n')
				continue
			# Extract relevent fields
			line_tissue = data[0]
			cluster_name = data[1]
			ss_name = data[4]
			ss_position = ss_name.split(':')[1]
			inlier_individuals = get_dictionary_from_comma_seperated_string(data[6])
			outlier_individuals = {data[5]: 1}
			if line_tissue != tissue_name:
				continue
			if ss_in_rna_seq_data(ss_position, junction_counts_object[cluster_name]) == False:
				continue
			outlier_ss_counts, outlier_cluster_counts = get_read_counts(ss_position, outlier_individuals, junction_counts_object[cluster_name], indi_arr)
			inlier_ss_counts, inlier_cluster_counts = get_read_counts(ss_position, inlier_individuals, junction_counts_object[cluster_name], indi_arr)
			orat = (float(outlier_ss_counts + 1)/float(outlier_cluster_counts+1))/(float(inlier_ss_counts+1)/float(inlier_cluster_counts+1))
			t.write(line + '\t' + str(outlier_ss_counts) + '\t' + str(outlier_cluster_counts) + '\t' + str(inlier_ss_counts) + '\t' + str(inlier_cluster_counts) + '\t' + str(orat) + '\n')
		f.close()

	t.close()



pvalue_outlier_threshold = float(sys.argv[1])
pvalue_inlier_threshold = float(sys.argv[2])
variant_bed_file = sys.argv[3]
splicing_outlier_dir = sys.argv[4]
splicing_outlier_suffix = sys.argv[5]
european_ancestry_individual_list = sys.argv[6]
gencode_gene_annotation_file = sys.argv[7]
cluster_info_file = sys.argv[8]
exon_file = sys.argv[9]
output_dir = sys.argv[10]
tissue_names_file = sys.argv[11]
filtered_cluster_dir = sys.argv[12]


# Extract dictionaries containing splice site information
cluster_to_ss_mapping = extract_mapping_from_cluster_id_to_junctions(cluster_info_file)
cluster_to_strand_mapping = get_cluster_to_strand_mapping(cluster_info_file, exon_file)
annotated_splice_sites = get_dictionary_of_annotated_splice_sites(exon_file)

# Extract array of tissue names
tissue_names = get_tissue_names(tissue_names_file)


# Part 1: Extract file containing cases where outlier (individual, cluster) has rare variant in nearby concensus site around splice site for that cluster
# Concensus sites include A-2, A-1, A+1, D-1, D+1, D+2,D+3,D+4,D+5,D+6
output_file = output_dir + 'tissue_by_tissue_outliers_with_rv_in_concensus_sites_outlier_individuals_' + str(pvalue_outlier_threshold) + '_inlier_individuals_' + str(pvalue_inlier_threshold) + '.txt'
extract_cases_where_outlier_individuals_has_concensus_variant(tissue_names, cluster_to_ss_mapping, cluster_to_strand_mapping, annotated_splice_sites, pvalue_outlier_threshold, pvalue_inlier_threshold, variant_bed_file, splicing_outlier_dir, splicing_outlier_suffix, european_ancestry_individual_list, output_file)


# Part 2: Add read count contingency table for each variant-jxn pair
output_file2 = output_dir + 'tissue_by_tissue_outliers_with_rv_in_concensus_sites_outlier_individuals_' + str(pvalue_outlier_threshold) + '_inlier_individuals_' + str(pvalue_inlier_threshold) + '_with_read_counts.txt'
add_read_count_contingency_table(output_file, output_file2, tissue_names, filtered_cluster_dir)




