import numpy as np 
import os
import sys
import pdb
import gzip


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
		pvalz = np.asarray(data[1:]).astype(float)
		# Filter out clusters that do not have any outlier individuals
		no_outliers = True
		for position, pvalue in enumerate(pvalz):
			indi = position_to_indi[position]
			if indi in individuals and pvalue < pvalue_threshold and np.isnan(pvalue) == False:
				no_outliers = False
		if no_outliers == True:
			continue
		cluster_id = data[0]
		# Add key (cluster_id) to cluster_struct object
		cluster_struct[cluster_id] = {}
		cluster_struct[cluster_id]['outlier_individuals'] = {}
		cluster_struct[cluster_id]['inlier_individuals'] = {}
		cluster_struct[cluster_id]['num_observed_samples'] = len(np.where(np.asarray(data[1:]) != 'NaN')[0])

		# Don't limit to most extreme outlier / cluster. Just take everyone that passes a threshold
		if enrichment_version == 'all':
			for position, pvalue in enumerate(pvalz):
				indi = position_to_indi[position]
				if indi in individuals and pvalue < pvalue_threshold and np.isnan(pvalue) == False:
					cluster_struct[cluster_id]['outlier_individuals'][indi] = 1
				if indi in individuals and pvalue >= pvalue_threshold and np.isnan(pvalue) == False:
					cluster_struct[cluster_id]['inlier_individuals'][indi] = 1
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

# Create mapping from cluster id to strand
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


# Print distance to splice site for outliers and non-outliers based on observed exon-exon junctions
def get_distance_from_variant_to_observed_splice_sites(cluster_struct, individuals, cluster_info_file, inlier_positions_file, outlier_positions_file, variant_bed_file, exon_file):
	# Create mapping from cluster id to junctions
	cluster_mapping = extract_mapping_from_cluster_id_to_junctions(cluster_info_file)
	# Create mapping from cluster id to strand
	cluster_to_strand_mapping = get_cluster_to_strand_mapping(cluster_info_file, exon_file)
	# Create dictionary of annotated splice sites
	# Where keys are "chr"$chrom_num"_"#splice_sites
	annotated_splice_sites = get_dictionary_of_annotated_splice_sites(exon_file)
	# Open output file handles
	t_outlier = open(outlier_positions_file, 'w')
	t_inlier = open(inlier_positions_file, 'w')
	t_outlier.write('distance\tsplice_site_type\tannotated_splice_site\tmajor_allele\tvariant_allele\tindividual\tchrom\tposition\n')
	t_inlier.write('distance\tsplice_site_type\tannotated_splice_site\tmajor_allele\tvariant_allele\tindividual\tchrom\tposition\n')
	# Stream variant bed file
	used = {}
	f = open(variant_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Parse line of variant bed file
		individual_id = data[3]
		chromosome_string = data[0]
		var_pos = int(data[1])
		cluster_id = data[10]
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
		# Ignore individuals not in individuals dictionary 
		if individual_id not in individuals:
			continue
		# Ignore variants mapped to clusters that we don't have outlier information for
		if cluster_id not in cluster_struct:
			continue
		# Individual is an outlier
		if individual_id in cluster_struct[cluster_id]['outlier_individuals']:
			# Compute distance from variant to nearest junction (positive values correspond to exons)
			distance, ss_type, ss_annotated = get_distance_to_nearest_splice_site(cluster_mapping[cluster_id], var_pos, chromosome_string, cluster_to_strand_mapping[cluster_id], annotated_splice_sites)
			t_outlier.write(str(distance) + '\t' + ss_type + '\t' + ss_annotated + '\t' + major_allele + '\t' + variant_allele + '\t' + individual_id + '\t' + chromosome_string + '\t' + str(var_pos) + '\n')
		# Individaul is an inlier
		if individual_id in cluster_struct[cluster_id]['inlier_individuals']:
			# Compute distance from variant to nearest junction (positive values correspond to exons)
			distance, ss_type, ss_annotated = get_distance_to_nearest_splice_site(cluster_mapping[cluster_id], var_pos, chromosome_string, cluster_to_strand_mapping[cluster_id], annotated_splice_sites)
			t_inlier.write(str(distance) + '\t' + ss_type + '\t' + ss_annotated + '\t' + major_allele + '\t' + variant_allele + '\t' + individual_id + '\t' + chromosome_string + '\t' + str(var_pos) + '\n')
	# CLose input and output file handles
	f.close()
	t_inlier.close()
	t_outlier.close()

# Create mapping from cluster id to genes
def extract_mapping_from_cluster_id_to_genes(cluster_info_file):
	# Initialize mapping
	mapping = {}
	# Initialize dictionary list of genes
	genes = {}
	# To skip header
	head_count = 0
	# Stream cluster info file
	f = open(cluster_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Parse line
		cluster_id = data[0]
		gene_arr = data[2].split(',')
		# Add cluster to mapping
		mapping[cluster_id] = gene_arr
		# Add gene to dictionary list
		for gene in gene_arr:
			genes[gene] = []
	return mapping, genes

def get_exon_positions_of_genes(gencode_gene_annotation_file, genes):
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
		if hits != 1:
			print('fundamental assumption error!!')
			pdb.set_trace()
		line_chrom_num = data[0]
		gene_part = data[2]  # gene,UTR,exon,etc
		# Skip genes not in our valid gene set
		if gene_name not in genes:
			continue
		if gene_part != 'exon':  # ignore other parts of the gene as 'gene' encomposes everything
			continue
		start = int(data[3])  # 5' gene start site (this is the min of all UTR,exons,etc)
		end = int(data[4])  # 3' gene end site (this is the max of all UTR,exons,etc)
		# Simple error checking
		if end < start:
			print('assumption error')
			pdb.set_trace()
		# Add exon to gene_mapping
		genes[gene_name].append((start,end))
	f.close()
	return genes

# Compute distance from variant to nearest junction (positive values correspond to exons)
def get_distance_to_nearest_exon_splice_site(gene_arr, genes, var_pos):
	# Initialize min distance to really high number (bigger than any distance possible in genomics)
	min_distance = 10000000000000000
	# Loop through genes
	for gene in gene_arr:
		# Extract list of exons for this gene
		exon_list = genes[gene]
		# Loop through exons
		for exon in exon_list:
			exon_start = exon[0]
			exon_end = exon[1]
			# Compute min distance from variant to splice site on this exon
			distance = min(abs(exon_start - var_pos), abs(exon_end - var_pos))
			# If not in exon
			if var_pos > exon_end or var_pos < exon_start:
				distance = distance*-1
			if abs(distance) <= abs(min_distance):
				min_distance = distance
	return distance

def check_if_variant_in_annotated_exon(var_pos, gene_arr, genes):
	booly = False
	for gene in gene_arr:
		exon_list = genes[gene]
		for exon in exon_list:
			exon_start = exon[0]
			exon_end = exon[1]
			if var_pos >= exon_start and var_pos <= exon_end:
				booly = True 
	return booly






#########################
# Command line args
#########################
rare_variant_dir = sys.argv[1]  # Contains rare variants
variant_position_enrichment_dir = sys.argv[2]  # Output dir
splicing_outlier_dir = sys.argv[3]  # Directory containing outlier calls
splicing_outlier_suffix = sys.argv[4]  # Suffix to files in splicing outlier dir
european_ancestry_individual_list = sys.argv[5]  # File containing list of all european ancestry individuals
gencode_gene_annotation_file = sys.argv[6]  # Gencode gene annotation file
cluster_info_file = sys.argv[7]  # File containing mapping from clusters to both genes and the composing junctions
pvalue_threshold = float(sys.argv[8])  # threshold for outlier calling
distance = sys.argv[9]  # Window size around splice sites that we consider rare variants for enrichment
exon_file = sys.argv[10]


# Outlier file for multi-tissue outliers
splicing_outlier_file = splicing_outlier_dir + 'cross_tissue' + splicing_outlier_suffix + '_emperical_pvalue.txt'
# Bed file for rare variants
variant_bed_file = rare_variant_dir + 'variant_cluster_only_bed_' + distance + '.txt'


# Use outlier file to get list of individuals we have RNA-seq for
individuals = extract_individuals_that_have_rna_and_are_european_ancestry(splicing_outlier_file, european_ancestry_individual_list)

# Extract object containing which individuals are outliers for which cluster
cluster_struct = extract_outliers(splicing_outlier_file, individuals, pvalue_threshold, 'all')

# Print distance to splice site for outliers and non-outliers based on observed exon-exon junctions
outlier_positions_file = variant_position_enrichment_dir + 'outlier_distance_to_observed_splice_site_distance_' + distance + '_pvalue_thresh_' + str(pvalue_threshold) + '.txt'  # For outliers
inlier_positions_file = variant_position_enrichment_dir + 'inlier_distance_to_observed_splice_site_distance_' + distance + '_pvalue_thresh_' + str(pvalue_threshold) + '.txt'  # For inliers
get_distance_from_variant_to_observed_splice_sites(cluster_struct, individuals, cluster_info_file, inlier_positions_file, outlier_positions_file, variant_bed_file, exon_file)

