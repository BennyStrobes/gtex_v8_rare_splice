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



# Add RV calls to cluster_struct object
def extract_rare_variants(variant_bed_file, cluster_struct, individuals):
	# Stream variant file
	f = open(variant_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Extract relevent fields
		indi = data[0]
		cluster_id = data[6]
		var_pos = data[2]
		# Add RV in each tissue that the individual has RNA-seq
		# Don't have RNA-seq for this individual in this tissue
		if indi not in individuals:
			continue
		# Cluster_id not in this tissue
		if cluster_id not in cluster_struct:
			continue
		cluster_struct[cluster_id]['rv_individuals'][indi] = var_pos
	f.close()
	return cluster_struct

def get_cluster_to_ss_mapping(cluster_info_file):
	mapping = {}
	head_count = 0
	f = open(cluster_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_id = data[0]
		jxns = data[1].split(',')
		start_sites = []
		end_sites = []
		for jxn in jxns:
			start_site = int(jxn.split(':')[1])
			end_site = int(jxn.split(':')[2])
			start_sites.append(start_site)
			end_sites.append(end_site)
		mapping[cluster_id] = {}
		mapping[cluster_id]['start_sites'] = np.asarray(start_sites)
		mapping[cluster_id]['end_sites'] = np.asarray(end_sites)
	f.close()
	return mapping

def get_position_relative_to_ss(var_pos, ss):
	start_sites = ss['start_sites']
	end_sites = ss['end_sites']
	min_distance = 1000000000000

	for ss in start_sites:
		disty = ss - var_pos
		if abs(disty) < abs(min_distance):
			min_distance = disty
	for ss in end_sites:
		disty = var_pos - ss
		if abs(disty) < abs(min_distance):
			min_distance = disty
	return str(min_distance)

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






cluster_to_ss_mapping = get_cluster_to_ss_mapping(cluster_info_file)

cluster_to_strand_mapping = get_cluster_to_strand_mapping(cluster_info_file, exon_file)

# In each tissue, extract list of individuals that we have RNA-seq for AND Have WGS and are european ancestry
# Extract list of individuals for this tissue
# Use outlier file to get list of individuals we have RNA-seq for
outlier_file = splicing_outlier_dir + tissue + splicing_outlier_suffix + '_emperical_pvalue.txt'
individuals = extract_individuals_that_have_rna_and_are_european_ancestry(outlier_file, european_ancestry_individual_list)


# Extract outliers in this tissue and save results in cluster_struct
outlier_file = splicing_outlier_dir + tissue + splicing_outlier_suffix + '_emperical_pvalue.txt'
cluster_struct = extract_outliers(outlier_file, individuals, pvalue_threshold, enrichment_version)

# Add RV calls to cluster_struct object
cluster_struct = extract_rare_variants(variant_bed_file, cluster_struct, individuals)

#################
# print to output
####################
# Open output file handle and print header
t = open(clusters_to_plot_file, 'w')
t.write('cluster_id\tvariant_position\toutlier_individual\tinlier_individuals\tdistance_to_ss\tstrand\n')

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
			var_pos = rv_indi[outlier_individual]
			position_relative_to_ss = get_position_relative_to_ss(int(var_pos), cluster_to_ss_mapping[cluster_id])
			strand = cluster_to_strand_mapping[cluster_id]
			inlier_string = ','.join(cluster_struct[cluster_id]['inlier_individuals'].keys())
			# Print to output file
			t.write(cluster_id + '\t' + var_pos + '\t' + outlier_individual + '\t' + inlier_string + '\t' + position_relative_to_ss + '\t' + strand + '\n')
t.close()
