import numpy as np 
import os
import sys
import pdb







# Extract vector of tissue names
def get_tissue_array(tissue_names_file):
	arr = []
	f = open(tissue_names_file)
	for line in f:
		line = line.rstrip()
		arr.append(line)
	f.close()
	return arr



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
		cluster_struct[cluster_id]['rv_individuals'] = {}
		# Don't limit to most extreme outlier / cluster. Just take everyone that passes a threshold
		if enrichment_version == 'all':
			pvalz = np.asarray(data[1:]).astype(float)
			for position, pvalue in enumerate(pvalz):
				indi = position_to_indi[position]
				##$$$$$##
				# if indi in individuals and pvalue < pvalue_threshold and np.isnan(pvalue) == False:
				##$$$$$##
				if indi in individuals and pvalue > pvalue_threshold and np.isnan(pvalue) == False:
					cluster_struct[cluster_id]['outlier_individuals'][indi] = 1
	f.close()
	return cluster_struct

# Add RV calls to cluster_struct object
def extract_rare_variants(variant_bed_file, cluster_struct, individuals, tissues):
	# Stream variant file
	f = open(variant_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Extract relevent fields
		indi = data[0]
		cluster_id = data[6]
		# Add RV in each tissue that the individual has RNA-seq
		for tissue in tissues:
			# Don't have RNA-seq for this individual in this tissue
			if indi not in individuals[tissue]:
				continue
			# Cluster_id not in this tissue
			if cluster_id not in cluster_struct[tissue]:
				continue
			cluster_struct[tissue][cluster_id]['rv_individuals'][indi] = 1
	f.close()
	return cluster_struct

# Do enrichment analysis in each tissue
def enrichment_analysis(tissue, cluster_struct, output_handle, num_samp):
	# Keep track of the following quantities
	a = 0 #outliers with rare variant
	b = 0 #outliers
	c = 0 #non outliers with rare variant
	d = 0 #non outliers

	# Loop through clusters
	for cluster_id in sorted(cluster_struct.keys()):
		# dictionary containing all individuals that have a RV for this cluster
		rv_indi = cluster_struct[cluster_id]['rv_individuals']
		# dictionary containing all individuals that have a RV for this cluster
		outlier_indi = cluster_struct[cluster_id]['outlier_individuals']

		# Skip Clusters that have no outliers
		if len(outlier_indi) == 0:
			continue
		count = 0 #keep track of number of outlier individuals we have for this gene only
		count_rv = 0 #keep track of number of individuals that are outliers and have rv for this gene only
		for indi_o in outlier_indi.keys():
			count = count + 1
			if indi_o in rv_indi:
				count_rv = count_rv + 1
				a = a + 1
			b = b + 1
		num_non_outliers = num_samp - count
		temp_c = len(rv_indi) - count_rv
		temp_d =num_non_outliers
		c = c + temp_c
		d = d + temp_d
	# Compute odds ratios
	num_frac = float(a)/float(b)
	den_frac = float(c)/float(d)
	odds_ratio = (num_frac/den_frac)
	output_handle.write(tissue + '\t' + str(a) + '\t' + str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(odds_ratio) + '\n')

	return output_handle



#############################
# Command Line Args
#############################
variant_bed_file = sys.argv[1]
output_root = sys.argv[2]
splicing_outlier_dir = sys.argv[3]
splicing_outlier_suffix = sys.argv[4]
european_ancestry_individual_list = sys.argv[5]
tissue_names_file = sys.argv[6]
enrichment_version = sys.argv[7]
pvalue_threshold = float(sys.argv[8])

# Extract vector of tissue names
tissues = get_tissue_array(tissue_names_file)

##$$$$$##
tissues = tissues[33:34]
##$$$$$##

# In each tissue, extract list of individuals that we have RNA-seq for AND Have WGS and are european ancestry
individuals = {}
# Extract list of individuals for each tissue
for tissue in tissues:
	# Use outlier file to get list of individuals we have RNA-seq for
	##$$$$$##
	outlier_file = splicing_outlier_dir + tissue + splicing_outlier_suffix + '_md.txt'
	#outlier_file = splicing_outlier_dir + tissue + splicing_outlier_suffix + '_emperical_pvalue.txt'
	##$$$$$##
	individuals[tissue] = extract_individuals_that_have_rna_and_are_european_ancestry(outlier_file, european_ancestry_individual_list)
	

# Initialize object to keep track of outliers and rvs
cluster_struct = {}
# Extract outliers in each tissue and save results in cluster_struct
for tissue in tissues:
	##$$$$$##
	outlier_file = splicing_outlier_dir + tissue + splicing_outlier_suffix + '_md.txt'
	#outlier_file = splicing_outlier_dir + tissue + splicing_outlier_suffix + '_emperical_pvalue.txt'
	##$$$$$##
	cluster_struct[tissue] = extract_outliers(outlier_file, individuals[tissue], pvalue_threshold, enrichment_version)

# Add RV calls to cluster_struct object
cluster_struct = extract_rare_variants(variant_bed_file, cluster_struct, individuals, tissues)

# Output File
output_file = output_root + '.txt'
output_handle = open(output_file, 'w')
# Do enrichment analysis in each tissue
for tissue in tissues:
	output_handle = enrichment_analysis(tissue, cluster_struct[tissue], output_handle, len(individuals[tissue]))


output_handle.close()

