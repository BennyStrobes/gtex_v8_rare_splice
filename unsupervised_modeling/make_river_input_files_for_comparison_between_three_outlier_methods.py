import numpy as np 
import os
import sys
import pdb


def get_total_expression_individuals(total_expression_outlier_file):
	indi = {}
	head_count = 0
	f = open(total_expression_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		indi[data[0]] = 1
	f.close()
	return indi

def get_ase_individuals(ase_outlier_file):
	indi = {}
	f = open(ase_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		indi[data[0]] = 1
	f.close()
	return indi

def get_splicing_individuals(splicing_outlier_file):
	indi = {}
	head_count = 0
	f = open(splicing_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for ele in data[1:]:
				indi[ele] = 1
	f.close()
	return indi

# Extract dictionary list of individuals used in all three methods
def get_list_of_individuals_used_in_all_methods(total_expression_outlier_file, ase_outlier_file, splicing_outlier_file):
	total_expression_individuals = get_total_expression_individuals(total_expression_outlier_file)
	ase_individuals = get_ase_individuals(ase_outlier_file)
	splicing_individuals = get_splicing_individuals(splicing_outlier_file)
	# Take union of the above three dictionaries
	union_individuals = {}
	for indi in splicing_individuals.keys():
		if indi in ase_individuals and indi in total_expression_individuals:
			union_individuals[indi] = 1
	return union_individuals

def get_total_expression_genes(total_expression_outlier_file):
	genes = {}
	head_count = 0
	f = open(total_expression_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes[data[1].split('.')[0]] = 1
	f.close()
	return genes

def get_ase_genes(ase_outlier_file):
	genes = {}
	f = open(ase_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		genes[data[1]] = 1
	f.close()
	return genes

def get_splicing_genes(splicing_outlier_file):
	f = open(splicing_outlier_file)
	head_count = 0
	genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble = data[0].split('.')[0]
		genes[ensamble] = 1
	f.close()
	return genes

# Extract dictionary list of genes used in all three methods
def get_list_of_genes_used_in_all_methods(total_expression_outlier_file, ase_outlier_file, splicing_outlier_file):
	splicing_genes = get_splicing_genes(splicing_outlier_file)
	#ase_genes = get_ase_genes(ase_outlier_file)
	total_expression_genes = get_total_expression_genes(total_expression_outlier_file)
	union_genes = {}
	for gene in splicing_genes.keys():
		#if gene in ase_genes and gene in total_expression_genes:
		if gene in total_expression_genes:
			union_genes[gene] = 1
	return union_genes

# For splicing outliers extract list of genes for which we have at least one outlier sample and extract pvalue threshold that gives us num_outlier_samples
def get_splicing_outlier_info(individuals, genes, num_outlier_samples, splicing_outlier_file):
	# Take first pass to extract array of pvalues (spanning all gene-individual pairs)
	pvalue_array = []
	head_count = 0
	f = open(splicing_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if head_count == 0:
			head_count = head_count + 1
			indiz = np.asarray(data[1:])
			# Get array of columns that map to individuals we are using
			good_columns = []
			for i,indi in enumerate(indiz):
				if indi in individuals:
					good_columns.append(i)
			# ERROR CHECKING
			new_indiz = indiz[good_columns]
			for indi in new_indiz:
				if indi not in individuals:
					print('assumption error')
					pdb.set_trace()
			continue
		# Skip genes not in genes dicti
		ensamble = data[0].split('.')[0]
		if ensamble not in genes:
			continue
		# get pvalue vector
		pvalues = np.asarray(data[1:]).astype(float)
		# Remove columns not in individuals
		pvalues_filtered = pvalues[good_columns]
		# Add to global pvalue_array
		for pvalue in pvalues_filtered:
			if np.isnan(pvalue):
				continue
			pvalue_array.append(pvalue)
	f.close()
	# Sort pvalue array
	sorted_pvalue_array = sorted(pvalue_array)
	# Get pvalue threshold
	if sorted_pvalue_array[num_outlier_samples-1] == sorted_pvalue_array[num_outlier_samples]:
		print('assumption error')
		pdb.set_trace()
	pvalue_threshold = (sorted_pvalue_array[num_outlier_samples-1] + sorted_pvalue_array[num_outlier_samples])/2.0
	# Take second pass through getting list of splicing outlier genes
	splicing_outlier_genes = {}
	head_count = 0
	f = open(splicing_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Skip genes not in genes dicti
		ensamble = data[0].split('.')[0]
		if ensamble not in genes:
			continue
		# get pvalue vector
		pvalues = np.asarray(data[1:]).astype(float)
		# Remove columns not in individuals
		pvalues_filtered = pvalues[good_columns]
		for pvalue in pvalues_filtered:
			if np.isnan(pvalue):
				continue
			# Check if (individual, gene) is an outlier. If so, add gene to dictionary
			if pvalue < pvalue_threshold:
				splicing_outlier_genes[ensamble] = 1
	f.close()
	return splicing_outlier_genes, pvalue_threshold


# For total expression outliers extract list of genes for which we have at least one outlier sample and extract pvalue threshold that gives us num_outlier_samples
def get_total_expression_outlier_info(individuals, genes, num_outlier_samples, total_expression_outlier_file):
	# Take first pass to extract array of zscores (spanning all gene-individual pairs)
	zscore_array = []
	head_count = 0
	f = open(total_expression_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		indi_id = data[0]
		ensamble_id = data[1].split('.')[0]
		if indi_id not in individuals:
			continue
		if ensamble_id not in genes:
			continue
		zscore = abs(float(data[4]))
		zscore_array.append(zscore)
	f.close()
	sorted_zscore_array = sorted(zscore_array,reverse=True)
	if sorted_zscore_array[num_outlier_samples-1] == sorted_zscore_array[num_outlier_samples]:
		print('assumption error')
		pdb.set_trace()
	zscore_threshold = (sorted_zscore_array[num_outlier_samples-1] + sorted_zscore_array[num_outlier_samples])/2.0
	total_expression_outlier_genes = {}
	f = open(total_expression_outlier_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		indi_id = data[0]
		ensamble_id = data[1].split('.')[0]
		if indi_id not in individuals:
			continue
		if ensamble_id not in genes:
			continue
		zscore = abs(float(data[4]))
		if zscore > zscore_threshold:
			total_expression_outlier_genes[ensamble_id] = 1
	f.close()
	return total_expression_outlier_genes, zscore_threshold

def get_ase_outlier_info(individuals, genes, num_outlier_samples, ase_outlier_file):
	pvalue_array = []
	f = open(ase_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		indi_id = data[0]
		ensamble_id = data[1]
		if indi_id not in individuals:
			continue
		if ensamble_id not in genes:
			continue
		pvalue_array.append(float(data[2]))
	f.close()
	# Sort pvalue array
	sorted_pvalue_array = sorted(pvalue_array)
	# Get pvalue threshold
	if sorted_pvalue_array[num_outlier_samples-1] == sorted_pvalue_array[num_outlier_samples]:
		print('assumption error')
		pdb.set_trace()
	pvalue_threshold = (sorted_pvalue_array[num_outlier_samples-1] + sorted_pvalue_array[num_outlier_samples])/2.0
	ase_outlier_genes = {}
	f = open(ase_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		indi_id = data[0]
		ensamble_id = data[1]
		if indi_id not in individuals:
			continue
		if ensamble_id not in genes:
			continue
		pvalue = float(data[2])
		if pvalue < pvalue_threshold:
			ase_outlier_genes[ensamble_id] = 1
	f.close()
	return ase_outlier_genes, pvalue_threshold

def union_of_three_dictionaries(dicti1, dicti2, dicti3):
	union = {}
	for key  in dicti1.keys():
		union[key] = 1
	for key  in dicti2.keys():
		union[key] = 1	
	for key  in dicti3.keys():
		union[key] = 1
	return union

def get_splicing_outliers(genes, individuals, splicing_pvalue_threshold, splicing_outlier_file):
	outlier_dicti = {}
	head_count = 0
	f = open(splicing_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if head_count == 0:
			head_count = head_count + 1
			indiz = np.asarray(data[1:])
			# Get array of columns that map to individuals we are using
			good_columns = []
			for i,indi in enumerate(indiz):
				if indi in individuals:
					good_columns.append(i)
			# ERROR CHECKING
			new_indiz = indiz[good_columns]
			for indi in new_indiz:
				if indi not in individuals:
					print('assumption error')
					pdb.set_trace()
			continue
		# Skip genes not in genes dicti
		ensamble = data[0].split('.')[0]
		if ensamble not in genes:
			continue
		# get pvalue vector
		pvalues = np.asarray(data[1:]).astype(float)
		# Remove columns not in individuals
		pvalues_filtered = pvalues[good_columns]
		# Add to global pvalue_array
		for i,pvalue in enumerate(pvalues_filtered):
			indi_id = new_indiz[i]
			sample_name = indi_id + '_' + ensamble
			if pvalue < splicing_pvalue_threshold:
				outlier_dicti[sample_name] = 1
			else:
				outlier_dicti[sample_name] = 0
	f.close()
	return outlier_dicti

def get_total_expression_outliers(genes, individuals, total_expression_zscore_threshold, total_expression_outlier_file):
	outlier_dicti = {}
	head_count = 0
	f = open(total_expression_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		indi_id = data[0]
		ensamble_id = data[1].split('.')[0]
		if indi_id not in individuals:
			continue
		if ensamble_id not in genes:
			continue
		zscore = abs(float(data[4]))
		sample_name = indi_id + '_' + ensamble_id
		if zscore > total_expression_zscore_threshold:
			outlier_dicti[sample_name] = 1
	f.close()
	for indi in individuals.keys():
		for gene in genes.keys():
			sample_name = indi + '_' + gene
			if sample_name not in outlier_dicti:
				outlier_dicti[sample_name] = 0
	return outlier_dicti

def get_ase_outliers(genes, individuals, ase_pvalue_threshold, ase_outlier_file):
	outlier_dicti = {}
	f = open(ase_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		indi_id = data[0]
		ensamble_id = data[1]
		if indi_id not in individuals:
			continue
		if ensamble_id not in genes:
			continue
		pvalue = float(data[2])
		sample_name = indi_id + '_' + ensamble_id
		if pvalue <  ase_pvalue_threshold:
			outlier_dicti[sample_name] = 1
	f.close()
	for indi in individuals.keys():
		for gene in genes.keys():
			sample_name = indi + '_' + gene
			if sample_name not in outlier_dicti:
				outlier_dicti[sample_name] = 0
	return outlier_dicti


def print_outlier_output_file(outlier_dicti, genomic_annotation_file, output_file):
	t = open(output_file, 'w')
	f = open(genomic_annotation_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'Zscore\n')
			continue
		indi_id = data[0]
		gene_id = data[1].split('.')[0]
		sample_name = indi_id + '_' + gene_id
		if sample_name in outlier_dicti:
			t.write(line + '\t' + str(outlier_dicti[sample_name]) + '\n')
	f.close()
	t.close()
	remove_no_variance_features(output_file)
	add_n2_pairs_column(output_file.split('.')[0] + '_features_filter.txt')


def add_n2_pairs_column(input_file):
	dicti = {}
	output_file = input_file.split('.')[0] + '_N2_pairs.txt'
	f = open(input_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header = line
			continue
		liner_short = '\t'.join(data[1:-1])
		liner_long = '\t'.join(data)
		if liner_short not in dicti:
			dicti[liner_short] = []
		dicti[liner_short].append(liner_long)
	f.close()
	t = open(output_file,'w')
	t.write(header + '\tN2pair\n')
	for key in dicti.keys():
		arr = dicti[key]
		if len(arr) == 1:
			liner_long = arr[0]
			t.write(liner_long + '\t' + 'NA' + '\n')
	n2_pair_counter = 0
	for key in dicti.keys():
		arr = dicti[key]
		if len(arr) > 1:
			n2_pair_counter = n2_pair_counter + 1
			for liner_long in arr[:2]:
				t.write(liner_long + '\t' + str(n2_pair_counter) + '\n')
	t.close()

def remove_no_variance_features(input_file):
	output_file = input_file.split('.')[0] + '_features_filter.txt'
	f = open(input_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for i, ele in enumerate(data):
				arr.append({})
			continue
		for i,ele in enumerate(data):
			arr[i][ele] = 1
	f.close()
	good_columns = []
	for i, dicti in enumerate(arr):
		if len(dicti) > 1:
			good_columns.append(i)
	t = open(output_file, 'w')
	head_count = 0
	f = open(input_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		new_data = np.asarray(data)[good_columns]
		t.write('\t'.join(new_data) + '\n')
	f.close()
	t.close()



genomic_annotation_file = sys.argv[1]
total_expression_outlier_file = sys.argv[2]
ase_outlier_file = sys.argv[3]
splicing_outlier_file = sys.argv[4]
unsupervised_learning_input_dir = sys.argv[5]
num_outlier_samples = int(sys.argv[6])

# Extract dictionary list of individuals used in all three methods
individuals = get_list_of_individuals_used_in_all_methods(total_expression_outlier_file, ase_outlier_file, splicing_outlier_file)

# Extract dictionary list of genes used in all three methods
genes = get_list_of_genes_used_in_all_methods(total_expression_outlier_file, ase_outlier_file, splicing_outlier_file)

# For splicing outliers extract list of genes for which we have at least one outlier sample and extract pvalue threshold that gives us num_outlier_samples
splicing_outlier_genes, splicing_pvalue_threshold = get_splicing_outlier_info(individuals, genes, num_outlier_samples, splicing_outlier_file)

# For total expression outliers extract list of genes for which we have at least one outlier sample and extract zscore threshold that gives us num_outlier_samples
total_expression_outlier_genes, total_expression_zscore_threshold = get_total_expression_outlier_info(individuals, genes, num_outlier_samples, total_expression_outlier_file)

# For ASE outliers extract list of genes for which we have at least one outlier sample and extract pvalue threshold that gives us num_outlier_samples
ase_outlier_genes, ase_pvalue_threshold = get_ase_outlier_info(individuals, genes, num_outlier_samples, ase_outlier_file)

# Take union of outlier genes from each of three methods
outlier_genes = union_of_three_dictionaries(splicing_outlier_genes, total_expression_outlier_genes, ase_outlier_genes)

# Extract dictionary of (sample, gene) pairs for each class that are valid. Value of dictionary is 1 if outlier and 0 if inlier
splicing_outliers = get_splicing_outliers(outlier_genes, individuals, splicing_pvalue_threshold, splicing_outlier_file)
splicing_outliers_unique = get_splicing_outliers(splicing_outlier_genes, individuals, splicing_pvalue_threshold, splicing_outlier_file)

total_expression_outliers = get_total_expression_outliers(outlier_genes, individuals, total_expression_zscore_threshold, total_expression_outlier_file)
total_expression_outliers_unique = get_total_expression_outliers(total_expression_outlier_genes, individuals, total_expression_zscore_threshold, total_expression_outlier_file)

ase_outliers = get_ase_outliers(outlier_genes, individuals, ase_pvalue_threshold, ase_outlier_file)
ase_outliers_unique = get_ase_outliers(ase_outlier_genes, individuals, ase_pvalue_threshold, ase_outlier_file)



# Print outliers
splicing_output_file = unsupervised_learning_input_dir + 'splicing_outliers_' + str(num_outlier_samples) + '_comparison_between_methods.txt'
print_outlier_output_file(splicing_outliers, genomic_annotation_file, splicing_output_file)

# Print outliers
total_expression_output_file = unsupervised_learning_input_dir + 'total_expression_outliers_' + str(num_outlier_samples) + '_comparison_between_methods.txt'
print_outlier_output_file(total_expression_outliers, genomic_annotation_file, total_expression_output_file)

# Print outliers
ase_output_file = unsupervised_learning_input_dir + 'ase_outliers_' + str(num_outlier_samples) + '_comparison_between_methods.txt'
print_outlier_output_file(ase_outliers, genomic_annotation_file, ase_output_file)


#####UNIQUE
# Print outliers
splicing_output_file = unsupervised_learning_input_dir + 'splicing_unique_outliers_' + str(num_outlier_samples) + '_comparison_between_methods.txt'
print_outlier_output_file(splicing_outliers_unique, genomic_annotation_file, splicing_output_file)

# Print outliers
total_expression_output_file = unsupervised_learning_input_dir + 'total_expression_unique_outliers_' + str(num_outlier_samples) + '_comparison_between_methods.txt'
print_outlier_output_file(total_expression_outliers_unique, genomic_annotation_file, total_expression_output_file)

# Print outliers
ase_output_file = unsupervised_learning_input_dir + 'ase_unique_outliers_' + str(num_outlier_samples) + '_comparison_between_methods.txt'
print_outlier_output_file(ase_outliers_unique, genomic_annotation_file, ase_output_file)

