import numpy as np 
import os
import sys
import pdb
import scipy.stats
import gzip


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
	head_count = 0
	if ase_outlier_file.endswith('.gz'):
		f = gzip.open(ase_outlier_file)
	else:
		f = open(ase_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for ele in data[1:]:
				indi[ele] = 1
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
	if ase_outlier_file.endswith('.gz'):
		f = gzip.open(ase_outlier_file)
	else:
		f = open(ase_outlier_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes[data[0]] = 1
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
	ase_genes = get_ase_genes(ase_outlier_file)
	total_expression_genes = get_total_expression_genes(total_expression_outlier_file)
	union_genes = {}
	for gene in splicing_genes.keys():
		if gene in ase_genes and gene in total_expression_genes:
			union_genes[gene] = 1
	return union_genes

def convert_na_to_nan(data_input):
	data_output = []
	for ele in data_input:
		if ele == 'NA':
			data_output.append('NaN')
		else:
			data_output.append(ele)
	return data_output


# For splicing outliers extract list of genes for which we have at least one outlier sample and extract pvalue threshold that gives us num_outlier_samples
def get_ase_outlier_info(individuals, genes, pvalue_threshold, splicing_outlier_file):
	splicing_outlier_genes = {}
	head_count = 0
	if splicing_outlier_file.endswith('.gz'):
		f = gzip.open(splicing_outlier_file)
	else:
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
		pvalues = np.asarray(convert_na_to_nan(data[1:])).astype(float)
		# Remove columns not in individuals
		pvalues_filtered = pvalues[good_columns]
		for pvalue in pvalues_filtered:
			if np.isnan(pvalue):
				continue
			# Check if (individual, gene) is an outlier. If so, add gene to dictionary
			if pvalue < pvalue_threshold:
				splicing_outlier_genes[ensamble] = 1
	f.close()
	return splicing_outlier_genes


# For splicing outliers extract list of genes for which we have at least one outlier sample and extract pvalue threshold that gives us num_outlier_samples
def get_splicing_outlier_info(individuals, genes, pvalue_threshold, splicing_outlier_file):
	splicing_outlier_genes = {}
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
		pvalues = np.asarray(convert_na_to_nan(data[1:])).astype(float)
		# Remove columns not in individuals
		pvalues_filtered = pvalues[good_columns]
		for pvalue in pvalues_filtered:
			if np.isnan(pvalue):
				continue
			# Check if (individual, gene) is an outlier. If so, add gene to dictionary
			if pvalue < pvalue_threshold:
				splicing_outlier_genes[ensamble] = 1
	f.close()
	return splicing_outlier_genes


# For total expression outliers extract list of genes for which we have at least one outlier sample and extract pvalue threshold that gives us num_outlier_samples
def get_total_expression_outlier_info(individuals, genes, pvalue_threshold, total_expression_outlier_file):
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
		pvalue = scipy.stats.norm.sf(zscore)*2
		if pvalue < pvalue_threshold:
			total_expression_outlier_genes[ensamble_id] = 1
	f.close()
	return total_expression_outlier_genes


def union_of_three_dictionaries(dicti1, dicti2, dicti3):
	union = {}
	for key  in dicti1.keys():
		union[key] = 1
	for key  in dicti2.keys():
		union[key] = 1
	for key in dicti3.keys():
		union[key] = 1
	return union


def intersection_of_three_dictionaries(dicti1, dicti2, dicti3):
	intersection = {}
	for key  in dicti1.keys():
		if key in dicti2 and key in dicti3:
			intersection[key] = 1
	return intersection

def get_ase_outliers(genes, individuals, splicing_outlier_file):
	outlier_dicti = {}
	head_count = 0
	if splicing_outlier_file.endswith('.gz'):
		f = gzip.open(splicing_outlier_file)
	else:
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
		pvalues = np.asarray(convert_na_to_nan(data[1:])).astype(float)
		# Remove columns not in individuals
		pvalues_filtered = pvalues[good_columns]
		# Add to global pvalue_array
		for i,pvalue in enumerate(pvalues_filtered):
			indi_id = new_indiz[i]
			sample_name = indi_id + '_' + ensamble
			outlier_dicti[sample_name] = pvalue
	f.close()
	return outlier_dicti

def get_splicing_outliers(genes, individuals, splicing_outlier_file):
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
		pvalues = np.asarray(convert_na_to_nan(data[1:])).astype(float)
		# Remove columns not in individuals
		pvalues_filtered = pvalues[good_columns]
		# Add to global pvalue_array
		for i,pvalue in enumerate(pvalues_filtered):
			indi_id = new_indiz[i]
			sample_name = indi_id + '_' + ensamble
			outlier_dicti[sample_name] = pvalue
	f.close()
	return outlier_dicti

def get_total_expression_outliers(genes, individuals, total_expression_outlier_file):
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
		zscore = float(data[4])
		pvalue = scipy.stats.norm.sf(abs(zscore))*2
		sample_name = indi_id + '_' + ensamble_id
		if zscore < 0:
			pvalue = pvalue*-1
		outlier_dicti[sample_name] = pvalue
	f.close()
	for indi in individuals.keys():
		for gene in genes.keys():
			sample_name = indi + '_' + gene
			if sample_name not in outlier_dicti:
				outlier_dicti[sample_name] = float('NaN')
	return outlier_dicti



def print_outlier_output_file(outlier_dicti, genomic_annotation_file, output_file, gene_individual_to_variant_mapping_file, seeder):
	t = open(output_file, 'w')
	f = open(genomic_annotation_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'pvalue\n')
			continue
		indi_id = data[0]
		gene_id = data[1].split('.')[0]
		sample_name = indi_id + '_' + gene_id
		if sample_name in outlier_dicti:
			t.write(line + '\t' + str(outlier_dicti[sample_name]) + '\n')
	f.close()
	t.close()
	remove_no_variance_features(output_file)
	add_n2_pairs_column(output_file.split('.tx')[0] + '_features_filter.txt',gene_individual_to_variant_mapping_file, seeder)
	os.system('rm ' + output_file)
	os.system('rm ' + output_file.split('.tx')[0] + '_features_filter.txt')


def print_outlier_output_file_no_nan(outlier_dicti1, outlier_dicti2, outlier_dicti3, genomic_annotation_file, output_file, gene_individual_to_variant_mapping_file, seeder):
	t = open(output_file, 'w')
	f = open(genomic_annotation_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'pvalue\n')
			continue
		indi_id = data[0]
		gene_id = data[1].split('.')[0]
		sample_name = indi_id + '_' + gene_id
		if sample_name in outlier_dicti1 and sample_name in outlier_dicti2 and sample_name in outlier_dicti3:
			if np.isnan(outlier_dicti1[sample_name]) == False and np.isnan(outlier_dicti2[sample_name]) == False and np.isnan(outlier_dicti3[sample_name]) == False:
				t.write(line + '\t' + str(outlier_dicti1[sample_name]) + '\n')
	f.close()
	t.close()
	remove_no_variance_features(output_file)
	add_n2_pairs_column(output_file.split('.tx')[0] + '_features_filter.txt',gene_individual_to_variant_mapping_file, seeder)
	os.system('rm ' + output_file)
	os.system('rm ' + output_file.split('.tx')[0] + '_features_filter.txt')


def add_n2_pairs_column(input_file, gene_individual_to_variant_mapping_file, seeder):
	f = open(gene_individual_to_variant_mapping_file)
	test_to_variant_mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		test_name = data[0] + '_' + data[1]
		variants = data[1] + '_' + ','.join(sorted(data[2].split(',')))
		if test_name in test_to_variant_mapping:
			print('assumption error')
			pdb.set_trace()
		test_to_variant_mapping[test_name] = variants
	f.close()
	dicti = {}
	output_file = input_file.split('.tx')[0] + '_N2_pairs.txt'
	f = open(input_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header = line
			continue
		test_name = data[0] + '_' + data[1]
		variants_mapped = test_to_variant_mapping[test_name]
		liner_long = '\t'.join(data)
		if variants_mapped not in dicti:
			dicti[variants_mapped] = []
		dicti[variants_mapped].append(liner_long)
	f.close()
	np.random.seed(seeder)
	t = open(output_file,'w')
	t.write(header + '\tN2pair\n')
	for key in sorted(dicti.keys()):
		arr = dicti[key]
		if len(arr) == 1:
			liner_long = arr[0]
			t.write(liner_long + '\t' + 'NA' + '\n')
	n2_pair_counter = 0
	for key in sorted(dicti.keys()):
		arr = dicti[key]
		num_indis = len(arr)
		if num_indis > 1:
			n2_pair_counter = n2_pair_counter + 1
			for liner_long in np.random.choice(arr,size=2,replace=False):
				t.write(liner_long + '\t' + str(n2_pair_counter) + '\n')
	t.close()

def remove_no_variance_features(input_file):
	output_file = input_file.split('.tx')[0] + '_features_filter.txt'
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

def merge_three_files(input_file1, input_file2, input_file3, output_file):
	# Loop through input file1
	# Keep track of feature vectors as well as pvalue from outlier status
	file_1_feature_vectors = []
	file_1_pvalues = []
	f = open(input_file1)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		pvalue = data[-2]
		feature_vector = '\t'.join(data[0:-2])
		file_1_pvalues.append(pvalue)
		file_1_feature_vectors.append(feature_vector)
	f.close()
	file_2_feature_vectors = []
	file_2_pvalues = []
	f = open(input_file2)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		pvalue = data[-2]
		feature_vector = '\t'.join(data[0:-2])
		file_2_pvalues.append(pvalue)
		file_2_feature_vectors.append(feature_vector)
	f.close()


	t = open(output_file,'w')
	f = open(input_file3)
	head_count = 0
	counter = -1
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(data[0:-2]) + '\tsplicing_pvalue\ttotal_expression_pvalue\tase_pvalue\t' + data[-1] + '\n')
			continue
		pvalue = data[-2]
		feature_vector = '\t'.join(data[0:-2])
		counter = counter + 1
		if feature_vector != file_1_feature_vectors[counter] or feature_vector != file_2_feature_vectors[counter]:
			print('assumption errror')
			pdb.set_trace()
		t.write(feature_vector + '\t' + file_1_pvalues[counter] + '\t' + file_2_pvalues[counter] + '\t' + pvalue + '\t' + data[-1] + '\n')
	f.close()
	t.close()

# Make output file with no tissue specific genomic annotations
def remove_tissue_specific_annotations(input_file, output_file):
	f = open(input_file)
	t = open(output_file, 'w')
	for line in f:
		line = line.rstrip()
		data = line.split()
		t.write('\t'.join(data[0:-58]) + '\t' + data[-4] + '\t' + data[-3] + '\t' + data[-2] + '\t' + data[-1] + '\n')
	f.close()
	t.close()

def randomly_filter_training_instances(input_file, standard_watershed_file, output_file, N, seeder):
	np.random.seed(3)
	# First get current number of samples
	f = open(input_file)
	head_count = 0
	num_training = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[-1] == 'NA': # is a training instance
			num_training = num_training + 1
	f.close()
	if N > num_training:
		N = num_training
	# Get dictionary list of valid column names
	valid_column_names = {}
	head_count = 0
	f = open(standard_watershed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for ele in data:
				valid_column_names[ele] = 1
			header=data
			continue
	f.close()
	# Then randomly subset to N training samples
	randomly_selected_samples_arr = np.random.choice(range(num_training),size=N,replace=False)
	randomly_selected_samples_dict = {}
	for ele in randomly_selected_samples_arr:
		randomly_selected_samples_dict[ele] = 1
	f = open(input_file)
	t = open(output_file, 'w')
	head_count = 0
	line_counter = 0
	for line in f:
		line = line.rstrip()
		data = np.asarray(line.split())
		if head_count == 0:
			head_count = head_count + 1
			valid_column_positions = []
			for i,ele in enumerate(data):
				if ele in valid_column_names:
					valid_column_positions.append(i)
			valid_column_positions = np.asarray(valid_column_positions)
			if np.array_equal(header,data[valid_column_positions]) == False:
				print('assumption error; cant continue')
				pdb.set_trace()
			t.write('\t'.join(data[valid_column_positions]) + '\n')
			continue
		if data[-1] != 'NA':
			filtered_line = data[valid_column_positions]
			t.write('\t'.join(filtered_line) + '\n')
		if line_counter in randomly_selected_samples_dict:
			filtered_line = data[valid_column_positions]
			t.write('\t'.join(filtered_line) + '\n')
		line_counter = line_counter + 1
	f.close()
	t.close()
	return


def merge_with_standard_watershed_file(input_file, standard_watershed_file, output_file, max_number_of_training_instances, seeder):
	np.random.seed(seeder)
	# First get current number of samples
	f = open(input_file)
	head_count = 0
	num_training = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[-1] == 'NA': # is a training instance
			num_training = num_training + 1
	f.close()
	if max_number_of_training_instances > num_training:
		max_number_of_training_instances = num_training
	# Get dictionary list of valid column names
	valid_column_names = {}
	n2_pair_lines = []
	head_count = 0
	f = open(standard_watershed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for ele in data:
				valid_column_names[ele] = 1
			header=data
			continue
		# Skip lines that are not N2 pairs
		if data[-1] == 'NA':
			continue
		n2_pair_lines.append(line)
	f.close()
	# Then randomly subset to N training samples
	randomly_selected_samples_arr = np.random.choice(range(num_training),size=max_number_of_training_instances,replace=False)
	randomly_selected_samples_dict = {}
	for ele in randomly_selected_samples_arr:
		randomly_selected_samples_dict[ele] = 1
	# Print output file
	f = open(input_file)
	t = open(output_file, 'w')
	head_count = 0
	line_counter = 0
	for line in f:
		line = line.rstrip()
		data = np.asarray(line.split())
		if head_count == 0:
			head_count = head_count + 1
			valid_column_positions = []
			for i,ele in enumerate(data):
				if ele in valid_column_names:
					valid_column_positions.append(i)
			valid_column_positions = np.asarray(valid_column_positions)
			if np.array_equal(header,data[valid_column_positions]) == False:
				print('assumption error; cant continue')
				pdb.set_trace()
			t.write('\t'.join(data[valid_column_positions]) + '\n')
			continue
		if data[-1] == 'NA':
			if line_counter in randomly_selected_samples_arr:
				filtered_line = data[valid_column_positions]
				t.write('\t'.join(filtered_line) + '\n')
		line_counter = line_counter + 1
	f.close()
	# Add N2 pair lines from standard watershed file
	for n2_pair_line in n2_pair_lines:
		t.write(n2_pair_line + '\n')
	t.close()
	return


genomic_annotation_file = sys.argv[1]
total_expression_outlier_file = sys.argv[2]
ase_outlier_file = sys.argv[3]
splicing_outlier_file = sys.argv[4]
unsupervised_learning_input_dir = sys.argv[5]
pvalue = float(sys.argv[6])
gene_individual_to_variant_mapping_file = sys.argv[7]
seeder = int(sys.argv[8])

'''
# Extract dictionary list of individuals used in all three methods
individuals = get_list_of_individuals_used_in_all_methods(total_expression_outlier_file, ase_outlier_file, splicing_outlier_file)
# Extract dictionary list of genes used in all three methods
genes = get_list_of_genes_used_in_all_methods(total_expression_outlier_file, ase_outlier_file, splicing_outlier_file)


# For ase outliers extract list of genes for which we have at least one outlier sample and extract pvalue threshold that gives us num_outlier_samples
ase_outlier_genes = get_ase_outlier_info(individuals, genes, pvalue, ase_outlier_file)
# For splicing outliers extract list of genes for which we have at least one outlier sample and extract pvalue threshold that gives us num_outlier_samples
splicing_outlier_genes = get_splicing_outlier_info(individuals, genes, pvalue, splicing_outlier_file)
# For total expression outliers extract list of genes for which we have at least one outlier sample and extract zscore threshold that gives us num_outlier_samples
total_expression_outlier_genes = get_total_expression_outlier_info(individuals, genes, pvalue, total_expression_outlier_file)

# Take union of outlier genes from each of three methods
outlier_genes = union_of_three_dictionaries(splicing_outlier_genes, total_expression_outlier_genes, ase_outlier_genes)

# Extract dictionary of (sample, gene) pairs for each class that are valid. Value of dictionary is 1 if outlier and 0 if inlier
splicing_outliers = get_splicing_outliers(outlier_genes, individuals, splicing_outlier_file)
total_expression_outliers = get_total_expression_outliers(outlier_genes, individuals, total_expression_outlier_file)
ase_outliers = get_ase_outliers(outlier_genes, individuals, ase_outlier_file)



# Print outliers
ase_output_file = unsupervised_learning_input_dir + 'ase_outliers_' + str(pvalue) + '_genes_union_between_te_ase_splicing.txt'
print_outlier_output_file(ase_outliers, genomic_annotation_file, ase_output_file, gene_individual_to_variant_mapping_file, seeder)
splicing_output_file = unsupervised_learning_input_dir + 'splicing_outliers_' + str(pvalue) + '_genes_union_between_te_ase_splicing.txt'
print_outlier_output_file(splicing_outliers, genomic_annotation_file, splicing_output_file, gene_individual_to_variant_mapping_file, seeder)
total_expression_output_file = unsupervised_learning_input_dir + 'total_expression_outliers_' + str(pvalue) + '_genes_union_between_te_ase_splicing.txt'
print_outlier_output_file(total_expression_outliers, genomic_annotation_file, total_expression_output_file, gene_individual_to_variant_mapping_file, seeder)

# Merge three outlier files
merged_output_file = unsupervised_learning_input_dir + 'merged_outliers_' + str(pvalue) + '_genes_union_between_te_ase_splicing_features_filter_N2_pairs.txt'
merge_three_files(splicing_output_file.split('.tx')[0] + '_features_filter_N2_pairs.txt', total_expression_output_file.split('.tx')[0] + '_features_filter_N2_pairs.txt', ase_output_file.split('.tx')[0] + '_features_filter_N2_pairs.txt', merged_output_file)







# Print outliers
splicing_output_file = unsupervised_learning_input_dir + 'fully_observed_splicing_outliers_' + str(pvalue) + '_genes_union_between_te_ase_splicing.txt'
print_outlier_output_file_no_nan(splicing_outliers, total_expression_outliers,ase_outliers, genomic_annotation_file, splicing_output_file, gene_individual_to_variant_mapping_file, seeder)
total_expression_output_file = unsupervised_learning_input_dir + 'fully_observed_total_expression_outliers_' + str(pvalue) + '_genes_union_between_te_ase_splicing.txt'
print_outlier_output_file_no_nan(total_expression_outliers, splicing_outliers, ase_outliers, genomic_annotation_file, total_expression_output_file, gene_individual_to_variant_mapping_file, seeder)
ase_output_file = unsupervised_learning_input_dir + 'fully_observed_ase_outliers_' + str(pvalue) + '_genes_union_between_te_ase_splicing.txt'
print_outlier_output_file_no_nan(ase_outliers, total_expression_outliers, splicing_outliers, genomic_annotation_file, ase_output_file, gene_individual_to_variant_mapping_file, seeder)


# Merge three outlier files
merged_output_file = unsupervised_learning_input_dir + 'fully_observed_merged_outliers_' + str(pvalue) + '_genes_union_between_te_ase_splicing_features_filter_N2_pairs.txt'
merge_three_files(splicing_output_file.split('.tx')[0] + '_features_filter_N2_pairs.txt', total_expression_output_file.split('.tx')[0] + '_features_filter_N2_pairs.txt', ase_output_file.split('.tx')[0] + '_features_filter_N2_pairs.txt', merged_output_file)

'''
# Make output file with no tissue specific genomic annotations
merged_no_tissue_anno_output_file = unsupervised_learning_input_dir + 'fully_observed_merged_outliers_' + str(pvalue) + '_genes_union_between_te_ase_splicing_features_filter_no_tissue_anno_N2_pairs.txt'
#remove_tissue_specific_annotations(merged_output_file, merged_no_tissue_anno_output_file)

# Randomly subset the previous file (merged_no_tissue_anno_output_file) to N training instances
merged_no_tissue_anno_same_n2_pairs_output_file = unsupervised_learning_input_dir + 'fully_observed_merged_outliers_' + str(pvalue) + '_genes_union_between_te_ase_splicing_features_filter_no_tissue_anno_same_N2_pairs_as_standard_450000_' + str(seeder) + '.txt'
# Make sure it has the same features as standard watershed file
standard_watershed_file = unsupervised_learning_input_dir + 'fully_observed_merged_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_no_tissue_anno_N2_pairs_' + str(seeder) + '.txt'

# Use same N2 pairs as was used in standard watershed analysis (p=.01).. see standard_watershed_file
max_number_of_training_instances=450000
merge_with_standard_watershed_file(merged_no_tissue_anno_output_file, standard_watershed_file, merged_no_tissue_anno_same_n2_pairs_output_file, max_number_of_training_instances, seeder)





###### OLD (NO LONGER USED)
#randomly_filter_training_instances(merged_no_tissue_anno_output_file, standard_watershed_file, merged_no_tissue_anno_random_subset_output_file, seeder)
