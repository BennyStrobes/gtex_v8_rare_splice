import numpy as np 
import os
import sys
import pdb


def get_gene_name(stringer):
	gene_name = 'nuller'
	fields = stringer.split(';')
	for field in fields:
		field_info = field.split(' ')
		if field_info[0] == 'gene_id':
			gene_name = field_info[1].split('"')[1]
	if gene_name == 'nuller':
		print('get_gene_name assumption erroro')
	if gene_name.startswith('ENSG') == False:
		print('get_gene_name assumption erroro')
	return gene_name

def extract_list_of_valid_ensamble_ids(gene_list):
	valid_genes = {} # initilize list
	# Stream input file
	f = open(gene_list)
	for line in f:
		line = line.rstrip()
		data = line.split()
		valid_genes[data[0]] = 1
	f.close()
	return valid_genes



def get_annotated_splice_sites(gencode_gene_annotation_file, gene_list):
	annotated_splice_sites = {}
	valid_genes = extract_list_of_valid_ensamble_ids(gene_list)
	f = open(gencode_gene_annotation_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if line.startswith('#'):  # ignore header lines
			continue
		# error checking
		if len(data) != 9:
			print('length error in gencode file')
			pdb.set_trace()
		# Extract relevent fields
		line_chrom_num = data[0]
		gene_part = data[2]
		gene_name = get_gene_name(data[8])
		start = data[3]
		end = data[4]
		if gene_part != 'exon':
			continue
		if gene_name not in valid_genes:
			continue
		annotated_splice_sites[line_chrom_num + '_' + start] = 1
		annotated_splice_sites[line_chrom_num + '_' + end] = 1
	f.close()
	return annotated_splice_sites


#filter junction file to contain only junctions that are annotated
def filter_junction_file(tissue_specific_junction_file, annotated_junction_file, annotated_splice_sites):
	# Input file
	f = open(tissue_specific_junction_file)
	# Output file
	t = open(annotated_junction_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		jxn_identifier = data[0]
		chrom_num = jxn_identifier.split(':')[0]
		start = jxn_identifier.split(':')[1]
		end = jxn_identifier.split(':')[2]
		ss1 = chrom_num + '_' + start
		ss2 = chrom_num + '_' + end
		if ss1 in annotated_splice_sites or ss2 in annotated_splice_sites:
			t.write(line + '\n')
	f.close()
	t.close()


#call outliers junction by junction
def call_outliers_jxn_level(junction_object,input_file,outlier_calling_file_jxn_level,ratio):
	f = open(input_file)
	t = open(outlier_calling_file_jxn_level,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0: #header
			head_count = head_count + 1
			indiz = np.asarray(data[1:])
			t.write('JXNS\t' + '\t'.join(indiz) + '\n')
			continue
		id_info = data[0].split(':')
		line_chrom_num = id_info[0]
		start = id_info[1]
		end = id_info[2]
		intron_id = line_chrom_num + '_' + start + '_' + end
		#two splice sites for intron
		ss1 = line_chrom_num + '_' + start
		ss2 = line_chrom_num + '_' + end
		#Ignore introns that have no annotated splice sites
		if ss1 not in junction_object and ss2 not in junction_object:
			continue
		counts = np.asarray(data[1:]).astype(float)
		if ss1 in junction_object and ss2 not in junction_object:
			normalizer_counts = junction_object[ss1]
		elif ss1 not in junction_object and ss2 in junction_object:
			normalizer_counts = junction_object[ss2]
		elif ss1 in junction_object and ss2 in junction_object:
			normalizer_counts = np.maximum(junction_object[ss1],junction_object[ss2])
		else:
			print('erororroor')
			pdb.set_trace()

		normalized_counts = counts/normalizer_counts

		#if normalizer_count support == 0, then set normalized_value to zero
		for ele in np.where(normalizer_counts ==0)[0]:
			normalized_counts[ele] = 0
		max_index = np.argmax(normalized_counts)
		max_value = normalized_counts[max_index]

		#normalized value in extreme junction is more than twice the second largest value
		pvalues = np.ones(len(indiz))
		if len(np.where(normalized_counts > (max_value/ratio))[0]) ==1:
			#has more than 5% read support from annotated junctions and has more than 2 read support
			if normalized_counts[max_index] > .05 and counts[max_index] >= 2:
				pvalues[max_index] = 0
		t.write(data[0] + '\t' + '\t'.join(pvalues.astype(str)) + '\n')
	t.close()
	f.close()

def initialize_ss_in_jxn_object(jxn_object,ss,indiz):
	jxn_object[ss] = np.zeros(len(indiz))
	return jxn_object

#First create object that is dicti[chrom + jxn_pos][indi] = arr where arr contains read counts of every annotated intron overlapping this junction
def create_junction_object(input_file, annotated_splice_sites):
	f = open(input_file)
	head_count = 0
	jxn_object = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0: #header
			head_count = head_count + 1
			indiz = np.asarray(data[1:])
			continue
		id_info = data[0].split(':')
		line_chrom_num = id_info[0]
		start = id_info[1]
		end = id_info[2]
		intron_id = line_chrom_num + '_' + start + '_' + end
		#two splice sites for intron
		ss1 = line_chrom_num + '_' + start
		ss2 = line_chrom_num + '_' + end
		#Skip un-annotated intron
		if ss1 not in annotated_splice_sites or ss2 not in annotated_splice_sites:
			continue
		#Read counts overlapping this exon exon junction
		counts = np.asarray(data[1:]).astype(float)
		#Initialize ss in jxn_object if that ss has never been observed
		if ss1 not in jxn_object:
			jxn_object = initialize_ss_in_jxn_object(jxn_object,ss1,indiz)
		if ss2 not in jxn_object:
			jxn_object = initialize_ss_in_jxn_object(jxn_object,ss2,indiz)
		#Add sample specific counts to jxn_object[ss]
		jxn_object[ss1] = np.maximum(jxn_object[ss1],counts)
		jxn_object[ss2] = np.maximum(jxn_object[ss2],counts)
	return jxn_object

def heuristic_filtering(annotated_junction_file ,annotated_splice_sites, ratio, outlier_calls_file):
	#First create object that is dicti[chrom + jxn_pos]= arr where arr contains vector of len(num_samples), where each element is the maximum number of reads that known exon-exon junctions overlap that intron
	junction_object = create_junction_object(annotated_junction_file, annotated_splice_sites)
	#call outliers junction by junction
	call_outliers_jxn_level(junction_object, annotated_junction_file, outlier_calls_file, ratio)


tissue_type = sys.argv[1]
tissue_specific_junction_file = sys.argv[2]
output_root = sys.argv[3]
gencode_gene_annotation_file = sys.argv[4]
gene_list = sys.argv[5]

# Get list of annotated splice sites
annotated_splice_sites = get_annotated_splice_sites(gencode_gene_annotation_file, gene_list)

# Filter junction file to contain only junctions with at least one annotated splice site
annotated_junction_file = output_root + 'annotated_junctions_only.txt'
filter_junction_file(tissue_specific_junction_file, annotated_junction_file, annotated_splice_sites)


ratio=2
#Run main analysis for heuristic filtering (cummings et al).
#We will call outliers if:
######Normalized count is > .05
######Normalized count for junction is maximum in the sample
######Normalized count for junction,sample is twice that of any other sample
outlier_calls_file = output_root + 'outlier_calls.txt'
heuristic_filtering(annotated_junction_file ,annotated_splice_sites, ratio, outlier_calls_file)