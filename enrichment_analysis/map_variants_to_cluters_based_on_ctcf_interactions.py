import numpy as np 
import os
import sys
import pdb


def extract_list_of_positions_to_liftover(ctcf_interaction_file):
	f = open(ctcf_interaction_file)
	head_count = 0
	listy = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		pos1 = data[0] + '_' + data[1]
		pos2 = data[0] + '_' + data[2]
		pos3 = data[3] + '_' +  data[4]
		pos4 = data[3] + '_' + data[5]
		pos5 = data[6].split('_')[0] + '_' + data[6].split('_')[1]
		pos6 = data[6].split('_')[0] + '_' + data[6].split('_')[2]
		listy.append(pos1)
		listy.append(pos2)
		listy.append(pos3)
		listy.append(pos4)
		listy.append(pos5)
		listy.append(pos6)
	f.close()
	return np.unique(listy)

def make_liftover_input_bed_file(positions_to_liftover, liftover_input_bed_file):
	t = open(liftover_input_bed_file, 'w')
	for position in positions_to_liftover:
		chrom_string = position.split('_')[0]
		pos_start = int(position.split('_')[1])
		t.write('chr' + chrom_string + '\t' + str(pos_start) + '\t' + str(pos_start + 1) + '\n')
	t.close()

#Some jxns were not able to be mapped with liftover (lacked confidence). So first extract those unmapped jxns
def get_unmapped_positions(temporary_missing_file):
	f = open(temporary_missing_file)
	unmapped_jxns = {}
	for line in f:
		if line.startswith('#'):
			continue
		line = line.rstrip()
		data = line.split()
		jxn_id = data[0] + ':' + data[1] 
		unmapped_jxns[jxn_id] = 1
	return unmapped_jxns

def get_liftover_mapping(liftover_input_bed_file, liftover_output_bed_file, unmapped_positions):
	mapping = {}
	f = open(liftover_input_bed_file)
	t = open(liftover_output_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		namer = data[0] + ':' + data[1]
		if namer in unmapped_positions:
			continue
		hg38_line = t.next()
		hg38_data = hg38_line.rstrip().split()

		mapping[data[0].split('hr')[1] + '_' + data[1]] = hg38_data[0].split('hr')[1] + '_' + hg38_data[1]
	f.close()
	return mapping

#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file,output_file,missing_file,liftover_directory):
	stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg19ToHg38.over.chain.gz ' + output_file + ' ' + missing_file
	os.system(stringer)

def update_ctcf_interaction_file(ctcf_interaction_file, ctcf_interaction_file_hg38, liftover_mapping):
	f = open(ctcf_interaction_file)
	t = open(ctcf_interaction_file_hg38, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		pos1 = data[0] + '_' + data[1]
		pos2 = data[0] + '_' + data[2]
		pos3 = data[3] + '_' +  data[4]
		pos4 = data[3] + '_' + data[5]
		pos5 = data[6].split('_')[0] + '_' + data[6].split('_')[1]
		pos6 = data[6].split('_')[0] + '_' + data[6].split('_')[2]
		if pos1 in liftover_mapping and pos2 in liftover_mapping and pos3 in liftover_mapping and pos4 in liftover_mapping and pos5 in liftover_mapping and pos6 in liftover_mapping:
			t.write(data[0] + '\t')
			new_pos1 = liftover_mapping[pos1]
			new_pos2 = liftover_mapping[pos2]
			new_pos3 = liftover_mapping[pos3]
			new_pos4 = liftover_mapping[pos4]
			new_pos5 = liftover_mapping[pos5]
			new_pos6 = liftover_mapping[pos6]
			t.write(new_pos1.split('_')[1] + '\t')
			t.write(new_pos2.split('_')[1] + '\t')
			t.write(data[3] + '\t')
			t.write(new_pos3.split('_')[1] + '\t')
			t.write(new_pos4.split('_')[1] + '\t')
			exon_name = data[6].split('_')[0] + '_' + new_pos5.split('_')[1] + '_' + new_pos6.split('_')[1] + '_' + data[6].split('_')[-1]
			t.write(exon_name + '\t' + '\t'.join(data[7:]) + '\n')
	f.close()
	t.close()


def liftover_file_to_hg38(ctcf_interaction_file, liftover_directory, output_dir, ctcf_interaction_file_hg38):
	positions_to_liftover = extract_list_of_positions_to_liftover(ctcf_interaction_file)
	# make input bed file for liftover
	liftover_input_bed_file = output_dir + 'liftover_hg19_to_hg38_ctcf_interaction_input.txt'
	make_liftover_input_bed_file(positions_to_liftover, liftover_input_bed_file)
	liftover_output_bed_file = output_dir  + 'temp_hg38.bed' #temporary liftover output file
	temporary_missing_file = output_dir  + 'temp_liftover_missing.bed' #temporary liftover missing values file
	run_liftover(liftover_input_bed_file,liftover_output_bed_file,temporary_missing_file,liftover_directory)

	unmapped_positions = get_unmapped_positions(temporary_missing_file)

	liftover_mapping = get_liftover_mapping(liftover_input_bed_file, liftover_output_bed_file, unmapped_positions)
	update_ctcf_interaction_file(ctcf_interaction_file, ctcf_interaction_file_hg38, liftover_mapping)


def get_exon_to_cluster_mapping(ctcf_interaction_file_hg38, cluster_info_file):
	ss_to_cluster = {}
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
		for jxn in jxns:
			ss1 = jxn.split(':')[0] + ':' + str(int(jxn.split(':')[1]))
			ss2 = jxn.split(':')[0] + ':' +  str(int(jxn.split(':')[2]))
			if ss1 not in ss_to_cluster:
				ss_to_cluster[ss1] = cluster_id
			else:
				if ss_to_cluster[ss1] != cluster_id:
					print('assumption error')
			if ss2 not in ss_to_cluster:
				ss_to_cluster[ss2] = cluster_id
			else:
				if ss_to_cluster[ss2] != cluster_id:
					print('assumption error')
	f.close()
	exon_to_cluster = {}
	f = open(ctcf_interaction_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		exon_id = data[6]
		ss1 = 'chr' + exon_id.split('_')[0] + ':' + str(int(exon_id.split('_')[1]) - 1)
		ss2 = 'chr' + exon_id.split('_')[0] + ':' + str(int(exon_id.split('_')[2]))
		if ss1 in ss_to_cluster and ss2 in ss_to_cluster:
			if exon_id not in exon_to_cluster:
				exon_to_cluster[exon_id] = []
			exon_to_cluster[exon_id].append(ss_to_cluster[ss1])
			exon_to_cluster[exon_id].append(ss_to_cluster[ss2])
		elif ss1 in ss_to_cluster and ss2 not in ss_to_cluster:
			if exon_id not in exon_to_cluster:
				exon_to_cluster[exon_id] = []
			exon_to_cluster[exon_id].append(ss_to_cluster[ss1])
		elif ss1 not in ss_to_cluster and ss2 in ss_to_cluster:
			if exon_id not in exon_to_cluster:
				exon_to_cluster[exon_id] = []
			exon_to_cluster[exon_id].append(ss_to_cluster[ss2])

	f.close()
	print(len(exon_to_cluster))
	pdb.set_trace()



def get_exon_to_gene_mapping(ctcf_interaction_file_hg38, cluster_info_file):
	valid_genes = {}
	head_count = 0
	f = open(cluster_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_id = data[0]
		genes = data[2].split(',')
		for gene in genes:
			valid_genes[gene] = 1
	f.close()
	exon_to_gene = {}
	f = open(ctcf_interaction_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		exon_id = data[6]
		ensamble_id = data[8]
		if ensamble_id not in valid_genes:
			continue
		if exon_id in exon_to_gene:
			if exon_to_gene[exon_id] != ensamble_id:
				print('assumption error')
				pdb.set_trace()
		exon_to_gene[exon_id] = ensamble_id
	f.close()
	return exon_to_gene

def fill_in_chromosome(chrom, start, end, ensamble_id):
	if end < start:
		print('assumption error')
		pdb.set_trace()
	for pos in range(start, end):
		if chrom[pos] != 'NULL' and chrom[pos] != ensamble_id:
			print('assumption error')
			pdb.set_trace()
		chrom[pos] = ensamble_id
	return chrom

def make_chromosome(ctcf_interaction_file_hg38, chrom_num, exon_to_gene):
	# Initialize
	chrom = ['NULL']*259250621
	chrom_string = str(chrom_num)
	# Stream cluster file
	head_count = 0
	f = open(ctcf_interaction_file_hg38)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[8]
		exon_id = data[6]
		chromer = data[0]
		if chromer != chrom_string:
			continue
		if exon_id not in exon_to_gene:
			continue
		if exon_to_gene[exon_id] != ensamble_id:
			print('assumption error')
		start = int(data[1])
		end = int(data[2])
		chrom = fill_in_chromosome(chrom,start,end, ensamble_id)
	return chrom

# Stream variant bed file. For each variant on the current chromosome, try to map to a cluster
def stream_variant_bed_file(cluster_chromosome, variant_bed_file, t, chrom_num):
	chrom_string = 'chr' + str(chrom_num)
	# used to skip header
	head_count = 0
	# Stream variant file
	f = open(variant_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
		line_chrom_num = data[0]
		# Throw out variants not on current chromosome
		if line_chrom_num != chrom_string:
			continue
		# Make sure MAF is not >= .01
		maf = float(data[4])
		if maf >= .01:
			print('maf assumption error!')
			continue
		var_pos = int(data[1])  # Position of variant
		# Get list of cluster_ids that overlap the variant position
		overlapping_clusters = cluster_chromosome[var_pos].split(',')
		# Print one line for every variant-cluster pair
		# If no clusters map, put one line and a NULL where cluster_id goes
		for cluster_id in overlapping_clusters:
			if cluster_id != 'NULL':
				t.write(line + '\t' + cluster_id + '\n')
	f.close()
	return t

variant_bed_file = sys.argv[1]
cluster_info_file = sys.argv[2]
ctcf_interaction_file = sys.argv[3]
liftover_directory = sys.argv[4]
output_dir = sys.argv[5]


ctcf_interaction_file_hg38 = output_dir + 'ctcf_interaction_file_hg38.txt'
#liftover_file_to_hg38(ctcf_interaction_file, liftover_directory, output_dir, ctcf_interaction_file_hg38)

# Map exons to clusters
#exon_to_clusters = get_exon_to_cluster_mapping(ctcf_interaction_file_hg38, cluster_info_file)
exon_to_gene = get_exon_to_gene_mapping(ctcf_interaction_file_hg38, cluster_info_file)

t = open(output_dir + 'ctcf_interaction_variants_mapped_to_genes.txt', 'w')

# Loop through autosomal chromsosomes
for chrom_num in range(1,23):
	print(chrom_num)
	# Initialize object to keep track of which BP on the chromosome correspond to which cluster
	cluster_chromosome = {}
	# Fill in the object for the current chromosome
	cluster_chromosome = make_chromosome(ctcf_interaction_file_hg38, chrom_num, exon_to_gene)
	t = stream_variant_bed_file(cluster_chromosome, variant_bed_file, t, chrom_num)
