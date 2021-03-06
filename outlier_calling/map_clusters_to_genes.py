import numpy as np
import os
import sys
import pdb
import gzip

def extract_tissues(tissue_list_input_file):
    f = open(tissue_list_input_file)
    arr = []
    for line in f:
        line = line.rstrip()
        arr.append(line)
    return arr


#Fill in chromosome object from chromosome[start:end+1] with the addition of this gene name
def add_gene_to_chromosome_object(chromosome, start, end, gene_name):
    for pos in range(start, end+1):
        if chromosome[pos] == 'NULL':  # No genes already mapped to this position
            chromosome[pos] = gene_name
        else:  # at least one gene already mapped to this position
            chromosome[pos] = chromosome[pos] + ',' + gene_name
    return chromosome

#Fill in chromosome object from chromosome[start:end+1] with the addition of this gene name
def add_position_to_chromosome_object(chromosome, pos, gene_name):
    if chromosome[pos] == 'NULL':  # No genes already mapped to this position
        chromosome[pos] = gene_name
    else:  # at least one gene already mapped to this position
        chromosome[pos] = chromosome[pos] + ',' + gene_name
    return chromosome

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


#Create an array of lenth(chromosome) [in bp]. Value of an element corresponds to a list of genes that overlap that basepair. This is used for efficient searching.
def make_chromosome_with_gene_names(chrom_num, gencode_hg19_gene_annotation_file, valid_genes):
    #initialize chromosome array
    chromosome = ['NULL']*259250621
    #loop through gencode file
    counter = 0
    f = open(gencode_hg19_gene_annotation_file)
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if line.startswith('#'):  # ignore header lines
            continue
        # Error Checking
        if len(data) != 9:
            print('length error in gencode file')
            pdb.set_trace()
        counter = counter + 1
        #print(counter)
        # Extract relevent fields
        line_chrom_num = data[0]
        gene_part = data[2]
        gene_name = get_gene_name(data[8])
        start = int(data[3])  # 5' gene start site (this is the min of all UTR,exons,etc)
        end = int(data[4])  # 3' gene end site (this is the max of all UTR,exons,etc)
        # limit to chromosome of interest
        if line_chrom_num != 'chr' + chrom_num:
            continue
        # Only consider lines for whole gene
        if gene_part != 'gene':
            continue
        # ignore genes not in our list
        if gene_name not in valid_genes:
            continue
        # Error checking
        if end < start:
            print('order error in gencode file')
            pdb.set_trace()
        #Fill in chromosome object from chromosome[start:end+1] with the addition of this gene name
        chromosome = add_gene_to_chromosome_object(chromosome, start, end, gene_name)
        #NOTE: ENSAMBLE Id's never come up more than once when gene_part == 'gene'
    f.close()
    return chromosome

#Create an array of lenth(chromosome) [in bp]. Value of an element corresponds to a list of genes that overlap that basepair. This is used for efficient searching.
def make_chromosome_with_gene_names_with_exon(chrom_num, gencode_hg19_gene_annotation_file, valid_genes):
    #initialize chromosome array
    chromosome = ['NULL']*259250621
    parts = {}
    #loop through gencode file
    f = open(gencode_hg19_gene_annotation_file)
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if line.startswith('#'):  # ignore header lines
            continue
        # Error Checking
        if len(data) != 9:
            print('length error in gencode file')
            pdb.set_trace()
        # Extract relevent fields
        line_chrom_num = data[0]
        gene_part = data[2]
        gene_name = get_gene_name(data[8])
        start = int(data[3])  # 5' gene start site (this is the min of all UTR,exons,etc)
        end = int(data[4])  # 3' gene end site (this is the max of all UTR,exons,etc)
        # Skip genes not in our valid gene set
        if gene_name not in valid_genes:
            continue
        if line_chrom_num != 'chr' + chrom_num:  # limit to chromosome of interest
            continue
        if gene_part != 'exon':  # ignore other parts of the gene as 'gene' encomposes everything
            continue
        start = int(data[3])  # 5' gene start site (this is the min of all UTR,exons,etc)
        end = int(data[4])  # 3' gene end site (this is the max of all UTR,exons,etc)
        # Ignore genes not in our list
        if gene_name not in valid_genes:
            continue
        #Fill in chromosome object from chromosome[start:end+1] with the addition of this gene name
        chromosome = add_position_to_chromosome_object(chromosome, start, gene_name)
        chromosome = add_position_to_chromosome_object(chromosome, end, gene_name)
        #NOTE: ENSAMBLE Id's never come up more than once when gene_part == 'gene'
    f.close()
    return chromosome


def get_unique_genes(stringer):
    info = stringer.split(',')
    return ','.join(np.unique(info))


#Loop through jxn file. For each jxn, see what genes overlap is 5' ss and 3' ss. And the union of the genes to cluster_to_genes_mapping for the jxns associated cluster
def align_clusters_to_genes(chromosome, chrom_num, cluster_to_genes_mapping, input_file_name):
    f = open(input_file_name)
    count = 0  # header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if count == 0:  # ignore header
            count = count + 1
            continue
        junction_info = data[0].split(':')
        line_chrom_num = junction_info[0]
        if line_chrom_num != 'chr' + chrom_num:
            continue
        start = int(junction_info[1])
        end = int(junction_info[2])
        cluster_id = junction_info[3]
        genes_at_start = chromosome[start]
        genes_at_end = chromosome[end]
        if genes_at_start == 'NULL' and genes_at_end == 'NULL':  # Both ends of the jxn overlap zero genes
            continue
        elif genes_at_start == 'NULL':
            gene_string = genes_at_end
        elif genes_at_end == 'NULL':
            gene_string = genes_at_start
        else:
            gene_string = genes_at_start + ',' + genes_at_end
        if cluster_id not in cluster_to_genes_mapping:
            cluster_to_genes_mapping[cluster_id] = get_unique_genes(gene_string)
        else:
            cluster_to_genes_mapping[cluster_id] = get_unique_genes(gene_string + ',' + cluster_to_genes_mapping[cluster_id])
    return cluster_to_genes_mapping


#Perform mapping analysis on each chromosome seperately
#At each chromosome update the cluster_to_genes_mapping for clusters on that chromosome
def map_clusters_on_chromosome(chrom_num, cluster_to_genes_mapping, clusters_filter_output_dir, input_suffix, tissues, gencode_hg19_gene_annotation_file, valid_genes):
    #Create an array of lenth(chromosome) [in bp]. Value of an element corresponds to a list of genes that overlap that basepair. This is used for efficient searching.
    chromosome = make_chromosome_with_gene_names(chrom_num, gencode_hg19_gene_annotation_file, valid_genes)
    #Loop through jxn file. For each jxn, see what genes overlap is 5' ss and 3' ss. And the union of the genes to cluster_to_genes_mapping for the jxns associated cluster
    for tissue in tissues:
        input_file_name = clusters_filter_output_dir + tissue + input_suffix
        cluster_to_genes_mapping = align_clusters_to_genes(chromosome, chrom_num, cluster_to_genes_mapping, input_file_name)
    return cluster_to_genes_mapping

#Perform mapping analysis on each chromosome seperately
#At each chromosome update the cluster_to_genes_mapping for clusters on that chromosome
def map_clusters_on_chromosome_with_exon(chrom_num, cluster_to_genes_mapping, clusters_filter_output_dir, input_suffix, tissues, gencode_hg19_gene_annotation_file, valid_genes):
    #Create an array of lenth(chromosome) [in bp]. Value of an element corresponds to a list of genes that overlap that basepair. This is used for efficient searching.
    chromosome = make_chromosome_with_gene_names_with_exon(chrom_num, gencode_hg19_gene_annotation_file, valid_genes)
    #Loop through jxn file. For each jxn, see what genes overlap is 5' ss and 3' ss. And the union of the genes to cluster_to_genes_mapping for the jxns associated cluster
    for tissue in tissues:
        input_file_name = clusters_filter_output_dir + tissue + input_suffix
        cluster_to_genes_mapping = align_clusters_to_genes(chromosome, chrom_num, cluster_to_genes_mapping, input_file_name)
    return cluster_to_genes_mapping


def print_helper_map_clusters_to_genes(input_file_name, output_file_name, cluster_to_genes_mapping):
    f = open(input_file_name)
    t = open(output_file_name, 'w')
    count = 0  # header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if count == 0:
            count = count + 1
            t.write(line + '\n')
            continue
        junction_info = data[0].split(':')
        cluster_id = junction_info[3]
        if cluster_id not in cluster_to_genes_mapping:
            continue
        new_jxn_id = data[0] + ':' + cluster_to_genes_mapping[cluster_id]
        t.write(new_jxn_id + '\t' + '\t'.join(data[1:]) + '\n')
    t.close()

def map_clusters_to_genes_with_exons(tissues, gencode_hg19_file, clusters_filter_output_dir, input_suffix, output_suffix, valid_genes):
    #Main dictionary to keep track of cluster --> genes.
    #Key will be leafcutter cluster id. Value is string gene1,gene2,gene3,... Or 'NULL' if no gene is mapped
    cluster_to_genes_mapping = {}
    #Perform mapping analysis on each chromosome seperately
    #At each chromosome update the cluster_to_genes_mapping for clusters on that chromosome
    for chrom_num in range(1, 23):
        print(chrom_num)
        cluster_to_genes_mapping = map_clusters_on_chromosome_with_exon(str(chrom_num), cluster_to_genes_mapping, clusters_filter_output_dir, input_suffix, tissues, gencode_hg19_file, valid_genes)

    for tissue in tissues:
        print(tissue)
        input_file = clusters_filter_output_dir + tissue + input_suffix
        output_file = clusters_filter_output_dir + tissue + output_suffix
        print_helper_map_clusters_to_genes(input_file, output_file, cluster_to_genes_mapping)


def map_clusters_to_genes(tissues, gencode_hg19_file, clusters_filter_output_dir, input_suffix, output_suffix, valid_genes):
    #Main dictionary to keep track of cluster --> genes.
    #Key will be leafcutter cluster id. Value is string gene1,gene2,gene3,... Or 'NULL' if no gene is mapped
    cluster_to_genes_mapping = {}
    #Perform mapping analysis on each chromosome seperately
    #At each chromosome update the cluster_to_genes_mapping for clusters on that chromosome
    for chrom_num in range(1, 23):
        print(chrom_num)
        cluster_to_genes_mapping = map_clusters_on_chromosome(str(chrom_num), cluster_to_genes_mapping, clusters_filter_output_dir, input_suffix, tissues, gencode_hg19_file, valid_genes)

    for tissue in tissues:
        print(tissue)
        input_file = clusters_filter_output_dir + tissue + input_suffix
        output_file = clusters_filter_output_dir + tissue + output_suffix
        print_helper_map_clusters_to_genes(input_file, output_file, cluster_to_genes_mapping)


def report_cluster_info(tissues, clusters_filter_output_dir, output_suffix, output_file):
    clusters = {}  # clusters[cluster_id]['jxns'] yields array of all jxn ids that are mapped to that cluster. clusters[cluster_id]['genes'] contains an array of all ensemble id genes that the cluster is mapped to
    for tissue in tissues:
        file_name = clusters_filter_output_dir + tissue + output_suffix  # File that was just created in 'map_clusters_to_genes'
        f = open(file_name)
        count = 0
        for line in f:
            if count == 0:
                count = count + 1
                continue
            line = line.rstrip()
            data = line.split()

            junction_info = data[0].split(':')
            cluster_id = junction_info[3]
            jxn_name = junction_info[0] + ':' + junction_info[1] + ':' + junction_info[2]
            genes = junction_info[4].split(',')
            if cluster_id not in clusters:  # Cluster id has never been seen before
                clusters[cluster_id] = {}
                clusters[cluster_id]['jxns'] = []
                clusters[cluster_id]['jxns'].append(jxn_name)
                clusters[cluster_id]['genes'] = []
                for gene in genes:
                    clusters[cluster_id]['genes'].append(gene)
            else:  # Cluster id has been seen before. No need to initialize
                clusters[cluster_id]['jxns'].append(jxn_name)
                for gene in genes:
                    clusters[cluster_id]['genes'].append(gene)
    # Print mapping from clusters --> (jxns, genes) to an output file
    t = open(output_file, 'w')
    t.write('cluster_id\tjxns\tgenes\n')
    for cluster in sorted(clusters.keys()):
        genes = np.unique(clusters[cluster]['genes'])
        jxns = np.unique(clusters[cluster]['jxns'])
        t.write(cluster + '\t' + ','.join(jxns) + '\t' + ','.join(genes) + '\n')

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

# First extract dictionary list of ensamble ids used in our analysis
def extract_list_of_used_ensamble_ids(cluster_info_file):
    # Initialize output list
    ensamble_ids = {}
    # Used to skip header
    head_count = 0
    # Stream input file
    f = open(cluster_info_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        # Skip header
        if head_count == 0:
            head_count = head_count + 1
            continue
        # Get array of genes mapped to this cluster (use array b/c there can be more than one gene)
        gene_arr = data[2].split(',')
        for gene in gene_arr:
            ensamble_ids[gene] = 1
    f.close()
    return ensamble_ids

# Make text file summarizing exons
def make_text_file_summarizing_exons(cluster_info_file, gencode_hg19_file, exon_file):
    # First extract dictionary list of ensamble ids used in our analysis
    ensamble_ids = extract_list_of_used_ensamble_ids(cluster_info_file)
    # Open output file handle and print header
    t = open(exon_file, 'w')
    t.write('chr\tstart\tend\tstrand\tgene_name\n')

    f = open(gencode_hg19_file)
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if line.startswith('#'):  # ignore header lines
            continue
        # Error Checking
        if len(data) != 9:
            print('length error in gencode file')
            pdb.set_trace()
        # Extract relevent fields
        line_chrom_num = data[0]
        gene_part = data[2]
        gene_name = get_gene_name(data[8])
        strand = data[6]
        # Skip genes not in our valid gene set
        if gene_name not in ensamble_ids:
            continue
        if gene_part != 'exon':  # ignore other parts of the gene that are not the exon
            continue
        start = data[3]  # 5' gene start site (this is the min of all UTR,exons,etc)
        end = data[4]  # 3' gene end site (this is the max of all UTR,exons,etc)
        strand = data[6]
        t.write(line_chrom_num + '\t' + start + '\t' + end + '\t' + strand + '\t' + gene_name + '\n')
    t.close()
    f.close()

tissue_list_input_file = sys.argv[1]
clusters_filter_output_dir = sys.argv[2]  # Input dir and output dir
gencode_hg19_file = sys.argv[3]
gene_list = sys.argv[4]

tissues = extract_tissues(tissue_list_input_file)  # get array of tissue types

# Limit to genes used in other outlier analysis.
valid_genes = extract_list_of_valid_ensamble_ids(gene_list)

input_suffix = '_filtered_jxns_cross_tissue_clusters.txt'  # Suffix of input files
output_suffix = '_filtered_jxns_cross_tissue_clusters_full_gene_mapped.txt'  # Suffix of output files

map_clusters_to_genes(tissues, gencode_hg19_file, clusters_filter_output_dir, input_suffix, output_suffix, valid_genes)

# Make output file containing all clusters (generated in map_clusters_to_genes) as well as which junctions are mapped to which cluster, as well as which genes are mapped to which cluster
report_cluster_info(tissues, clusters_filter_output_dir, output_suffix, clusters_filter_output_dir + 'cluster_info_full_gene_mapped.txt')



input_suffix = '_filtered_jxns_cross_tissue_clusters.txt'  # Suffix of input files
output_suffix = '_filtered_jxns_cross_tissue_clusters_gene_mapped.txt'  # Suffix of output files

map_clusters_to_genes_with_exons(tissues, gencode_hg19_file, clusters_filter_output_dir, input_suffix, output_suffix, valid_genes)
report_cluster_info(tissues, clusters_filter_output_dir, output_suffix, clusters_filter_output_dir + 'cluster_info.txt')


# Make text file summarizing exons
cluster_info_file = clusters_filter_output_dir + 'cluster_info.txt'  # Input file (created in above lines)
exon_file = clusters_filter_output_dir + 'gencode_v26_exons.txt'  # Output file
make_text_file_summarizing_exons(cluster_info_file, gencode_hg19_file, exon_file)




