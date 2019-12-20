import numpy as np 
import os
import sys
import pdb
import pickle
import dirichlet_multinomial_glm as dm_glm
import scipy.stats

def extract_junction_count_matrix_for_cluster(desired_cluster_id, whole_blood_junction_file, variant_position):
	f = open(whole_blood_junction_file)
	jxn_positions = []
	jxn_counts = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			sample_names = np.asarray(data[1:])
			continue
		jxn_id = data[0]
		cluster_id = jxn_id.split(':')[3]
		if cluster_id != desired_cluster_id:
			continue
		counts = np.asarray(data[1:]).astype(int)
		jxn_counts.append(counts)
		jxn_start = int(jxn_id.split(':')[1])
		jxn_end = int(jxn_id.split(':')[2])
		disty = np.min((np.abs(variant_position-jxn_start), np.abs(variant_position-jxn_end)))
		jxn_positions.append(disty)
	jxn_counts = np.transpose(np.asmatrix(jxn_counts))
	if len(sample_names) != jxn_counts.shape[0]:
		print('assumption error')
		pdb.set_trace()
	if len(jxn_positions) != jxn_counts.shape[1]:
		print('assumption error')
		pdb.set_trace()
	return jxn_counts, sample_names, np.argmin(jxn_positions)

def extract_outlier_jxn_counts(jxn_counts, sample_names, outlier_samples):
	arr = []
	for i, sample_name in enumerate(sample_names):
		if sample_name in outlier_samples:
			arr.append(i)
	arr = np.asarray(arr)
	return jxn_counts[arr,:]

def simulate_inlier_counts(num_reads, num_replicates, alpha):
	inlier_counts = dm_glm.draw_samples_from_dm(alpha, num_reads, num_replicates)
	return inlier_counts

def simulate_outlier_counts(num_reads, num_replicates, efficiency, alpha, outlier_alpha):
	inlier_counts = dm_glm.draw_samples_from_dm(alpha, num_reads, num_replicates)
	outlier_counts = dm_glm.draw_samples_from_dm(outlier_alpha, num_reads, num_replicates)

	mixed_counts = ((1.0-efficiency)*inlier_counts) + (efficiency*outlier_counts)
	return np.round(mixed_counts)

def compute_chi_squared_pvalue(simulated_inlier_counts, simulated_outlier_counts):
	inlier_sum = np.sum(simulated_inlier_counts,axis=0)
	outlier_sum = np.sum(simulated_outlier_counts, axis=0)
	contingency_table = np.vstack((inlier_sum,outlier_sum))
	nonzero_elements = np.where(np.sum(contingency_table,axis=0)!=0)[1]
	test = scipy.stats.chi2_contingency(contingency_table + 1)
	return test[1]

def compute_fisher_exact_pvalue(simulated_inlier_counts, simulated_outlier_counts, nearby_jxn):
	inlier_sum = np.sum(simulated_inlier_counts,axis=0)
	outlier_sum = np.sum(simulated_outlier_counts, axis=0)

	contingency_table = np.zeros((2,2))
	contingency_table[0,0] = inlier_sum[0, nearby_jxn]
	contingency_table[1,0] = outlier_sum[0, nearby_jxn]
	contingency_table[0,1] = np.sum(inlier_sum) - inlier_sum[0, nearby_jxn]
	contingency_table[1,1] = np.sum(outlier_sum) - outlier_sum[0, nearby_jxn]
	contingency_table = contingency_table + 1
	test = scipy.stats.fisher_exact(contingency_table)
	return test[1]


def compute_dm_pvalue(simulated_inlier_counts, simulated_outlier_counts, DM_GLM):
	num_reps = simulated_inlier_counts.shape[0]
	labels = np.zeros(num_reps*2)
	labels[:num_reps] = 1
	counts = np.vstack((simulated_inlier_counts, simulated_outlier_counts))

	pvalue = dm_glm.dirichlet_multinomial_differential_splicing(counts, labels, DM_GLM)
	return pvalue


def run_power_analysis_for_one_cluster(cluster_id, outlier_samples, variant_position, whole_blood_junction_file, t):
	# Stan model to fit Dirichlet multinomial
	DM_GLM = pickle.load(open('/home-1/bstrobe1@jhu.edu/scratch/gtex_v8/rare_var/gtex_v8_rare_splice/outlier_calling/dm_glm_multi_conc_no_prior.pkl', 'rb'))
	# Load in Dirichlet multinomial
	jxn_counts, sample_names, nearby_jxn = extract_junction_count_matrix_for_cluster(cluster_id, whole_blood_junction_file, variant_position)
	outlier_jxn_counts = extract_outlier_jxn_counts(jxn_counts, sample_names, outlier_samples) + 1e-19

	# Fit Dirichlet multinomial to jxn coutns
	alpha = dm_glm.dirichlet_multinomial_fit(jxn_counts, DM_GLM)

	# Get alpha for outlier distribution
	outlier_multinomial_mean = np.squeeze(np.asarray(np.sum(outlier_jxn_counts,axis=0)))/float(np.sum(outlier_jxn_counts))
	outlier_alpha = outlier_multinomial_mean*np.sum(alpha)

	# Simulation settings
	num_reps = 40
	num_reads_arr = [30, 50, 100, 500]
	num_replicates_arr = [2, 3, 4]
	efficiency_arr = [.05, .1, .2, .3, .4, .5]
	for num_reads in num_reads_arr:
		for num_replicates in num_replicates_arr:
			for efficiency in efficiency_arr:
				for rep in range(num_reps):
					simulated_inlier_counts = simulate_inlier_counts(num_reads, num_replicates, alpha)
					simulated_outlier_counts = simulate_outlier_counts(num_reads, num_replicates, efficiency, alpha, outlier_alpha)

					contingency_pvalue = compute_chi_squared_pvalue(simulated_inlier_counts, simulated_outlier_counts)
					fisher_exact_pvalue = compute_fisher_exact_pvalue(simulated_inlier_counts, simulated_outlier_counts, nearby_jxn)
					#dm_pvalue = compute_dm_pvalue(simulated_inlier_counts, simulated_outlier_counts, DM_GLM)

					t.write(cluster_id + '\t' + str(num_reads) + '\t' + str(num_replicates) + '\t' + str(efficiency) + '\t' + str(contingency_pvalue) + '\t' + str(fisher_exact_pvalue) + '\n')
	return t






###################
# Input data
###################
cluster_info_file = sys.argv[1]
whole_blood_junction_file = sys.argv[2]
test_variant_file = sys.argv[3]
power_analysis_data_output_dir = sys.argv[4]

t = open(power_analysis_data_output_dir + 'simulation_results.txt','w')
t.write('cluster_id\tnum_reads\tnum_replicates\tefficiency\tchi_squared_pvalue\tfisher_exact_pvalue\n')

# Loop through clusters
f = open(test_variant_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	cluster_id = data[5]
	outlier_samples = data[6].split(',')
	variant_position = int(data[0].split(':')[1])
	t = run_power_analysis_for_one_cluster(cluster_id, outlier_samples, variant_position, whole_blood_junction_file, t)
