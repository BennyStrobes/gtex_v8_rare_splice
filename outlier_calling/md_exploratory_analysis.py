import numpy as np 
import os
import sys
import pdb
import pickle
import dirichlet_multinomial_glm as dm_glm










DM_GLM = pickle.load(open('/home-1/bstrobe1@jhu.edu/scratch/gtex_v8/rare_var/gtex_v8_rare_splice/outlier_calling/dm_glm_multi_conc.pkl', 'rb'))
junction_count_file = sys.argv[1]


X = np.loadtxt(junction_count_file, delimiter='\t')
num_samples = X.shape[0]
cov_mat = np.ones((X.shape[0], 1))
num_reads=20000
mahalanobis_distances, pvalues, alpha = dm_glm.run_dm_outlier_analysis(X, cov_mat, DM_GLM, num_reads)


mds = []
counts = []
for num_reads in range(1,100):
	md = dm_glm.compute_mahalanobis_distance([num_reads,0], alpha)
	mds.append(md)
	counts.append(num_reads)

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')

fig = plt.figure()
plt.plot(counts,mds, linewidth=.2, color='green')
plt.xlabel('N')
plt.ylabel('MD([N,0])')
fig.savefig('temp.png')
