import numpy as np
import pdb
import pystan
import pickle
from scipy.special import gamma, factorial, loggamma

###########################################################
#DM_GLM = pystan.StanModel(file = "dm_glm_multi_conc.stan")
# For quick loading
#f = open('dm_glm_multi_conc.pkl', 'wb')
#pickle.dump(DM_GLM, f)
#f.close()
#DM_GLM = pickle.load(open('dm_glm_multi_conc.pkl', 'rb'))
###########################################################


DC_GLM = pystan.StanModel(file = "dirichlet_categorical.stan")
#DC_GLM = pickle.load(open('dirichlet_categorical.pkl', 'rb'))

def correct_betas(beta_raw_object,beta_scale,K_input):
	corrected_betas = (beta_raw_object  - (1.0/K_input))*np.asmatrix(beta_scale)[0,0]
	return corrected_betas

def compute_alphas(betas,conc_param):
	term_a = (np.exp(betas)/np.sum(np.exp(betas)))
	alphas = []
	for i,ele in enumerate(conc_param):
		alphas.append(ele*term_a[i])
	return alphas

def compute_pmf(alphas, num_bins):
	pmf = np.zeros(num_bins)
	for bin_num in range(num_bins):
		bin_log_pmf = 0
		bin_log_pmf = bin_log_pmf + loggamma(sum(alphas)) - loggamma(sum(alphas) + 1)
		for bin_iter in range(num_bins):
			if bin_iter == bin_num:
				bin_log_pmf = bin_log_pmf + loggamma(1.0 + alphas[bin_iter]) - loggamma(alphas[bin_iter])
		pmf[bin_num] = np.exp(bin_log_pmf)
	return pmf

def dirichlet_categorical_fit(outlier_counts, inlier_counts, concShape, concRate):
	num_classes = outlier_counts.shape[0]
	num_bins = inlier_counts.shape[1]
	outlier_pmf = np.zeros(outlier_counts.shape)
	inlier_pmf = np.zeros(inlier_counts.shape)
	for class_num in range(num_classes):
		outlier_class_counts = outlier_counts[class_num,:]
		inlier_class_counts = inlier_counts[class_num,:]
		# Put data in dictionary (required input for pystan)
		data_outlier = dict(K=num_bins, y=outlier_class_counts, concShape = concShape,concRate = concRate)
		data_inlier = dict(K=num_bins, y=inlier_class_counts, concShape = concShape,concRate = concRate)
		# Optimize GLM using pystan
		op_outlier = DC_GLM.optimizing(data = data_outlier,verbose=False,iter=5000,seed=1)
		op_inlier = DC_GLM.optimizing(data = data_inlier,verbose=False,iter=5000,seed=1)
		# Extract betas
		beta_outlier = correct_betas(op_outlier['beta_raw'], op_outlier['beta_scale'], num_bins)
		beta_inlier = correct_betas(op_inlier['beta_raw'], op_inlier['beta_scale'], num_bins)
		#compute actual alpha that defines DM
		alphas_outlier = compute_alphas(beta_outlier, op_outlier['conc'])
		alphas_inlier = compute_alphas(beta_inlier, op_inlier['conc'])
		# Compute pmf of this distribution
		#pmf_outlier = compute_pmf(alphas_outlier, num_bins)
		#pmf_inlier = compute_pmf(alphas_inlier, num_bins)
		pmf_outlier = alphas_outlier/sum(alphas_outlier)
		pmf_inlier = alphas_inlier/sum(alphas_inlier)
		outlier_pmf[class_num, :] = pmf_outlier
		inlier_pmf[class_num, :] = pmf_inlier
	ret = dict(phi_inlier=inlier_pmf, phi_outlier=outlier_pmf)
	return ret

