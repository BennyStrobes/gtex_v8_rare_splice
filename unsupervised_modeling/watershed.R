args = commandArgs(trailingOnly=TRUE)
library(glmnet)
library(methods)
library(stats)
library(utils)
library(Biobase)
library(pROC)
library(ggplot2)
library(sigmoid)
library(Rcpp)
library(numDeriv)
library(lbfgs)
library(reshape)
library(grid)
library(PRROC)
library(RColorBrewer)
sourceCpp("crf_exact_updates.cpp")


get_discretized_outliers <- function(outlier_pvalues) {
	# initialize output
	outliers_discretized <- matrix(0,dim(outlier_pvalues)[1], dim(outlier_pvalues)[2])
	for (dimension in 1:ncol(outlier_pvalues)) {
		# Check if it is total expression
		if (as.character(colnames(outlier_pvalues)[dimension]) == "total_expression_pvalue") {
			under_expression = outlier_pvalues[,dimension] < 0
			log_pvalues = -log10(abs(outlier_pvalues[,dimension]) + 1e-6)
			log_pvalues[under_expression] = log_pvalues[under_expression]*-1
			discretized <- cut(log_pvalues,breaks=c(-6.01,-4,-2,-1,1,2,4,6.01))
		} else {
			log_pvalues = -log10(abs(outlier_pvalues[,dimension]) + 1e-6)
			discretized <- cut(log_pvalues, 7)
		}
		outliers_discretized[,dimension] = as.numeric(discretized)
	}
	colnames(outliers_discretized) = colnames(outlier_pvalues)
	return(outliers_discretized)
}



load_data <- function(input_file, number_of_dimensions, pvalue_threshold) {
	raw_data <- read.table(input_file, header=TRUE)
	# Get genomic features (first 2 columns are line identifiers and last (number_of_dimensions+1) columns are outlier status' and N2 pair
	feat <- raw_data[,3:(ncol(raw_data)-number_of_dimensions-1)]
	# sample name as SubjectID:GeneName
	rownames(feat) <- paste(raw_data[,"SubjectID"], ":", raw_data[,"GeneName"],sep="")
	# Pvalues of outlier status of a particular sample (rows) for a particular outlier type (columns)
	outlier_pvalues <- raw_data[,(ncol(raw_data)-number_of_dimensions):(ncol(raw_data)-1)]
	# sample name as SubjectID:GeneName
	rownames(outlier_pvalues) <- paste(raw_data[,"SubjectID"], ":", raw_data[,"GeneName"],sep="")
	# Convert outlier status into binary random variables
	outliers_binary <- ifelse(abs(outlier_pvalues)<=pvalue_threshold,1,0)
	outliers_binary[,2] <- ifelse(abs(outlier_pvalues[,2])<=.1,1,0)
	# Convert outlier status into discretized random variables
	outliers_discrete <- get_discretized_outliers(outlier_pvalues)
	# Extract array of N2 pairs
	N2_pairs=factor(raw_data[,"N2pair"], levels=unique(raw_data[,"N2pair"]))
	# Put all data into compact data structure
	data_input <- list(feat=as.matrix(feat), outlier_pvalues=as.matrix(outlier_pvalues),outliers_binary=as.matrix(outliers_binary),outliers_discrete=outliers_discrete, N2_pairs=N2_pairs)
	return(data_input)
}


initialize_phi<- function(num_bins,dim) {
  phi_outlier <- matrix(1,dim,num_bins)
  phi_inlier <- matrix(1,dim,num_bins)
  phi_inlier[,1] = .4
  phi_inlier[,2] = .1
  phi_inlier[,3] = .1
  phi_inlier[,4] = .1
  phi_inlier[,5] = .1
  phi_inlier[,6] = .1
  phi_inlier[,7] = .1

  phi_outlier[,1] = .05
  phi_outlier[,2] = .05
  phi_outlier[,3] = .1
  phi_outlier[,4] = .1
  phi_outlier[,5] = .2
  phi_outlier[,6] = .2
  phi_outlier[,7] = .3

  ####################
  # Total expression
  ####################
  phi_inlier[2,1] = .1
  phi_inlier[2,2] = .1
  phi_inlier[2,3] = .1
  phi_inlier[2,4] = .4
  phi_inlier[2,5] = .1
  phi_inlier[2,6] = .1
  phi_inlier[2,7] = .1

  phi_outlier[2,1] = .25
  phi_outlier[2,2] = .1
  phi_outlier[2,3] = .1
  phi_outlier[2,4] = .1
  phi_outlier[2,5] = .1
  phi_outlier[2,6] = .1
  phi_outlier[2,7] = .25


  phi_init <- list(inlier_component = phi_inlier, outlier_component = phi_outlier)
  return(phi_init)
}

initialize_model_params <- function(num_samples, num_genomic_features, number_of_dimensions, phi_init, beta_init, pseudoc, lambda, lambda_singleton, lambda_pair, inference_method, independent_variables) {
	if (independent_variables == "true") {
		pair_value = 0
	} else {
		pair_value = 1e-7
	}
	model_params <- list(theta_pair = matrix(pair_value,1, choose(number_of_dimensions, 2)), 
						 theta_singleton = beta_init[1,],
						 theta = beta_init[2:(num_genomic_features + 1),],
						 mu = matrix(.5, num_samples, number_of_dimensions),
						 mu_pairwise = matrix(.5, num_samples, choose(number_of_dimensions, 2)),
						 posterior = matrix(.5, num_samples, number_of_dimensions),
						 posterior_pairwise = matrix(.5, num_samples, choose(number_of_dimensions, 2)),
						 num_samples = num_samples,
						 num_genomic_features = num_genomic_features,
						 number_of_dimensions = number_of_dimensions,
						 phi = phi_init,
						 lambda = lambda,
						 lambda_singleton = lambda_singleton,
						 lambda_pair = lambda_pair,
						 pseudoc = pseudoc,
						 observed_data_log_likelihood = 0,
						 independent_variables = independent_variables,
						 inference_method = inference_method)

   return(model_params)
}


extract_all_binary_combinations <-function(number_of_dimensions) {
	l <- list()
	for (i in 1:number_of_dimensions) {
		l[[i]] <- c(0,1)
	}
	return(expand.grid(l))
}

extract_marginal_binary_combinations <- function(number_of_dimensions, dimension) {
	l <- list()
	for (i in 1:number_of_dimensions) {
		if (i != dimension) {
			l[[i]] <- c(0,1)
		} else {
			l[[i]] <- c(1)
		}
	}
	return(expand.grid(l))
}

extract_marginal_pairwise_binary_combinations <- function(number_of_dimensions, dimension1, dimension2) {
	l <- list()
	for (i in 1:number_of_dimensions) {
		if (i != dimension1 & i != dimension2) {
			l[[i]] <- c(0,1)
		} else {
			l[[i]] <- c(1)
		}
	}
	return(expand.grid(l))
}

un_normalized_crf_posterior <- function(zs, theta, theta_singleton, theta_pair, phi, feat_vec, discrete_outliers_vec, number_of_dimensions) {
	weight <- 0
	for (dimension in 1:number_of_dimensions) {
		weight <- weight + zs[dimension]*theta_singleton[dimension]
		weight <- weight + zs[dimension]*sum(theta[,dimension]*feat_vec)
		dimension_counter = 1
		for (dimension2 in dimension:number_of_dimensions) {
			if (dimension != dimension2) {
			    weight <- weight + zs[dimension]*zs[dimension2]*theta_pair[1, dimension_counter]
			    dimension_counter = dimension_counter + 1
			}
		}
		# Check if functional rv or not
		if (zs[dimension] == 1) { # FRV
			weight <- weight + log(phi$outlier_component[dimension, discrete_outliers_vec[dimension]])
		} else { # !FRV
			weight <- weight + log(phi$inlier_component[dimension, discrete_outliers_vec[dimension]])
		}
	}
	return(as.numeric(weight))

}

exact_posterior_normalization_constant <- function(feat_vec, discrete_outliers_vec, theta_singleton, theta_pair, theta, phi, number_of_dimensions) {
	all_binary_combinations_matrix <- extract_all_binary_combinations(number_of_dimensions)
	val <- 0
	for (combination_number in 1:nrow(all_binary_combinations_matrix)) {
		weight <- 0
		zs <- all_binary_combinations_matrix[combination_number,]
		weight <- un_normalized_crf_posterior(zs, theta, theta_singleton, theta_pair, phi, feat_vec, discrete_outliers_vec, number_of_dimensions)
		val <- val + exp(weight)
	}
	return(as.numeric(log(val)))
}

exact_posterior_prob <- function(zs, normalization_constant, feat_vec, discrete_outliers_vec, theta_singleton, theta_pair, theta, phi, number_of_dimensions) {
	return(as.numeric(exp(un_normalized_crf_posterior(zs, theta, theta_singleton, theta_pair, phi, feat_vec, discrete_outliers_vec, number_of_dimensions) - normalization_constant)))
}

exact_marginal_posterior_prob <- function(dimension, normalization_constant, feat_vec, discrete_outliers_vec, theta_singleton, theta_pair, theta, phi, number_of_dimensions) {
	marginal_binary_combinations_matrix <- extract_marginal_binary_combinations(number_of_dimensions, dimension)
	marginal_prob <- 0
	for (combination_number in 1:nrow(marginal_binary_combinations_matrix)) { 
		zs <- marginal_binary_combinations_matrix[combination_number,]
		prob <- exact_posterior_prob(zs, normalization_constant, feat_vec, discrete_outliers_vec, theta_singleton, theta_pair, theta, phi, number_of_dimensions)
		marginal_prob <- marginal_prob + prob
	}
	return(as.numeric(marginal_prob))
}

exact_marginal_pairwise_posterior_prob <- function(dimension1, dimension2, normalization_constant, feat_vec, discrete_outliers_vec, theta_singleton, theta_pair, theta, phi, number_of_dimensions) {
	marginal_binary_combinations_matrix <- extract_marginal_pairwise_binary_combinations(number_of_dimensions, dimension1, dimension2)
	marginal_prob <- 0
	for (combination_number in 1:nrow(marginal_binary_combinations_matrix)) { 
		zs <- marginal_binary_combinations_matrix[combination_number,]
		prob <- exact_posterior_prob(zs, normalization_constant, feat_vec, discrete_outliers_vec, theta_singleton, theta_pair, theta, phi, number_of_dimensions)
		marginal_prob <- marginal_prob + prob
	}
	return(as.numeric(marginal_prob))
}


update_marginal_posterior_probabilities_exact_inference <- function(feat, discrete_outliers, model_params) {
	# Outer loop going through samples
	#for (n in 1:model_params$num_samples) {
	for (n in 1:600) {
		# Compute normalization constant for this sample
		normalization_constant <- exact_posterior_normalization_constant(feat[n,], discrete_outliers[n,], model_params$theta_singleton, model_params$theta_pair, model_params$theta, model_params$phi, model_params$number_of_dimensions)
		# Loop through dimensions
		pairwise_counter <- 1
		for (dimension in 1:model_params$number_of_dimensions) {
			model_params$posterior[n,dimension] <- exact_marginal_posterior_prob(dimension, normalization_constant, feat[n,], discrete_outliers[n,], model_params$theta_singleton, model_params$theta_pair, model_params$theta, model_params$phi, model_params$number_of_dimensions)
			for (dimension2 in dimension:model_params$number_of_dimensions) {
				if (dimension != dimension2) {
					model_params$posterior_pairwise[n, pairwise_counter] <- exact_marginal_pairwise_posterior_prob(dimension, dimension2, normalization_constant, feat[n,], discrete_outliers[n,], model_params$theta_singleton, model_params$theta_pair, model_params$theta, model_params$phi, model_params$number_of_dimensions)
					pairwise_counter <- pairwise_counter + 1
				}
			}
		}

	}
	return(model_params)
}

update_marginal_posterior_probabilities <- function(feat, discrete_outliers, model_params) {
	if (model_params$inference_method == "exact") {
		#model_params <- update_marginal_posterior_probabilities_exact_inference(feat, discrete_outliers, model_params)
		posterior_list <- update_marginal_probabilities_exact_inference_cpp(feat, discrete_outliers, model_params$theta_singleton, model_params$theta_pair, model_params$theta, model_params$phi$inlier_component, model_params$phi$outlier_component, model_params$number_of_dimensions, choose(model_params$number_of_dimensions, 2), TRUE)
		model_params$posterior = posterior_list$probability
		model_params$posterior_pairwise = posterior_list$probability_pairwise
	}
	return(model_params)
}

# Extract gradient variable vector
# First model_params$number_of_dimensions terms are intercepts for each dimension
# Next there are model_params$number_of_dimensions chunks of length $number_of_genomic_features (each chunk is that dimension's beta)
# Next there are model_params$number_of_dimensions choose 2 theta_pairs
extract_gradient_variable_vector <- function(model_params) {
	# Initialize vector
	x <- c()
	# Add theta_singleton (intercepts)
	x <- c(x, model_params$theta_singleton)
	# Add theta for each dimension (betas)
	for (dimension in 1:model_params$number_of_dimensions) {
		x <- c(x, model_params$theta[, dimension])
	}
	# Add theta_pair (edges between unobserved nodes)
	x <- c(x, model_params$theta_pair[1,])
	return(x)
}

# Calculate gradient of crf likelihood (fxn formatted to be used in LBFGS)
compute_exact_crf_gradient_for_lbfgs <- function(x, feat, discrete_outliers, posterior, posterior_pairwise, phi, lambda, lambda_pair, lambda_singleton, independent_variables) {
	# Extract relevent scalers describing data
	num_genomic_features <- ncol(feat)
	num_samples <- nrow(feat)
	number_of_dimensions <- ncol(discrete_outliers)

	# Get crf coefficients back into inference format
	theta_singleton <- x[1:number_of_dimensions]
	theta <- matrix(0,num_genomic_features,number_of_dimensions)
	for (dimension in 1:number_of_dimensions) {
		theta[,dimension] <- x[(number_of_dimensions + 1 + num_genomic_features*(dimension-1)):(number_of_dimensions + num_genomic_features*(dimension))]
	}
	theta_pair <- matrix(0,1, choose(number_of_dimensions, 2))
	theta_pair[1,] <- x[(number_of_dimensions + (number_of_dimensions*num_genomic_features) + 1):length(x)]

	# Compute expected value of the CRFs (mu)
	mu_list <- update_marginal_probabilities_exact_inference_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi$inlier_component, phi$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), FALSE)
	mu <- mu_list$probability
	mu_pairwise <- mu_list$probability_pairwise

	# Gradient of singleton terms (intercepts)
	grad_singleton <- (colSums(posterior) - colSums(mu))*(1/nrow(posterior)) - 2*lambda_singleton*theta_singleton

	# Gradient of theta terms (betas)
	theta_vec <- x[(number_of_dimensions+1):(length(x)-choose(number_of_dimensions, 2))]
	grad_theta <- c()
	for (dimension in 1:number_of_dimensions) {
		temp_grad <- colSums(posterior[,dimension]*feat) - colSums(mu[,dimension]*feat)
		grad_theta <- c(grad_theta, temp_grad)
	}
	grad_theta <- grad_theta*(1/nrow(posterior)) - 2*lambda*theta_vec

	# Gradient of theta pair terms (edges)
	if (independent_variables == "true") {
		grad_pair <- numeric(nrow(posterior_pairwise))
	} else {
		grad_pair <- (colSums(posterior_pairwise) - colSums(mu_pairwise))*(1/nrow(posterior_pairwise)) - 2*lambda_pair*theta_pair[1,]
	}

	# Merge all gradients
	grad <- c(grad_singleton, grad_theta, grad_pair)
	return(-grad)
}

compute_exact_crf_likelihood_for_lbfgs <- function(x, feat, discrete_outliers, posterior, posterior_pairwise, phi, lambda, lambda_pair, lambda_singleton, independent_variables) {
	# Extract relevent scalers describing data
	num_genomic_features <- ncol(feat)
	num_samples <- nrow(feat)
	number_of_dimensions <- ncol(discrete_outliers)

	# Get crf coefficients back into inference format
	theta_singleton <- x[1:number_of_dimensions]
	theta <- matrix(0,num_genomic_features,number_of_dimensions)
	for (dimension in 1:number_of_dimensions) {
		theta[,dimension] <- x[(number_of_dimensions + 1 + num_genomic_features*(dimension-1)):(number_of_dimensions + num_genomic_features*(dimension))]
	}
	theta_pair <- matrix(0,1, choose(number_of_dimensions, 2))
	theta_pair[1,] <- x[(number_of_dimensions + (number_of_dimensions*num_genomic_features) + 1):length(x)]

	# Compute likelihood in cpp function
	log_likelihood <- compute_crf_likelihood_exact_inference_cpp(posterior, posterior_pairwise, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi$inlier_component, phi$outlier_component, number_of_dimensions, lambda, lambda_pair, lambda_singleton)
	return(-log_likelihood)
}


# Compute MAP estimates of the coefficients defining the conditional random field (CRF)
map_crf <- function(feat, discrete_outliers, model_params) {
	# Extract gradient variable vector
	# First model_params$number_of_dimensions terms are intercepts for each dimension
	# Next there are model_params$number_of_dimensions chunks of length $number_of_genomic_features (each chunk is that dimension's beta)
	# Next there are model_params$number_of_dimensions choose 2 theta_pairs
	gradient_variable_vec <- extract_gradient_variable_vector(model_params)


	# Run LBFGS (https://cran.r-project.org/web/packages/lbfgs/lbfgs.pdf) using our gradient and likelihood functions.
	lbfgs_output <- lbfgs(compute_exact_crf_likelihood_for_lbfgs, compute_exact_crf_gradient_for_lbfgs, gradient_variable_vec, feat=feat, discrete_outliers=discrete_outliers, posterior=model_params$posterior, posterior_pairwise=model_params$posterior_pairwise, phi=model_params$phi, lambda=model_params$lambda, lambda_pair=model_params$lambda_pair, lambda_singleton=model_params$lambda_singleton, independent_variables=model_params$independent_variables)
	
	# Check to make sure LBFGS converged OK
	if (lbfgs_output$convergence != 0) {
		print(paste0("LBFGS optimazation on CRF did not converge. It reported convergence error of: ", lbfgs_output$convergence))
		print(lbfgs_output$message)
	}

	# Get optimized crf coefficients back into model_params format
	model_params$theta_singleton <- lbfgs_output$par[1:model_params$number_of_dimensions]
	for (dimension in 1:model_params$number_of_dimensions) {
		model_params$theta[,dimension] <- lbfgs_output$par[(model_params$number_of_dimensions + 1 + ncol(feat)*(dimension-1)):(model_params$number_of_dimensions + ncol(feat)*(dimension))]
	}
 	model_params$theta_pair[1,] <- lbfgs_output$par[(model_params$number_of_dimensions + (model_params$number_of_dimensions*ncol(feat)) + 1):length(lbfgs_output$par)]

	# Calculate gradient of crf likelihood (fxn formatted to be used in LBFGS)
	#calc_grad <- compute_exact_crf_gradient_for_lbfgs(gradient_variable_vec, feat=feat, discrete_outliers=discrete_outliers, posterior=model_params$posterior, posterior_pairwise=model_params$posterior_pairwise, phi=model_params$phi, lambda=model_params$lambda, lambda_pair=model_params$lambda_pair, lambda_singleton=model_params$lambda_singleton)
	# Calculate likelihood of crf likelihood (fxn formatted to be used in LBFGS)
	#likelihood <- compute_exact_crf_likelihood_for_lbfgs(gradient_variable_vec, feat=feat, discrete_outliers=discrete_outliers, posterior=model_params$posterior, posterior_pairwise=model_params$posterior_pairwise, phi=model_params$phi, lambda=model_params$lambda, lambda_pair=model_params$lambda_pair, lambda_singleton=model_params$lambda_singleton)
	# Calculate numerically  derived gradient to make sure they match
	#numerical_grad <- grad(compute_exact_crf_likelihood_for_lbfgs, gradient_variable_vec, feat=feat, discrete_outliers=discrete_outliers, posterior=model_params$posterior, posterior_pairwise=model_params$posterior_pairwise, phi=model_params$phi, lambda=model_params$lambda, lambda_pair=model_params$lambda_pair, lambda_singleton=model_params$lambda_singleton)
	return(model_params)
}

# Compute MAP estimates of the coefficients defined by P(outlier_status| FR)
map_phi <- function(discrete_outliers, model_params) {
	num_bins = 7
	# Initialize output matrices
	phi_outlier <- matrix(1,model_params$number_of_dimensions, num_bins)	
	phi_inlier <- matrix(1,model_params$number_of_dimensions, num_bins)
	# Count number of times we fall into each bin
	for (bin_number in 1:num_bins) {
    	phi_outlier[,bin_number] <- colSums(((discrete_outliers==bin_number)*model_params$posterior),na.rm=TRUE)
    	phi_inlier[,bin_number] <- colSums(((discrete_outliers==bin_number)*(1-model_params$posterior)),na.rm=TRUE)
  	}
	# Add prior
	phi_outlier <- phi_outlier + model_params$pseudoc
	phi_inlier <- phi_inlier + model_params$pseudoc
	# Normalize
	phi_outlier <- phi_outlier/rowSums(phi_outlier)
	phi_inlier <- phi_inlier/rowSums(phi_inlier)

	# Add to model_params
	model_params$phi$outlier_component <- phi_outlier
	model_params$phi$inlier_component <- phi_inlier
	return(model_params)
}

integratedEM <- function(feat, discrete_outliers, phi_init, beta_init, pseudoc, lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method, independent_variables, output_root) {
	model_params <- initialize_model_params(dim(feat)[1], dim(feat)[2], number_of_dimensions, phi_init, beta_init, pseudoc, lambda, lambda_singleton, lambda_pair, inference_method, independent_variables)



	#saveRDS(model_params, "/home-1/bstrobe1@jhu.edu/scratch/gtex_v8/rare_var/gtex_v8_rare_splice/unsupervised_modeling/model_params.rds")
	saveRDS(feat, "/home-1/bstrobe1@jhu.edu/scratch/gtex_v8/rare_var/gtex_v8_rare_splice/unsupervised_modeling/feat.rds")
	saveRDS(discrete_outliers, "/home-1/bstrobe1@jhu.edu/scratch/gtex_v8/rare_var/gtex_v8_rare_splice/unsupervised_modeling/discrete_outiers.rds")


	#################
	# Start loop here
	##################


	for (iter in 1:101) {
		################ E Step
		model_params <- update_marginal_posterior_probabilities(feat, discrete_outliers, model_params)


		################## Compute observed data log likelihood
		observed_data_log_likelihood <- compute_exact_observed_data_log_likelihood_cpp(feat, discrete_outliers, model_params$theta_singleton, model_params$theta_pair, model_params$theta, model_params$phi$inlier_component, model_params$phi$outlier_component, model_params$number_of_dimensions, choose(model_params$number_of_dimensions, 2))
		model_params$observed_data_log_likelihood = observed_data_log_likelihood
		print('########################')
		print(model_params$theta)
		print(model_params$theta_singleton)
		print(model_params$theta_pair)
		print(model_params$phi)
		print(paste0("ITERATION ", iter))
		print(paste0("observed data log likelihood: ", observed_data_log_likelihood))
		saveRDS(model_params, paste0(output_root, "_model_params_iteration_",iter,".rds"))
		#  Keep track of previous iteration's parameters in order to check for convergence
		phi_old <- model_params$phi
		beta_old <- model_params$beta
		theta_singleton_old <- model_params$theta_singleton
		theta_pair_old <- model_params$theta_pair


		################ M Step
		# Compute MAP estimates of the coefficients defining the conditional random field (CRF)
		model_params <- map_crf(feat, discrete_outliers, model_params)
		# Compute MAP estimates of the coefficients defined by P(outlier_status| FR)
		model_params <- map_phi(discrete_outliers, model_params)
	}
	return(model_params)
}


make_posterior_predictions_object_exact_inference <- function(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier_component, phi_outlier_component, number_of_dimensions) {
	prediction_output <- compute_all_exact_posterior_predictions_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier_component, phi_outlier_component, number_of_dimensions)
	predictions_list <- list()
	for (num_sample in 1:nrow(feat)) {
		sample_list <- list()
		for (combination_number in 1:nrow(prediction_output$combination)) {
			sample_list[[paste(prediction_output$combination[combination_number,]+1,collapse=" ")]] = prediction_output$probability[num_sample, combination_number]
		}
		predictions_list[[num_sample]] = sample_list
	}
	return(predictions_list)

}

compute_output_probabilities <- function(discrete_outliers_test2, predictions_object) {
	pvalues <- c()
	for (sample_num in 1:nrow(discrete_outliers_test2)) {
		if (sum(discrete_outliers_test2[sample_num,] == numeric(ncol(discrete_outliers_test2)) + 1) != ncol(discrete_outliers_test2)) {
			prob <- predictions_object[[sample_num]][[paste(discrete_outliers_test2[sample_num,], collapse=" ")]]
			pvalues <- c(pvalues,prob)
		}
	}
	return(pvalues)
}

compute_accuracy <- function(discrete_outliers_test2, predictions_object) {
	correct_count = 0
	total_count = 0
	for (sample_num in 1:nrow(discrete_outliers_test2)) {
		gold_standard <- paste(discrete_outliers_test2[sample_num,], collapse=" ")
		max_value = -1
		for (combo_num in 1:length(predictions_object[[sample_num]])) {
			val = predictions_object[[sample_num]][[combo_num]]
			if (val > max_value) {
				max_value = val
				name = labels(predictions_object[[sample_num]][combo_num])
			}
			if (name == gold_standard) {
				correct_count = correct_count + 1
			}
			total_count = total_count + 1
		}
	}
	return(correct_count/total_count)
}

make_confusion_matrix <- function(predictions_object, binary_outliers_test2) {
	confusion_matrix <- matrix(0,8,8)
	rownames(confusion_matrix) = labels(predictions_object[[1]])
	colnames(confusion_matrix) = labels(predictions_object[[1]])
	for (sample_num in 1:nrow(binary_outliers_test2)) {
		gold_standard <- paste(binary_outliers_test2[sample_num,], collapse=" ")
		gold_index <- which(labels(predictions_object[[sample_num]]) == gold_standard)
		max_value = -1
		max_index = -1
		for (combo_num in 1:length(predictions_object[[sample_num]])) {
			val = predictions_object[[sample_num]][[combo_num]]
			if (val > max_value) {
				max_value = val
				max_index <- combo_num
			}
		}
		confusion_matrix[gold_index, max_index] = confusion_matrix[gold_index,max_index] + 1
	}
	confusion_matrix = confusion_matrix/rowSums(confusion_matrix)
	return(confusion_matrix)
}

visualize_confusion_matrix <- function(correlation_matrix, output_file) {
	myPalette <- colorRampPalette(brewer.pal(9,'Blues'), space = "Lab")
    melted_corr <- melt(correlation_matrix)

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1, levels=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"))
    melted_corr$X2 <- factor(melted_corr$X2, levels=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"))
    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) #+ scale_fill_gradient(low="grey",high="plum2")
    #heatmap <- heatmap + scale_fill_distiller()
    #heatmap <- heatmap + scale_fill_brewer(values = brewer.pal(3,"RdPu"))
    #heatmap <- heatmap + scale_fill_distiller(palette = "Blues", direction=1)
    #heatmap <- heatmap + scale_fill_gradient(high = "#132B43", low = "white")
    heatmap <- heatmap + scale_fill_gradientn(colours = myPalette(10)[1:6], limits = c(0,1), values = c(0, 0.005, 0.05, 0.1, 0.5, 1))
    heatmap <- heatmap + theme(text = element_text(size=12),axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11), axis.text.x = element_text(angle = 0, vjust=.5)) 
    heatmap <- heatmap + labs(x = "Observed Class", y = "Predicted Class",fill="")

    heatmap <- heatmap + scale_x_discrete(breaks=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"),labels=c("0 0 0", "1 0 0", "0 1 0", "0 0 1", "1 0 1","0 1 1", "1 1 0", "1 1 1"))
    heatmap <- heatmap + scale_y_discrete(breaks=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"),labels=c("0 0 0", "1 0 0", "0 1 0", "0 0 1", "1 0 1","0 1 1", "1 1 0", "1 1 1"))

    ggsave(heatmap, file=output_file, width=15, height=11, unit="cm")



}



roc_analysis <- function(data_input, number_of_dimensions, phi_init, costs, pseudoc, inference_method, output_root, independent_variables) {
	# Load in all data (training and test)
	feat_all <- data_input$feat
	discrete_outliers_all <- data_input$outliers_discrete
	binary_outliers_all <- data_input$outliers_binary
	N2_pairs <- data_input$N2_pairs

	# Extract training data
	feat_train <- feat_all[is.na(N2_pairs),]
	discrete_outliers_train <- discrete_outliers_all[is.na(N2_pairs),]
	binary_outliers_train <-  binary_outliers_all[is.na(N2_pairs),]
	# Extract Test data
  	feat_test <- rbind(feat_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], feat_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	discrete_outliers_test1 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	discrete_outliers_test2 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])
  	binary_outliers_test1 <- rbind(binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	binary_outliers_test2 <- rbind(binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])

  	real_valued_outliers_test1 <- -log10(abs(rbind(data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])) + 1e-7)


	## Standardize Features
	mean_feat <- apply(feat_all, 2, mean)
	sd_feat <- apply(feat_all, 2, sd)
 	feat_all <- scale(feat_all, center=mean_feat, scale=sd_feat)
 	feat_train <- scale(feat_train, center=mean_feat, scale=sd_feat)
 	feat_test <- scale(feat_test, center=mean_feat, scale=sd_feat)

 	# Initialize betas
 	num_features <- dim(feat_all)[2] + 1  # Plus 1 comes from the intercept
 	beta_init <- matrix(0, num_features, number_of_dimensions)
 	# Initialize matrix keeping track of GAM posterior probabilities in test data
 	gam_posteriors <- matrix(0, dim(feat_test)[1], number_of_dimensions)
 	# Initialize coefficient vectors
 	# Loop through outlier types
 	for (dimension in 1:number_of_dimensions) {
 		# Training binary outlier status for dimension #dimension
 		outlier_status <- binary_outliers_train[,dimension] 
 		# Train GAM for this outlier type
 		logisticCV <- cv.glmnet(feat_train, as.vector(outlier_status), lambda=costs, family="binomial", alpha=0, nfolds=10)
		# Get beta intercept corresponding to the lambda (l2 penalty) that does best in n-fold cross-validation
		temp_beta <- logisticCV$glmnet.fit$beta[,logisticCV$lambda == logisticCV$lambda.min]
		intercept <- logisticCV$glmnet.fit$a0[logisticCV$lambda == logisticCV$lambda.min]
		# Fill in values into beta_init
		beta_init[1,dimension] <- intercept  # Add intercept
		beta_init[2:num_features, dimension] <- temp_beta  # Add coefficient vector
		print(paste0("Optimal lambda for dimension ", dimension, " is ", logisticCV$lambda.min))

		## Compute a P(FR | G) for Test data
  		gam_post_prob_test <- predict(logisticCV, feat_test, s="lambda.min", type="response")
  		gam_posteriors[,dimension] <- gam_post_prob_test
  		# SAVE FOR LATER
  		#test_outlier_status <- discrete_outliers_test2[,dimension] - 1
  		#GAM.roc <- roc(test_outlier_status, as.numeric(gam_post_prob_test))
  		#print(GAM.roc)
 	}


  	## Train RIVER on training data
  	lambda_singleton <- 0
  	lambda_pair <- 0
  	# lambda <- logisticCV$lambda.min
  	lambda <- 0
  	#emModel <- integratedEM(feat_train, discrete_outliers_train, phi_init, beta_init, pseudoc, lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method, independent_variables, output_root)
	#saveRDS(emModel, paste0(output_root, "_model_params.rds"))
	emModel <- readRDS(paste0(output_root, "_model_params.rds"))


 	# Get posteriors on test data
  	posterior_info_test <- update_marginal_probabilities_exact_inference_cpp(feat_test, discrete_outliers_test1, emModel$theta_singleton, emModel$theta_pair, emModel$theta, emModel$phi$inlier_component, emModel$phi$outlier_component, emModel$number_of_dimensions, choose(emModel$number_of_dimensions, 2), TRUE)
  	posterior_prob_test <- posterior_info_test$probability
  	posterior_pairwise_prob_test <- posterior_info_test$probability_pairwise

  	predictions_object <- make_posterior_predictions_object_exact_inference(feat_test, discrete_outliers_test1, emModel$theta_singleton, emModel$theta_pair, emModel$theta, emModel$phi$inlier_component, emModel$phi$outlier_component, emModel$number_of_dimensions)
	confusion_matrix <- make_confusion_matrix(predictions_object, binary_outliers_test2+1)  	
	visualize_confusion_matrix(confusion_matrix, paste0(output_root, "confusion_heatmap.pdf"))

  	roc_object_across_dimensions <- list()
  	# Loop through dimensions
  	for (dimension in 1:number_of_dimensions) {
  		# Name of dimension
  		dimension_name <- strsplit(colnames(data_input$outliers_binary)[dimension],"_pval")[[1]][1]
  		# Pseudo gold standard
  		test_outlier_status <- binary_outliers_test2[,dimension]
  		# river predictions
  		# roc_obj <- roc(test_outlier_status, posterior_prob_test[,dimension])
  		roc_obj <- roc.curve(scores.class0 = posterior_prob_test[,dimension][test_outlier_status==1], scores.class1 = posterior_prob_test[,dimension][test_outlier_status==0], curve = T)
  		pr_obj <- pr.curve(scores.class0 = posterior_prob_test[,dimension][test_outlier_status==1], scores.class1 = posterior_prob_test[,dimension][test_outlier_status==0], curve = T)
  		# Predictions with only RNA
  		#rna_only_roc_obj <- roc(test_outlier_status, real_valued_outliers_test1[,dimension])
  		rna_only_roc_obj <- roc.curve(scores.class0 = real_valued_outliers_test1[,dimension][test_outlier_status==1], scores.class1 = real_valued_outliers_test1[,dimension][test_outlier_status==0], curve = T)
  		rna_only_pr_obj <- pr.curve(scores.class0 = real_valued_outliers_test1[,dimension][test_outlier_status==1], scores.class1 = real_valued_outliers_test1[,dimension][test_outlier_status==0], curve = T)

  		# predictions with only genomic annotations
  		#gam_roc_obj <- roc(test_outlier_status, gam_posteriors[,dimension])
   		gam_roc_obj <- roc.curve(scores.class0 = gam_posteriors[,dimension][test_outlier_status==1], scores.class1 = gam_posteriors[,dimension][test_outlier_status==0], curve = T)
   		gam_pr_obj <- pr.curve(scores.class0 = gam_posteriors[,dimension][test_outlier_status==1], scores.class1 = gam_posteriors[,dimension][test_outlier_status==0], curve = T)


		evaROC <-	
		 list(watershed_sens=roc_obj$curve[,2],
              watershed_spec=1-roc_obj$curve[,1],
         	  watershed_auc=roc_obj$auc,
         	  watershed_pr_auc=pr_obj$auc.integral,
         	  watershed_recall=pr_obj$curve[,1],
         	  watershed_precision=pr_obj$curve[,2],
         	  GAM_sens=gam_roc_obj$curve[,2],
              GAM_spec=1-gam_roc_obj$curve[,1],
              GAM_auc=gam_roc_obj$auc,
         	  GAM_pr_auc=gam_pr_obj$auc.integral,
         	  GAM_recall=gam_pr_obj$curve[,1],
         	  GAM_precision=gam_pr_obj$curve[,2],
         	  rna_only_pr_auc=rna_only_pr_obj$auc.integral,
         	  rna_only_recall=rna_only_pr_obj$curve[,1],
         	  rna_only_precision=rna_only_pr_obj$curve[,2],
              rna_only_sens=rna_only_roc_obj$curve[,2],
              rna_only_spec=1-rna_only_roc_obj$curve[,1],
              rna_only_auc=rna_only_roc_obj$auc)


		 roc_object_across_dimensions[[dimension]] <- list(name=dimension_name, evaROC=evaROC)

  	}

  	#predictions_object <- make_posterior_predictions_object_exact_inference(feat_test, discrete_outliers_test1, emModel$theta_singleton, emModel$theta_pair, emModel$theta, emModel$phi$inlier_component, emModel$phi$outlier_component, emModel$number_of_dimensions)
  	#predictions_object_ind <- make_posterior_predictions_object_exact_inference(feat_test, discrete_outliers_test1, emModel_ind$theta_singleton, emModel_ind$theta_pair, emModel_ind$theta, emModel_ind$phi$inlier_component, emModel_ind$phi$outlier_component, emModel_ind$number_of_dimensions)

  	#output_probabilities <- compute_output_probabilities(discrete_outliers_test2, predictions_object)
  	#output_probabilities_ind <- compute_output_probabilities(discrete_outliers_test2, predictions_object_ind)

  	#accuracy <- compute_accuracy(discrete_outliers_test2, predictions_object)
  	#accuracy_ind <- compute_accuracy(discrete_outliers_test2, predictions_object_ind)

  	# Get posteriors on test data
  	#posterior_info_test <- update_marginal_probabilities_exact_inference_cpp(feat_test, discrete_outliers_test1, emModel$theta_singleton, emModel$theta_pair, emModel$theta, emModel$phi$inlier_component, emModel$phi$outlier_component, emModel$number_of_dimensions, choose(emModel$number_of_dimensions, 2), TRUE)
  	#posterior_prob_test <- posterior_info_test$probability
  	#posterior_pairwise_prob_test <- posterior_info_test$probability_pairwise

	# Get posteriors on test data
  	#posterior_info_test_ind <- update_marginal_probabilities_exact_inference_cpp(feat_test, discrete_outliers_test1, emModel_ind$theta_singleton, emModel_ind$theta_pair, emModel_ind$theta, emModel_ind$phi$inlier_component, emModel_ind$phi$outlier_component, emModel_ind$number_of_dimensions, choose(emModel_ind$number_of_dimensions, 2), TRUE)
  	#posterior_prob_test_ind <- posterior_info_test_ind$probability
  	#posterior_pairwise_prob_test_ind <- posterior_info_test_ind$probability_pairwise

 	# Get posteriors on test data
  	#posterior_info_test <- update_marginal_probabilities_exact_inference_cpp(feat_test, discrete_outliers_test1, emModel$theta_singleton, emModel$theta_pair, emModel$theta, emModel$phi$inlier_component, emModel$phi$outlier_component, emModel$number_of_dimensions, choose(emModel$number_of_dimensions, 2), TRUE)
  	#posterior_prob_test <- posterior_info_test$probability
  	#posterior_pairwise_prob_test <- posterior_info_test$probability_pairwise

  	#for (dimension in 1:number_of_dimensions) {
  	#	test_outlier_status <- binary_outliers_test2[,dimension]
  	#	roc_obj <- roc(test_outlier_status, posterior_prob_test[,dimension])
  	#	print(colnames(data_input$outliers_binary)[dimension])
  	#	print(roc_obj$auc)
  	#	rna_only_roc_obj <- roc(test_outlier_status, real_valued_outliers_test1[,dimension])
  	#	print(rna_only_roc_obj$auc)
  	#}

  	# Loop through dimensions
  	#dimension <- 1
  	#test_outlier_status <- binary_outliers_test2[,dimension]
  	#splice_roc <- roc(test_outlier_status, posterior_prob_test[,dimension])
  	#splice_roc_ind <- roc(test_outlier_status, posterior_prob_test_ind[,dimension])

  	 # Loop through dimensions
  	#dimension <- 2
  	#test_outlier_status <- binary_outliers_test2[,dimension] - 1
  	#te_roc <- roc(test_outlier_status, posterior_prob_test[,dimension])
  	#te_roc_ind <- roc(test_outlier_status, posterior_prob_test_ind[,dimension])
  	return(roc_object_across_dimensions)

}

# x is 1-specificity (false positive rate)
# y is sensitivity  (true positive rate)
plot_pr_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	precision <- c()
	recall <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		precision <- c(precision, dimension_roc_object$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_precision)))
		# independent
		precision <- c(precision, dimension_roc_object_ind$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object_ind$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("river", length(dimension_roc_object_ind$evaROC$watershed_precision)))
		# rna only
		precision <- c(precision, dimension_roc_object$evaROC$rna_only_precision)
		recall <- c(recall, dimension_roc_object$evaROC$rna_only_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$rna_only_precision)))
		prediction_type <- c(prediction_type, rep("rna only", length(dimension_roc_object$evaROC$rna_only_precision)))
		# GAM
		precision <- c(precision, dimension_roc_object$evaROC$GAM_precision)
		recall <- c(recall, dimension_roc_object$evaROC$GAM_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$GAM_precision)))
		prediction_type <- c(prediction_type, rep("GAM", length(dimension_roc_object$evaROC$GAM_precision)))
	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type))
  

  	plotter <- ggplot(data=df, aes(x=recall, y=precision, colour=prediction_type)) + geom_line() + facet_wrap( ~ outlier_type, ncol=3) +
                labs(x="Recall", y="Precision", colour="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10))

	ggsave(plotter, file=output_file,width = 26,height=11,units="cm")
}



# x is 1-specificity (false positive rate)
# y is sensitivity  (true positive rate)
plot_roc_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	tpr <- c()
	fpr <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		tpr <- c(tpr, dimension_roc_object$evaROC$watershed_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$watershed_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_sens)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_sens)))
		# independent
		tpr <- c(tpr, dimension_roc_object_ind$evaROC$watershed_sens)
		fpr <- c(fpr, 1-dimension_roc_object_ind$evaROC$watershed_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_sens)))
		prediction_type <- c(prediction_type, rep("river", length(dimension_roc_object_ind$evaROC$watershed_sens)))
		# rna only
		tpr <- c(tpr, dimension_roc_object$evaROC$rna_only_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$rna_only_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$rna_only_sens)))
		prediction_type <- c(prediction_type, rep("rna only", length(dimension_roc_object$evaROC$rna_only_sens)))
		# GAM
		tpr <- c(tpr, dimension_roc_object$evaROC$GAM_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$GAM_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$GAM_sens)))
		prediction_type <- c(prediction_type, rep("GAM", length(dimension_roc_object$evaROC$GAM_sens)))
	}
	df <- data.frame(tpr=tpr, fpr=fpr, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type))
  

  	plotter <- ggplot(data=df, aes(x=fpr, y=tpr, colour=prediction_type)) + geom_line() + facet_wrap( ~ outlier_type, ncol=3) +
  				geom_abline(slope=1) +
                labs(x="False positive rate", y="True positive rate", colour="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=10),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10))

	ggsave(plotter, file=output_file,width = 26,height=11,units="cm")
}



# x is 1-specificity (false positive rate)
# y is sensitivity  (true positive rate)
plot_roc_comparison_curve2 <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	tpr <- c()
	fpr <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		tpr <- c(tpr, dimension_roc_object$evaROC$watershed_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$watershed_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_sens)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_sens)))
		# GAM
		tpr <- c(tpr, dimension_roc_object$evaROC$GAM_sens)
		fpr <- c(fpr, 1-dimension_roc_object$evaROC$GAM_spec)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$GAM_sens)))
		prediction_type <- c(prediction_type, rep("GAM", length(dimension_roc_object$evaROC$GAM_sens)))
	}
	df <- data.frame(tpr=tpr, fpr=fpr, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("watershed","GAM")))
  

  	plotter <- ggplot(data=df, aes(x=fpr, y=tpr, colour=outlier_type, linetype=prediction_type)) + geom_line() +
  				geom_abline(slope=1) +
                labs(x="False positive rate", y="True positive rate", colour="",linetype="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="right") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	ggsave(plotter, file=output_file,width = 15,height=11,units="cm")
}

plot_pr_comparison_curve2 <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	precision <- c()
	recall <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		precision <- c(precision, dimension_roc_object$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_precision)))

		# GAM
		precision <- c(precision, dimension_roc_object$evaROC$GAM_precision)
		recall <- c(recall, dimension_roc_object$evaROC$GAM_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$GAM_precision)))
		prediction_type <- c(prediction_type, rep("GAM", length(dimension_roc_object$evaROC$GAM_precision)))
	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("watershed","GAM")))
  

  	plotter <- ggplot(data=df, aes(x=recall, y=precision, colour=outlier_type, linetype=prediction_type)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="",linetype="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="right") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	ggsave(plotter, file=output_file,width = 19,height=11,units="cm")
}

plot_pr_comparison_curve3 <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
	precision <- c()
	recall <- c()
	outlier_type <- c()
	prediction_type <- c()
	for (dimension in 1:number_of_dimensions) {
		dimension_roc_object <- roc_object[[dimension]]
		dimension_name <- dimension_roc_object$name
		dimension_roc_object_ind <- roc_object_ind[[dimension]]
		# Tied watershed
		precision <- c(precision, dimension_roc_object$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("watershed", length(dimension_roc_object$evaROC$watershed_precision)))

		# Indepdent
		precision <- c(precision, dimension_roc_object_ind$evaROC$watershed_precision)
		recall <- c(recall, dimension_roc_object_ind$evaROC$watershed_recall)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_precision)))
		prediction_type <- c(prediction_type, rep("river", length(dimension_roc_object_ind$evaROC$watershed_precision)))
	}
	df <- data.frame(precision, recall, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("watershed","river")))
  

  	plotter <- ggplot(data=df, aes(x=recall, y=precision, colour=outlier_type, linetype=prediction_type)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="",linetype="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="right") +
                theme(panel.spacing = unit(2, "lines")) +
                theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

	ggsave(plotter, file=output_file,width = 19,height=11,units="cm")
}












#########################################
# Command line arguments
#########################################

pvalue_threshold <- as.numeric(args[1])  # Threshold for calling outliers
input_file <- args[2]  # Watershed input file
stem <- args[3]  # Used in output files as a unique identifier for this run
watershed_run_dir <- args[4]  # Output directory
number_of_dimensions <- as.numeric(args[5])  # Dimensionality of space
inference_method <- args[6]  # Currently only open is "exact"

#####################
# Parameters
#####################
pseudoc=20
costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)
phi_init <- initialize_phi(7, number_of_dimensions) 



# Load in data
data_input <- load_data(input_file, number_of_dimensions, pvalue_threshold)

independent_variables = "false"
output_root <- paste0(watershed_run_dir, stem, "_independent_", independent_variables)
roc_object <- roc_analysis(data_input, number_of_dimensions, phi_init, costs, pseudoc, inference_method, output_root, independent_variables)


independent_variables = "true"
output_root <- paste0(watershed_run_dir, stem, "_independent_", independent_variables)
roc_object_ind <- roc_analysis(data_input, number_of_dimensions, phi_init, costs, pseudoc, inference_method, output_root, independent_variables)

saveRDS(roc_object, "roc_object.rds")
saveRDS(roc_object_ind, "roc_object_ind.rds")

output_file <- paste0(watershed_run_dir, stem, "_comparison_roc.pdf")
plot_roc_comparison_curve(roc_object, roc_object_ind, number_of_dimensions, output_file)


output_file <- paste0(watershed_run_dir, stem, "_comparison_pr.pdf")
plot_pr_comparison_curve(roc_object, roc_object_ind, number_of_dimensions, output_file)

output_file <- paste0(watershed_run_dir, stem, "_watershed_gam_comparison_roc.pdf")
plot_roc_comparison_curve2(roc_object, roc_object_ind, number_of_dimensions, output_file)

output_file <- paste0(watershed_run_dir, stem, "_watershed_gam_comparison_pr.pdf")
plot_pr_comparison_curve2(roc_object, roc_object_ind, number_of_dimensions, output_file)

output_file <- paste0(watershed_run_dir, stem, "_watershed_river_comparison_pr.pdf")
plot_pr_comparison_curve3(roc_object, roc_object_ind, number_of_dimensions, output_file)
