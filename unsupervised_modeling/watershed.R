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





load_watershed_data <- function(input_file, number_of_dimensions, pvalue_threshold) {
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



initialize_model_params <- function(num_samples, num_genomic_features, number_of_dimensions, phi_init, theta_pair_init, theta_singleton_init, theta_init, pseudoc, lambda, lambda_singleton, lambda_pair, inference_method, independent_variables) {

	model_params <- list(theta_pair = theta_pair_init, 
						 theta_singleton = theta_singleton_init,
						 theta = theta_init,
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




update_marginal_posterior_probabilities <- function(feat, discrete_outliers, model_params) {
	if (model_params$inference_method == "exact") {
		posterior_list <- update_marginal_probabilities_exact_inference_cpp(feat, discrete_outliers, model_params$theta_singleton, model_params$theta_pair, model_params$theta, model_params$phi$inlier_component, model_params$phi$outlier_component, model_params$number_of_dimensions, choose(model_params$number_of_dimensions, 2), TRUE)
	}
	return(posterior_list)
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

integratedEM <- function(feat, discrete_outliers, phi_init, theta_pair_init, theta_singleton_init, theta_init, pseudoc, lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method, independent_variables, output_root) {
	model_params <- initialize_model_params(dim(feat)[1], dim(feat)[2], number_of_dimensions, phi_init, theta_pair_init, theta_singleton_init, theta_init, pseudoc, lambda, lambda_singleton, lambda_pair, inference_method, independent_variables)


	#################
	# Start loop here
	##################


	for (iter in 1:101) {
		################ E Step
		expected_posteriors <- update_marginal_posterior_probabilities(feat, discrete_outliers, model_params)
		model_params$posterior = expected_posteriors$probability
		model_params$posterior_pairwise = expected_posteriors$probability_pairwise

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
		# saveRDS(model_params, paste0(output_root, "_model_params_iteration_",iter,".rds"))
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







