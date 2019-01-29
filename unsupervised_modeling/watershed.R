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
sourceCpp("crf_exact_updates.cpp")






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
	# Convert outlier status into discretized random variables
	outliers_discrete <- ifelse(outlier_pvalues<=pvalue_threshold,2,1)
	# Extract array of N2 pairs
	N2_pairs=factor(raw_data[,"N2pair"], levels=unique(raw_data[,"N2pair"]))
	# Put all data into compact data structure
	data_input <- list(feat=as.matrix(feat), outlier_pvalues=as.matrix(outlier_pvalues),outliers_discrete=as.matrix(outliers_discrete), N2_pairs=N2_pairs)
	return(data_input)
}


initialize_phi<- function(num_bins,dim) {
  phi_outlier <- matrix(1,dim,num_bins)
  phi_inlier <- matrix(1,dim,num_bins)
  phi_inlier[,1] = .7
  phi_inlier[,2] = .3


  phi_outlier[,1] = .05
  phi_outlier[,2] = .95

  phi_init <- list(inlier_component = phi_inlier, outlier_component = phi_outlier)

}

initialize_model_params <- function(num_samples, num_genomic_features, number_of_dimensions, phi_init, beta_init, pseudoc, lambda, lambda_singleton, lambda_pair, inference_method) {
	model_params <- list(theta_pair = matrix(1e-7,number_of_dimensions,number_of_dimensions), 
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
		for (dimension2 in dimension:number_of_dimensions) {
			if (dimension != dimension2) {
			    weight <- weight + zs[dimension]*zs[dimension2]*theta_pair[dimension, dimension2]
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
	for (n in 1:model_params$num_samples) {
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
		model_params <- update_marginal_posterior_probabilities_exact_inference(feat, discrete_outliers, model_params)
		binary_combinations <- update_marginal_probabilities_exact_inference_cpp(model_params$posterior, feat, discrete_outliers, model_params$theta_singleton, model_params$theta_pair, model_params$theta, model_params$phi$inlier_component, model_params$phi$outlier_component, model_params$number_of_dimensions)
	}
	return(model_params)
}

integratedEM <- function(feat, discrete_outliers, phi_init, beta_init, pseudoc, lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method) {
	model_params <- initialize_model_params(dim(feat)[1], dim(feat)[2], number_of_dimensions, phi_init, beta_init, pseudoc, lambda, lambda_singleton, lambda_pair, inference_method)



	saveRDS(model_params, "/home-1/bstrobe1@jhu.edu/scratch/gtex_v8/rare_var/gtex_v8_rare_splice/unsupervised_modeling/model_params.rds")
	saveRDS(feat, "/home-1/bstrobe1@jhu.edu/scratch/gtex_v8/rare_var/gtex_v8_rare_splice/unsupervised_modeling/feat.rds")
	saveRDS(discrete_outliers, "/home-1/bstrobe1@jhu.edu/scratch/gtex_v8/rare_var/gtex_v8_rare_splice/unsupervised_modeling/discrete_outiers.rds")

	print("Start")
	model_params <- update_marginal_posterior_probabilities(feat, discrete_outliers, model_params)
	print("END")
}



integratedEM_old <- function(Feat, Out, lambda_ideal,
                         phi_init, beta_init, num_bins, costs,pseudoc, mixture_component_init,
                         verbose, lambda, lambda_theta_singleton,lambda_theta_pair){
  #  temporarily reduce number of tissues for algorith development purposes
  #Feat <- cbind(matrix(1,dim(Feat)[1],1),Feat)
  # Out <- Out[,1:3]
  #beta_init <- append(1,beta_init,after=1)
  # phi_init$outlier_component <- phi_init$outlier_component[1:3,]
  # phi_init$inlier_component <- phi_init$inlier_component[1:3,]

  # row.has.na <- apply(Out, 1, function(x){any(is.na(x))})
  # Out <- Out[!row.has.na,]
  # Feat <- Feat[!row.has.na,]
  ###############

  model_params <- temp_initialization(Feat, Out, phi_init, beta_init, lambda, lambda_theta_singleton,lambda_theta_pair,mixture_component_init)


  steps <- 1
  maxIter <- 1000  
  converged <- 0
  for (iter in 1:maxIter) {
    if (verbose) {
      cat(' *** STREAM: EM step ',steps,'\n',sep="")
    }

    ## E-step:
    ## Compute expected posterior probabilities
    ##           given current parameters and data
    model_params <- getFuncRvFeat(Feat, model_params)
    model_params <- getFuncRvPosteriors(Out, model_params)
    expected_data_log_like <- compute_expected_complete_log_likelihood(Feat, Out, model_params)
    cat('    Current expected data log probability after E step: ', expected_data_log_like,'\n',sep='')
    if (verbose) {
      cat('E-step: complete', '\n',sep='')
    }
    ## M-step:

    #  Keep track of previous iteration's parameters in order to check for convergence
    phi_old <- model_params$phi
    beta_old <- model_params$beta
    theta_singleton_old <- model_params$theta_singleton
    theta_pair_old <- model_params$theta_pair

    # Maximum Likelihood Estimate (really Map...) of Phi
    model_params <- mle_phi(Out, model_params, pseudoc, num_bins)

    # Maximum Likelihood estimate (really MAP...) of Beta and Theta
    model_params <- mle_beta(model_params, Feat)


    # Compute observed log probability


    print(head(model_params$posterior))
    print(head(model_params$theta_pair))
    print(model_params$theta_singleton)
    print(model_params$beta)
    print(model_params$phi)


    # Print convergence info
    if (verbose) {
      cat('     M-step: norm(phi_inlier difference) = ',
          round(norm(matrix(model_params$phi$inlier_component)-matrix(phi_old$inlier_component)),4),
          ', norm(phi_outlier difference) = ',
          round(norm(matrix(model_params$phi$outlier_component)-matrix(phi_old$outlier_component)),4),
          ', norm(beta difference) = ',
          round(norm(matrix(model_params$beta)-matrix(beta_old)),4),
          ', norm(theta singleton difference) = ',
          round(norm(matrix(model_params$theta_singleton)-matrix(theta_singleton_old)),4),
          ', norm(theta pairs difference) = ',
          round(norm(matrix(model_params$theta_pair)-matrix(theta_pair_old)),4),
          " *** \n\n", sep="")
    }

    expected_data_log_like <- compute_expected_complete_log_likelihood(Feat, Out, model_params)
    cat('    Current expected data log probability after M step: ', expected_data_log_like,'\n',sep='')

    ## Check convergence
    if ((norm(matrix(model_params$beta) - matrix(beta_old)) < 5e-2) &
        (norm(model_params$phi$inlier_component - phi_old$inlier_component) < 5e-2) &
        (norm(model_params$phi$outlier_component - phi_old$outlier_component) < 1.6) &
        (norm(matrix(model_params$theta_singleton) - matrix(theta_singleton_old)) < 5e-2) &
        (norm(matrix(model_params$theta_pair) - matrix(theta_pair_old)) < .08)
        ) {
      converged <- 1
      break
    }
    steps <- steps + 1
  }

  if (converged == 1) {
    cat(" ::: EM iteration is terminated since it converges within a
        predefined tolerance (0.001) ::: \n\n\n",sep="")
  } else if ((converged == 0) && (iter == maxIter)) {
    cat(" ::: EM iteration is terminated since it reaches a
        predefined maximum value (1000) ::: \n\n\n",sep="")
  }

  median_observed_posterior <- calculate_median_observed_posterior(Out, model_params$posterior)
  median_posterior <- apply(model_params$posterior,1,median)
  list(model_params=model_params,
      posteriors=median_observed_posterior)
}



roc_analysis <- function(data_input, number_of_dimensions, phi_init, costs, pseudoc, inference_method, output_root) {
	# Load in all data (training and test)
	feat_all <- data_input$feat
	discrete_outliers_all <- data_input$outliers_discrete
	N2_pairs <- data_input$N2_pairs

	# Extract training data
	feat_train <- feat_all[is.na(N2_pairs),]
	discrete_outliers_train <- discrete_outliers_all[is.na(N2_pairs),]
	# Extract Test data
  	feat_test <- rbind(feat_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], feat_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	discrete_outliers_test1 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	discrete_outliers_test2 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])

	## Standardize Features
	mean_feat <- apply(feat_all, 2, mean)
	sd_feat <- apply(feat_all, 2, sd)
 	feat_all <- scale(feat_all, center=mean_feat, scale=sd_feat)
 	feat_train <- scale(feat_train, center=mean_feat, scale=sd_feat)
 	feat_test <- scale(feat_test, center=mean_feat, scale=sd_feat)

 	# Initialize beta
 	num_features <- dim(feat_all)[2] + 1  # Plus 1 comes from the intercept
 	beta_init <- matrix(0, num_features, number_of_dimensions)
 	# Initialize matrix keeping track of GAM posterior probabilities in test data
 	gam_posteriors <- matrix(0, dim(feat_test)[1], number_of_dimensions)
 	# Initialize coefficient vectors
 	# Loop through outlier types
 	for (dimension in 1:number_of_dimensions) {
 		# Training binary outlier status for dimension #dimension
 		outlier_status <- discrete_outliers_train[,dimension] - 1 # Minus 1 beause currently on {1,2} scale
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
  	emModel <- integratedEM(feat_train, discrete_outliers_train, phi_init, beta_init, pseudoc, logisticCV$lambda.min, lambda_singleton, lambda_pair, number_of_dimensions, inference_method)



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
#output_root <- paste0(watershed_run_dir, stem)
#pseudoc=50
#costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#phi_init <- initialize_phi(2, number_of_dimensions) 


# Load in data
#data_input <- load_data(input_file, number_of_dimensions, pvalue_threshold)


#roc_analysis(data_input, number_of_dimensions, phi_init, costs, pseudoc, inference_method, output_root)

#roc_analysis_driver(input_file, ZscoreThrd, output_root, dimensions, phi_init, num_bins, costs, verbose, pseudoc,lambda,lambda_theta_singleton,lambda_theta_pair)

