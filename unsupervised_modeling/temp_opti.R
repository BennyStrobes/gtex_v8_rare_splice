source("watershed.R")

initialize_genomic_annotation_variables_v_rand <- function(number_of_features, number_of_dimensions, independent_variables) {
	if (independent_variables == "true") {
		pair_value = 0
		theta_pair = matrix(pair_value,1, choose(number_of_dimensions, 2))
	} else if (independent_variables == "false") {
		pair_value = 4
		theta_pair = matrix(pair_value,1, choose(number_of_dimensions, 2))
	} else if (independent_variables == "false_geno") {
		pair_value = 1e-7
		theta_pair = matrix(0, number_of_features+1, choose(number_of_dimensions, 2))
		theta_pair[1,] <- numeric(choose(number_of_dimensions, 2)) + pair_value
	}

	beta_init = matrix(rnorm((number_of_features+1)*number_of_dimensions,sd=.1),number_of_features+1, number_of_dimensions)
	theta_singleton = beta_init[1,]
	theta = beta_init[2:(number_of_features + 1),]
	# Initialize vector
	x <- c()
	# Add theta_singleton (intercepts)
	x <- c(x, theta_singleton)
	# Add theta for each dimension (betas)
	for (dimension in 1:number_of_dimensions) {
		x <- c(x, theta[, dimension])
	}
	# Add theta_pair (edges between unobserved nodes)
	for (row_number in 1:(dim(theta_pair)[1])) {
		x <- c(x, theta_pair[row_number,])
	}
	return(x)
}


adadelta=function(grad_fxn, log_likelihood_fxn, x, feat, discrete_outliers, posterior, posterior_pairwise, phi, lambda, lambda_pair, lambda_singleton, mu_init, mu_pairwise_init, convergence_thresh, step_size, independent_variables, rho=.9, eps=1e-8) {
  #NOTE: rho=.7, eps=1e-5,bin_size=5000,thresh=.3 is working nicely
  #Note: length(x) == noptim + 2*nintegrate (1 for mean and one for standard deviation)
  #previous gradeint (will be length(x))
  exp_grad_old = numeric(length(x))
  exp_delta_x_old = numeric(length(x))
  progress=c()
  iterations = 1

  number_of_dimensions <- dim(discrete_outliers)[2]
  num_genomic_features <- dim(feat)[2]

  # Get crf coefficients back into inference format
  theta_singleton <- x[1:number_of_dimensions]
  theta <- matrix(0,num_genomic_features,number_of_dimensions)
  for (dimension in 1:number_of_dimensions) {
	theta[,dimension] <- x[(number_of_dimensions + 1 + num_genomic_features*(dimension-1)):(number_of_dimensions + num_genomic_features*(dimension))]
  }
  theta_pair <- matrix(x[(number_of_dimensions + (number_of_dimensions*num_genomic_features) + 1):length(x)], ncol=choose(number_of_dimensions, 2),byrow=TRUE)
  mu_list <- update_marginal_probabilities_vi_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi$inlier_component, phi$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), step_size, convergence_thresh, mu_init, FALSE)
  mu <- mu_list$probability
  mu_pairwise <- mu_list$probability_pairwise





  ##This is literally just adagrad for all (noptim + 2*nintegrate) parameters
  for (iterat in 1:10000) {
    #length(g) == length(x)
    g <- grad_fxn(x, feat=feat, discrete_outliers=discrete_outliers, posterior=posterior, posterior_pairwise=posterior_pairwise, phi=phi, lambda=lambda, lambda_pair=lambda_pair, lambda_singleton=lambda_singleton, mu_init=mu, mu_pairwise_init=mu_pairwise, convergence_thresh=convergence_thresh, step_size=step_size, independent_variables=independent_variables)
    exp_grad_old = rho*exp_grad_old + (1-rho)*g^2
    denominator = sqrt(exp_grad_old + eps)
    numerator = sqrt(exp_delta_x_old + eps)
    delta_x = -(numerator/denominator)*g 
    exp_delta_x_old = rho*exp_delta_x_old + (1-rho)*delta_x^2
    x = x + delta_x
    iterations = iterations + 1
    likelihood <- log_likelihood_fxn(x, feat=feat, discrete_outliers=discrete_outliers, posterior=posterior, posterior_pairwise=posterior_pairwise, phi=phi, lambda=lambda, lambda_pair=lambda_pair, lambda_singleton=lambda_singleton, mu_init=mu_init, mu_pairwise_init=mu_pairwise_init, convergence_thresh=convergence_thresh, step_size=step_size, independent_variables=independent_variables)
  	progress <- c(progress, likelihood)

  	# Get crf coefficients back into inference format
    theta_singleton <- x[1:number_of_dimensions]
    theta <- matrix(0,num_genomic_features,number_of_dimensions)
    for (dimension in 1:number_of_dimensions) {
	   theta[,dimension] <- x[(number_of_dimensions + 1 + num_genomic_features*(dimension-1)):(number_of_dimensions + num_genomic_features*(dimension))]
     }
     theta_pair <- matrix(x[(number_of_dimensions + (number_of_dimensions*num_genomic_features) + 1):length(x)], ncol=choose(number_of_dimensions, 2),byrow=TRUE)
     mu_list <- update_marginal_probabilities_vi_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi$inlier_component, phi$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), step_size, convergence_thresh, mu, FALSE)
     mu <- mu_list$probability
     mu_pairwise <- mu_list$probability_pairwise

  }
  list(x=x,log_prob=progress)
}

adagrad=function(grad_fxn, log_likelihood_fxn, x, feat, discrete_outliers, posterior, posterior_pairwise, phi, lambda, lambda_pair, lambda_singleton, mu_init, mu_pairwise_init, convergence_thresh, step_size, independent_variables, master_stepsize=1.0, eps=1e-6, iterations=10000) {


  number_of_dimensions <- dim(discrete_outliers)[2]
  num_genomic_features <- dim(feat)[2]

  # Get crf coefficients back into inference format
  theta_singleton <- x[1:number_of_dimensions]
  theta <- matrix(0,num_genomic_features,number_of_dimensions)
  for (dimension in 1:number_of_dimensions) {
	theta[,dimension] <- x[(number_of_dimensions + 1 + num_genomic_features*(dimension-1)):(number_of_dimensions + num_genomic_features*(dimension))]
  }
  theta_pair <- matrix(x[(number_of_dimensions + (number_of_dimensions*num_genomic_features) + 1):length(x)], ncol=choose(number_of_dimensions, 2),byrow=TRUE)
  mu_list <- update_marginal_probabilities_vi_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi$inlier_component, phi$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), step_size, convergence_thresh, mu_init, FALSE)
  mu <- mu_list$probability
  mu_pairwise <- mu_list$probability_pairwise


  #Note: length(x) == noptim + 2*nintegrate (1 for mean and one for standard deviation)
  #previous gradeint (will be length(x))
  historical_grad=0
  progress=list()
  ##This is literally just adagrad for all (noptim + 2*nintegrate) parameters
  for (i in 1:iterations) {
    #length(g) == length(x)
    g <- grad_fxn(x, feat=feat, discrete_outliers=discrete_outliers, posterior=posterior, posterior_pairwise=posterior_pairwise, phi=phi, lambda=lambda, lambda_pair=lambda_pair, lambda_singleton=lambda_singleton, mu_init=mu, mu_pairwise_init=mu_pairwise, convergence_thresh=convergence_thresh, step_size=step_size, independent_variables=independent_variables)
    #Note: length(g^2) == length(x). Therefor g^2 is just the dot product
    historical_grad=historical_grad+g^2
    x=x-master_stepsize*g/(eps+sqrt(historical_grad))
    #progress[[i]]=attr(g,"log_prob")#This is the new log posterior probability
    #print(x)
    likelihood <- log_likelihood_fxn(x, feat=feat, discrete_outliers=discrete_outliers, posterior=posterior, posterior_pairwise=posterior_pairwise, phi=phi, lambda=lambda, lambda_pair=lambda_pair, lambda_singleton=lambda_singleton, mu_init=mu_init, mu_pairwise_init=mu_pairwise_init, convergence_thresh=convergence_thresh, step_size=step_size, independent_variables=independent_variables)
  	print(likelihood)


  	# Get crf coefficients back into inference format
    theta_singleton <- x[1:number_of_dimensions]
    theta <- matrix(0,num_genomic_features,number_of_dimensions)
    for (dimension in 1:number_of_dimensions) {
	   theta[,dimension] <- x[(number_of_dimensions + 1 + num_genomic_features*(dimension-1)):(number_of_dimensions + num_genomic_features*(dimension))]
     }
     theta_pair <- matrix(x[(number_of_dimensions + (number_of_dimensions*num_genomic_features) + 1):length(x)], ncol=choose(number_of_dimensions, 2),byrow=TRUE)
     mu_list <- update_marginal_probabilities_vi_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi$inlier_component, phi$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), step_size, convergence_thresh, mu, FALSE)
     mu <- mu_list$probability
     mu_pairwise <- mu_list$probability_pairwise

  }

  list(x=x,log_prob=unlist(progress))
}


grad_desc=function(grad_fxn, log_likelihood_fxn, elbo_fxn, x, feat, discrete_outliers, posterior, posterior_pairwise, phi, lambda, lambda_pair, lambda_singleton, mu_init, mu_pairwise_init, convergence_thresh, step_size, independent_variables, master_stepsize=1, convergence_criteria=.000001) {


  print(lambda_pair)
  number_of_dimensions <- dim(discrete_outliers)[2]
  num_genomic_features <- dim(feat)[2]

  # Get crf coefficients back into inference format
  theta_singleton <- x[1:number_of_dimensions]
  theta <- matrix(0,num_genomic_features,number_of_dimensions)
  for (dimension in 1:number_of_dimensions) {
	theta[,dimension] <- x[(number_of_dimensions + 1 + num_genomic_features*(dimension-1)):(number_of_dimensions + num_genomic_features*(dimension))]
  }
  theta_pair <- matrix(x[(number_of_dimensions + (number_of_dimensions*num_genomic_features) + 1):length(x)], ncol=choose(number_of_dimensions, 2),byrow=TRUE)
  mu_list <- update_marginal_probabilities_vi_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi$inlier_component, phi$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), step_size, convergence_thresh, mu_init, FALSE)
  mu <- mu_list$probability
  mu_pairwise <- mu_list$probability_pairwise

  convergence = FALSE
  iterations = 1
  progress=list()
  prev_likelihood = 1000000000000
  while(convergence==FALSE) {
    g <- grad_fxn(x, feat=feat, discrete_outliers=discrete_outliers, posterior=posterior, posterior_pairwise=posterior_pairwise, phi=phi, lambda=lambda, lambda_pair=lambda_pair, lambda_singleton=lambda_singleton, mu_init=mu, mu_pairwise_init=mu_pairwise, convergence_thresh=convergence_thresh, step_size=step_size, independent_variables=independent_variables)
    x = x-master_stepsize*g
    

  	# Get crf coefficients back into inference format
    theta_singleton <- x[1:number_of_dimensions]
    theta <- matrix(0,num_genomic_features,number_of_dimensions)
    for (dimension in 1:number_of_dimensions) {
	   theta[,dimension] <- x[(number_of_dimensions + 1 + num_genomic_features*(dimension-1)):(number_of_dimensions + num_genomic_features*(dimension))]
     }
     theta_pair <- matrix(x[(number_of_dimensions + (number_of_dimensions*num_genomic_features) + 1):length(x)], ncol=choose(number_of_dimensions, 2),byrow=TRUE)
     mu_list <- update_marginal_probabilities_vi_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi$inlier_component, phi$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), step_size, convergence_thresh, mu, FALSE)
     mu <- mu_list$probability
     mu_pairwise <- mu_list$probability_pairwise




    likelihood <- log_likelihood_fxn(x, feat=feat, discrete_outliers=discrete_outliers, posterior=posterior, posterior_pairwise=posterior_pairwise, phi=phi, lambda=lambda, lambda_pair=lambda_pair, lambda_singleton=lambda_singleton, mu_init=mu, mu_pairwise_init=mu_pairwise, convergence_thresh=convergence_thresh, step_size=step_size, independent_variables=independent_variables)

    print(paste0("Likelihood: ",likelihood))
    print(paste0("ELBO: ", elbo))
    print(as.vector(x))
    if (abs(prev_likelihood - likelihood) < convergence_criteria) {
    	convergence = FALSE
    }
    prev_likelihood = likelihood
    iterations = iterations + 1
    #print(as.vector(x))
  }
 list(x=x,log_prob=unlist(progress))
}

set.seed(5)
lambda = .001
independent_variables="false"

feat <- readRDS("feat.rds")
out <- readRDS("out.rds")
out_pair <- readRDS("pair_out.rds")

number_of_dimensions <- dim(out)[2]
number_of_features <- dim(feat)[2]

gradient_variable_vec <- initialize_genomic_annotation_variables_v_rand(number_of_features, number_of_dimensions, independent_variables)

#if (independent_variables == "true") {
#	pair_value = 0
#	theta_pair = matrix(pair_value,1, choose(number_of_dimensions, 2))
#} else if (independent_variables == "false") {
#	pair_value = 1e-7
	######################
	#pair_value = 2
	######################
#	theta_pair = matrix(pair_value,1, choose(number_of_dimensions, 2))
#} else if (independent_variables == "false_geno") {
#	pair_value = 1e-7
#	theta_pair = matrix(0, number_of_features+1, choose(number_of_dimensions, 2))
#	theta_pair[1,] <- numeric(choose(number_of_dimensions, 2)) + pair_value
#}

# Create GAM parameters section
#beta_init = matrix(0,number_of_features+1, number_of_dimensions)
#################################
# beta_init = matrix(rnorm((number_of_features+1)*number_of_dimensions,sd=.01), number_of_features+1, number_of_dimensions)
#################################
#theta_singleton = beta_init[1,]
#theta = beta_init[2:(number_of_features + 1),]
#gam_parameters = list(theta_pair=theta_pair, theta_singleton=theta_singleton, theta=theta)


phi_placeholder <- initialize_phi(3, number_of_dimensions) 


# Inference using exact inference
#posterior_list <- update_marginal_probabilities_exact_inference_cpp(feat, out+1, gam_parameters$theta_singleton, gam_parameters$theta_pair, gam_parameters$theta, phi_placeholder$inlier_component, phi_placeholder$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), TRUE)
#print(head(posterior_list$probability))

# Inference using vi
mu_init <- matrix(.5, dim(feat)[1], number_of_dimensions)
mu_pairwise_init <- matrix(.5, dim(feat)[1], number_of_dimensions)

step_size = .5
convergence_thresh = 1e-8
#posterior_list <- update_marginal_probabilities_vi_cpp(feat, out+1, gam_parameters$theta_singleton, gam_parameters$theta_pair, gam_parameters$theta, phi_placeholder$inlier_component, phi_placeholder$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), step_size, convergence_thresh, mu_init, TRUE)

#print(head(posterior_list$probability))



# Gradients using exact inference
#num_grad <- grad(compute_exact_crf_likelihood_for_lbfgs, gradient_variable_vec, feat=feat, discrete_outliers=out, posterior=out, posterior_pairwise=out_pair, phi=phi_placeholder, lambda=lambda, lambda_pair=0, lambda_singleton=0, independent_variables=independent_variables)
actual_grad <- compute_exact_crf_gradient_for_lbfgs(gradient_variable_vec, feat=feat, discrete_outliers=out, posterior=out, posterior_pairwise=out_pair, phi=phi_placeholder, lambda=lambda, lambda_pair=0, lambda_singleton=0, independent_variables=independent_variables)


vi_grad <- compute_vi_crf_gradient_for_lbfgs(gradient_variable_vec, feat=feat, discrete_outliers=out, posterior=out, posterior_pairwise=out_pair, phi=phi_placeholder, lambda=lambda, lambda_pair=0, lambda_singleton=0, mu_init=mu_init, mu_pairwise_init=mu_pairwise_init, convergence_thresh=convergence_thresh, step_size=step_size, independent_variables=independent_variables)



# Extract relevent scalers describing data
num_genomic_features <- ncol(feat)
num_samples <- nrow(feat)
number_of_dimensions <- ncol(out)

# Get crf coefficients back into inference format
theta_singleton <- gradient_variable_vec[1:number_of_dimensions]
theta <- matrix(0,num_genomic_features,number_of_dimensions)
for (dimension in 1:number_of_dimensions) {
	theta[,dimension] <- gradient_variable_vec[(number_of_dimensions + 1 + num_genomic_features*(dimension-1)):(number_of_dimensions + num_genomic_features*(dimension))]
}
theta_pair <- matrix(gradient_variable_vec[(number_of_dimensions + (number_of_dimensions*num_genomic_features) + 1):length(gradient_variable_vec)], ncol=choose(number_of_dimensions, 2),byrow=TRUE)



#mu_list <- update_marginal_probabilities_vi_cpp(feat, out, theta_singleton, theta_pair, theta, phi_placeholder$inlier_component, phi_placeholder$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), step_size, convergence_thresh, mu_init, FALSE)
#mu_init <- mu_list$probability
#mu_pairwise_init <- mu_list$probability_pairwise

#num_vi_grad <- grad(compute_vi_crf_likelihood_for_lbfgs, gradient_variable_vec, feat=feat, discrete_outliers=out, posterior=out, posterior_pairwise=out_pair, phi=phi_placeholder, lambda=lambda, lambda_pair=0, lambda_singleton=0, mu_init=mu_init, mu_pairwise_init=mu_pairwise_init, convergence_thresh=convergence_thresh, step_size=step_size, independent_variables=independent_variables)

#lbfgs_output <- lbfgs(compute_exact_crf_likelihood_for_lbfgs, compute_exact_crf_gradient_for_lbfgs, gradient_variable_vec, feat=feat, discrete_outliers=out, posterior=out, posterior_pairwise=out_pair, phi=phi_placeholder, lambda=lambda, lambda_pair=0, lambda_singleton=0, independent_variables=independent_variables,invisible=1)

#print(lbfgs_output$value)
#print(as.vector(lbfgs_output$par))
#print(lbfgs_output$convergence)
#print(lbfgs_output$message)


#lbfgs_output <- lbfgs(compute_vi_crf_likelihood_for_lbfgs, compute_vi_crf_gradient_for_lbfgs, gradient_variable_vec,feat=feat, discrete_outliers=out, posterior=out, posterior_pairwise=out_pair, phi=phi_placeholder, lambda=lambda, lambda_pair=0, lambda_singleton=0,mu_init=mu_init, mu_pairwise_init=mu_pairwise_init, convergence_thresh=convergence_thresh, step_size=step_size, independent_variables=independent_variables,invisible=0)
#print(lbfgs_output$value)
#print(lbfgs_output$par)compute_vi_crf_gradient_for_lbfgs
#print(lbfgs_output$convergence)
#print(lbfgs_output$message)

adadelta_output <- grad_desc(compute_vi_crf_gradient_for_lbfgs, compute_vi_crf_likelihood_for_lbfgs, compute_vi_elbo, gradient_variable_vec, feat, out, out, out_pair, phi_placeholder, lambda, .1, 0, mu_init, mu_pairwise_init, convergence_thresh, step_size, independent_variables)

#print(as.vector(adadelta_output$x))