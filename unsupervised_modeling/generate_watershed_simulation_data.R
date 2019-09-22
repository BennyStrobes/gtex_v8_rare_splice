args = commandArgs(trailingOnly=TRUE)
source("watershed.R")

simulate_phi<- function(num_bins,dim) {
  phi_outlier <- matrix(1,dim,num_bins)
  phi_inlier <- matrix(1,dim,num_bins)
  phi_inlier[,1] = .99
  phi_inlier[,2] = .005
  phi_inlier[,3] = .005
 

  phi_outlier[,1] = .01
  phi_outlier[,2] = .29
  phi_outlier[,3] = .7


  ####################
  # Total expression
  ####################
  phi_inlier[2,1] = .005
  phi_inlier[2,2] = .99
  phi_inlier[2,3] = .005


  phi_outlier[2,1] = .5
  phi_outlier[2,2] = .01
  phi_outlier[2,3] = .49



  phi_init <- list(inlier_component = phi_inlier, outlier_component = phi_outlier)
  return(phi_init)
}




simulate_model_params <- function(number_of_dimensions, number_of_features) {
	# Simulate edge weights
	theta_pair = matrix(0,1, choose(number_of_dimensions, 2))
	theta_pair[1,1] = 2.0
	theta_pair[1,2] = 3.52
	theta_pair[1,3] = 2.0
	# Simulate intercepts
	beta_init = matrix(0,number_of_features+1, number_of_dimensions)
	theta_singleton = beta_init[1,]
	theta_singleton[1] = -5
	theta_singleton[2] = -5
	theta_singleton[3] = -5
	# Simulate genomic annotation coefficients
	old_model <- readRDS("/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/watershed_three_class_roc/fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_n2_pair_outlier_fraction_.01_binary_pvalue_threshold_.01_pseudocount_30_inference_exact_independent_false_roc_object2.rds")
	theta <- (old_model$model_params$theta[1:number_of_features,])*.1

	phi <- simulate_phi(3,number_of_dimensions)

	model_params <- list(theta_pair = theta_pair, 
						 theta_singleton = theta_singleton,
						 theta = theta,
						 phi=phi)
	return(model_params)
	
}


simulate_latent_variables <- function(model_params, feat_all, number_of_dimensions) {
	num_samples = dim(feat_all)[1]
	discrete_outliers_ph <- matrix(0, num_samples, 3)
	prediction_output = compute_all_exact_crf_posterior_predictions_cpp(feat_all, discrete_outliers_ph, model_params$theta_singleton, model_params$theta_pair, model_params$theta, model_params$phi$inlier_component, model_params$phi$outlier_component, number_of_dimensions)
	lv = matrix(0, num_samples, number_of_dimensions)
	probs = matrix(0, num_samples, 8)

	for (sample_num in 1:num_samples) {
		sample_prob = prediction_output$probability[sample_num,]

		selection = sample(1:length(sample_prob),size=1,prob=sample_prob)
		probs[sample_num,] = sample_prob
		lv[sample_num,] =prediction_output$combination[selection,]
	}
	print(summary(probs))
	print(summary(lv))
	return(lv)
}


categorical_to_pvalue <- function(categorical, dimension_num) {
	pvalue = 10
	if (dimension_num == 2) {
		if (categorical == 1) {
			pvalue = 10^(-runif(1, min = 1.0000000001, max = 6))
			pvalue = pvalue*-1
		} else if (categorical == 2) {
			pvalue = 10^(-runif(1, min = 0, max = 1))
			if (runif(1,min=0,max=1) < .5) {
				pvalue=pvalue*-1
			}
		} else {
			pvalue = 10^(-runif(1, min = 1.0000000001, max = 6))
		}
	} else {
		if (categorical == 1) {
			pvalue = 10^(-runif(1, min = 0, max = .9999999999999))
		} else if (categorical == 2) {
			pvalue = 10^(-runif(1, min = 1, max = 3.999999999))
		} else {
			pvalue = 10^(-runif(1, min =4, max = 6))
		}

	}
	return(pvalue)
}


simulate_expression_pvalues <- function(model_params, lvs) {
	num_samples = dim(lvs)[1]
	number_of_dimensions = dim(lvs)[2]
	expr <- matrix(0, num_samples, number_of_dimensions)
	for (sample_num in 1:num_samples) {
		for (dimension_num in 1:number_of_dimensions) {
			if (lvs[sample_num,dimension_num] == 1) {
				prob <- model_params$phi$outlier_component[dimension_num,]
			} else {
				prob <- model_params$phi$inlier_component[dimension_num,]
			}
			categorical = sample(1:length(prob),size=1,prob=prob)
			expr[sample_num, dimension_num] = categorical_to_pvalue(categorical, dimension_num)
		}
	}
	return(expr)

}


input_file = args[1]
output_file = args[2]

set.seed(1)


data <- load_watershed_data(input_file, 1, .01,.01)

number_of_dimensions <- 3
number_of_features <- dim(data$feat)[2]

# Generate parameters defining watershed
model_params <- simulate_model_params(number_of_dimensions, number_of_features)

# Standardize features
mean_feat <- apply(data$feat, 2, mean)
sd_feat <- apply(data$feat, 2, sd)
feat_all <- scale(data$feat, center=mean_feat, scale=sd_feat)

# Simulate latent variables given genomic annotations and watershed model params

lvs = simulate_latent_variables(model_params, feat_all, number_of_dimensions)


# Simulate expression p-values given the simulated latent variables
expr = simulate_expression_pvalues(model_params, lvs)


print(summary(abs(expr)))
pvalue_fraction <- .01
for (dimension_num in 1:number_of_dimensions) {
	ordered <- sort(abs(expr[,dimension_num]))
	max_val <- ordered[floor(length(ordered)*pvalue_fraction)]
    print(max_val)
}

# Stream input file
stop = FALSE
count = 0
f = file(input_file, "r")

next_line = readLines(f, n = 1)


sink(output_file)
counter = 0
while(!stop) {
    # Parse the line
    data = strsplit(next_line,'\t')[[1]]    

    if (counter == 0) {
    	counter = counter + 1
    	cat(paste(data[1:20],collapse="\t"))
    	cat(paste0("\t", "pheno_1_pvalue", "\t", "pheno_2_pvalue", "\t", "pheno_3_pvalue","\t", "N2pair", "\n"))
    } else {
    	cat(paste(data[1:20],collapse="\t"))
    	cat("\t")
    	cat(paste(expr[counter,], collapse="\t"))
    	cat("\t")
    	cat(paste0(data[22], "\n"))

    	counter = counter + 1
    }


    next_line = readLines(f, n = 1)

    if(length(next_line) == 0) {
        stop = TRUE
        close(f)
    }

}
# close output file handle
sink()