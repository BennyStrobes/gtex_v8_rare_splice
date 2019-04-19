args = commandArgs(trailingOnly=TRUE)
source("watershed.R")


initialize_phi<- function(num_bins,dim) {
  phi_outlier <- matrix(1,dim,num_bins)
  phi_inlier <- matrix(1,dim,num_bins)
  phi_inlier[,1] = .8
  phi_inlier[,2] = .1
  phi_inlier[,3] = .1
 

  phi_outlier[,1] = .01
  phi_outlier[,2] = .29
  phi_outlier[,3] = .7


  ####################
  # Total expression
  ####################
  phi_inlier[2,1] = .05
  phi_inlier[2,2] = .9
  phi_inlier[2,3] = .05


  phi_outlier[2,1] = .49
  phi_outlier[2,2] = .02
  phi_outlier[2,3] = .49



  phi_init <- list(inlier_component = phi_inlier, outlier_component = phi_outlier)
  return(phi_init)
}

 #######################################
# Extract ROC curves and precision recall curves for test set (in each dimension seperately) using:
#### 1. Watershed
#### 2. GAM
#### 3. RNA-only
#######################################
compute_roc_across_dimensions <- function(number_of_dimensions, dimension_labels, posterior_prob_test, real_valued_outliers_test1, gam_posteriors, binary_outliers_test2) {
  roc_object_across_dimensions <- list()
    # Loop through dimensions
    for (dimension in 1:number_of_dimensions) {
      # Name of dimension
      dimension_name <- strsplit(dimension_labels[dimension],"_pval")[[1]][1]
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
  return(roc_object_across_dimensions)
}

extract_pairwise_observed_labels <- function(binary_outliers_train) {
  nrow <- dim(binary_outliers_train)[1]
  num_dim <- dim(binary_outliers_train)[2]

  pairwise_mat <- matrix(0, nrow, choose(number_of_dimensions, 2))

  dimension_counter <- 1
  for (dim1 in 1:num_dim) {
    for (dim2 in dim1:num_dim) {
      if (dim1 != dim2) {
        pairwise_mat[,dimension_counter] <- binary_outliers_train[,dim1]*binary_outliers_train[,dim2]
        dimension_counter <- dimension_counter + 1
      }
    }
  }
  return(pairwise_mat)
}

initialize_genomic_annotation_variables <- function(number_of_features, number_of_dimensions, independent_variables) {
  if (independent_variables == "true") {
    pair_value = 0
    theta_pair = matrix(pair_value,1, choose(number_of_dimensions, 2))
  } else if (independent_variables == "false") {
    pair_value = 1e-7
    theta_pair = matrix(pair_value,1, choose(number_of_dimensions, 2))
  } else if (independent_variables == "false_geno") {
    pair_value = 1e-7
    theta_pair = matrix(0, number_of_features+1, choose(number_of_dimensions, 2))
    theta_pair[1,] <- numeric(choose(number_of_dimensions, 2)) + pair_value
  }

  beta_init = matrix(0,number_of_features+1, number_of_dimensions)
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

genomic_annotation_model_cv <- function(feat_train, binary_outliers_train, nfolds, costs, independent_variables) {
  number_of_dimensions <- dim(binary_outliers_train)[2]
  number_of_features <- dim(feat_train)[2]

  gradient_variable_vec <- initialize_genomic_annotation_variables(number_of_features, number_of_dimensions, independent_variables)

  if (independent_variables == "true") {
    pair_value = 0
    theta_pair = matrix(pair_value,1, choose(number_of_dimensions, 2))
  } else if (independent_variables == "false") {
    pair_value = 1e-7
    theta_pair = matrix(pair_value,1, choose(number_of_dimensions, 2))
  } else if (independent_variables == "false_geno") {
    pair_value = 1e-7
    theta_pair = matrix(0, number_of_features+1, choose(number_of_dimensions, 2))
    theta_pair[1,] <- numeric(choose(number_of_dimensions, 2)) + pair_value
  }

  # Create GAM parameters section
  beta_init = matrix(0,number_of_features+1, number_of_dimensions)
  theta_singleton = beta_init[1,]
  theta = beta_init[2:(number_of_features + 1),]
  gam_parameters = list(theta_pair=theta_pair, theta_singleton=theta_singleton, theta=theta)


  phi_placeholder <- initialize_phi(3, number_of_dimensions) 

  pairwise_binary_outliers_train <- extract_pairwise_observed_labels(binary_outliers_train)

  #Randomly shuffle the data
  random_shuffling_indices <- sample(nrow(feat_train))
  feat_train_shuff <- feat_train[random_shuffling_indices,]
  binary_outliers_train_shuff <- binary_outliers_train[random_shuffling_indices,]
  pairwise_binary_outliers_train_shuff <- pairwise_binary_outliers_train[random_shuffling_indices,]

  #Create nfolds equally size folds
  folds <- cut(seq(1,nrow(feat_train_shuff)),breaks=nfolds,labels=FALSE)

  avg_aucs <- c()
  for (cost_iter in 1:length(costs)) {
    lambda <- costs[cost_iter]
    #Perform nfolds-fold cross validation
    aucs <- c()
    for(i in 1:nfolds){
        #Segement your data by fold using the which() function 
        testIndexes <- which(folds==i,arr.ind=TRUE)
        feat_test_fold <- feat_train_shuff[testIndexes,]
        outliers_test_fold <- binary_outliers_train_shuff[testIndexes,]
        pairwise_outliers_test_fold <- pairwise_binary_outliers_train_shuff[testIndexes,]
        feat_train_fold <- feat_train_shuff[-testIndexes,]
        outliers_train_fold <- binary_outliers_train_shuff[-testIndexes,]
        pairwise_outliers_train_fold <- pairwise_binary_outliers_train_shuff[-testIndexes,]

        # num_grad <- grad(compute_exact_crf_likelihood_for_lbfgs, gradient_variable_vec, feat=feat_test_fold, discrete_outliers=outliers_test_fold, posterior=outliers_test_fold, posterior_pairwise=pairwise_outliers_test_fold, phi=phi_placeholder, lambda=lambda, lambda_pair=0, lambda_singleton=0, independent_variables=independent_variables)
        # actual_grad <- compute_exact_crf_gradient_for_lbfgs(gradient_variable_vec, feat=feat_test_fold, discrete_outliers=outliers_test_fold, posterior=outliers_test_fold, posterior_pairwise=pairwise_outliers_test_fold, phi=phi_placeholder, lambda=lambda, lambda_pair=0, lambda_singleton=0, independent_variables=independent_variables)



      lbfgs_output <- lbfgs(compute_exact_crf_likelihood_for_lbfgs, compute_exact_crf_gradient_for_lbfgs, gradient_variable_vec, feat=feat_train_fold, discrete_outliers=outliers_train_fold, posterior=outliers_train_fold, posterior_pairwise=pairwise_outliers_train_fold, phi=phi_placeholder, lambda=lambda, lambda_pair=0, lambda_singleton=0, independent_variables=independent_variables,invisible=1)
      # Check to make sure LBFGS converged OK
      if (lbfgs_output$convergence != 0) {
        print(paste0("LBFGS optimazation on CRF did not converge. It reported convergence error of: ", lbfgs_output$convergence))
        print(lbfgs_output$message)
      }


      # Get optimized crf coefficients back into model_params format
      gam_parameters$theta_singleton <- lbfgs_output$par[1:number_of_dimensions]
      for (dimension in 1:number_of_dimensions) {
        gam_parameters$theta[,dimension] <- lbfgs_output$par[(number_of_dimensions + 1 + ncol(feat_train)*(dimension-1)):(number_of_dimensions + ncol(feat_train)*(dimension))]
      }
      #gam_parameters$theta_pair[1,] <- lbfgs_output$par[(number_of_dimensions + (number_of_dimensions*ncol(feat_train)) + 1):length(lbfgs_output$par)]
      #theta_pair[1,] <- x[(number_of_dimensions + (number_of_dimensions*num_genomic_features) + 1):length(x)]
      gam_parameters$theta_pair <- matrix(lbfgs_output$par[(number_of_dimensions + (number_of_dimensions*ncol(feat_train)) + 1):length(lbfgs_output$par)], ncol=choose(number_of_dimensions, 2),byrow=TRUE)

      # Compute expected value of the CRFs (mu)
      mu_list_test <- update_marginal_probabilities_exact_inference_cpp(feat_test_fold, outliers_test_fold, gam_parameters$theta_singleton, gam_parameters$theta_pair, gam_parameters$theta, phi_placeholder$inlier_component, phi_placeholder$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), FALSE)
      mu_test <- mu_list_test$probability

      # Compute roc
      test_predictions <- as.vector(mu_test)
      test_labels <- as.vector(outliers_test_fold)
      roc_obj <- roc.curve(scores.class0 = test_predictions[test_labels==1], scores.class1 = test_predictions[test_labels==0], curve = T)
      auc <- roc_obj$auc
      aucs <- c(aucs, auc)
    }
    avg_aucs <- c(avg_aucs, mean(aucs))
  }
  # Get best (one with highest avg auc across folds) lambda 
  print(avg_aucs)
  best_index <- which(avg_aucs==max(avg_aucs))
  best_lambda <- costs[best_index]
  # Using best lambda, recompute GAM
  lbfgs_output <- lbfgs(compute_exact_crf_likelihood_for_lbfgs, compute_exact_crf_gradient_for_lbfgs, gradient_variable_vec, feat=feat_train_shuff, discrete_outliers=binary_outliers_train_shuff, posterior=binary_outliers_train_shuff, posterior_pairwise=pairwise_binary_outliers_train_shuff, phi=phi_placeholder, lambda=best_lambda, lambda_pair=0, lambda_singleton=0, independent_variables=independent_variables,invisible=1)
  # Get optimized crf coefficients back into model_params format
  gam_parameters$theta_singleton <- lbfgs_output$par[1:number_of_dimensions]
  for (dimension in 1:number_of_dimensions) {
    gam_parameters$theta[,dimension] <- lbfgs_output$par[(number_of_dimensions + 1 + ncol(feat_train)*(dimension-1)):(number_of_dimensions + ncol(feat_train)*(dimension))]
  }
  #gam_parameters$theta_pair[1,] <- lbfgs_output$par[(number_of_dimensions + (number_of_dimensions*ncol(feat_train)) + 1):length(lbfgs_output$par)]
  gam_parameters$theta_pair <- matrix(lbfgs_output$par[(number_of_dimensions + (number_of_dimensions*ncol(feat_train)) + 1):length(lbfgs_output$par)], ncol=choose(number_of_dimensions, 2),byrow=TRUE)

  return(list(lambda=best_lambda, gam_parameters=gam_parameters))
}



######################
# Command Line args
######################
pvalue_fraction <- as.numeric(args[1])
number_of_dimensions <- as.numeric(args[2])
inference_method <- args[3]
pseudocount <- as.numeric(args[4])
fully_observed_input_file <- args[5]
all_variants_input_file <- args[6]
output_stem <- args[7]

#######################################
## Load in data
#######################################
# Load in data for all variants
# This is the data we are making predictions on
predictions_data <- load_watershed_data(all_variants_input_file, number_of_dimensions, pvalue_fraction)
mean_feat <- apply(predictions_data$feat, 2, mean)
sd_feat <- apply(predictions_data$feat, 2, sd)
predictions_feat <- scale(predictions_data$feat, center=mean_feat, scale=sd_feat) 
predictions_discretized_outliers <- predictions_data$outliers_discrete
predictions_pvalues_outliers <- predictions_data$outlier_pvalues


# Load in fully observed data. Note that pvalue fraction is just here as  dummy variable (as we do not use binary outlier calls for this analysis)
# This is the data we are training on
training_data <- load_watershed_data(fully_observed_input_file, number_of_dimensions, pvalue_fraction)  
training_feat <- scale(training_data$feat, center=mean_feat, scale=sd_feat)
training_discretized_outliers <- training_data$outliers_discrete
training_binary_outliers <- training_data$outliers_binary

#######################################
## Train Watershed model on fully observed training data
#######################################
#watershed_object <- readRDS("/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/watershed_three_class_roc/fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_roc_object.rds")
#watershed_ind_object <- readRDS("/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/unsupervised_modeling/watershed_three_class_roc/fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_roc_object_ind.rds")
phi_init <- initialize_phi(3, number_of_dimensions) 
costs= c(.1, .01, 1e-3, 1e-4)
nfolds <- 2
lambda_singleton <- 0
lambda_pair <- 0

independent_variables="false"
gam_data <- genomic_annotation_model_cv(training_feat, training_binary_outliers, nfolds, costs, independent_variables)
watershed_model <- integratedEM(training_feat, training_discretized_outliers, phi_init, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta, pseudocount, gam_data$lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method, independent_variables)
saveRDS(watershed_model, file=paste0(output_stem, "_watershed_model.rds"))

independent_variables="true"
gam_data <- genomic_annotation_model_cv(training_feat, training_binary_outliers, nfolds, costs, independent_variables)
river_model <- integratedEM(training_feat, training_discretized_outliers, phi_init, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta, pseudocount, gam_data$lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method, independent_variables)
saveRDS(river_model, file=paste0(output_stem, "_river_model.rds"))


#######################################
## Compute posterior probabilities based on fitted models
#######################################
# Watershed model
watershed_posterior_list <- update_marginal_probabilities_exact_inference_cpp(predictions_feat, predictions_discretized_outliers, watershed_model$theta_singleton, watershed_model$theta_pair, watershed_model$theta, watershed_model$phi$inlier_component, watershed_model$phi$outlier_component, watershed_model$number_of_dimensions, choose(watershed_model$number_of_dimensions, 2), TRUE)
watershed_posteriors <- watershed_posterior_list$probability
# Watershed-independent model
river_posterior_list <- update_marginal_probabilities_exact_inference_cpp(predictions_feat, predictions_discretized_outliers, river_model$theta_singleton, river_model$theta_pair, river_model$theta, river_model$phi$inlier_component, river_model$phi$outlier_component, river_model$number_of_dimensions, choose(river_model$number_of_dimensions, 2), TRUE)
river_posteriors <- river_posterior_list$probability

#######################################
## Save predictions to output file
#######################################

posterior_mat <- cbind(rownames(predictions_feat), predictions_pvalues_outliers, river_posteriors, watershed_posteriors)
colnames(posterior_mat) = c("sample_names","splicing_outlier_pvalue", "total_expression_outlier_pvalue", "ase_outlier_pvalue", "splicing_river_posterior", "total_expression_river_posterior", "ase_river_posterior", "splicing_watershed_posterior", "total_expression_watershed_posterior", "ase_watershed_posterior")

write.table(posterior_mat,file=paste0(output_stem,"_posteriors.txt"), sep="\t", quote=FALSE, row.names=FALSE)
