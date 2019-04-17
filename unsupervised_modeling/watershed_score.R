args = commandArgs(trailingOnly=TRUE)
source("watershed.R")
source("watershed_roc.R")


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
nfolds <- 5
lambda_singleton <- 0
lambda_pair <- 0

independent_variables="false"
gam_data <- genomic_annotation_model_cv(training_feat, training_binary_outliers, nfolds, costs, independent_variables)
watershed_model <- integratedEM(training_feat, training_discretized_outliers, phi_init, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta, pseudoc, gam_data$lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method, independent_variables)
saveRDS(watershed_model, file=paste0(output_stem, "_watershed_model.rds"))

independent_variables="true"
gam_data <- genomic_annotation_model_cv(training_feat, training_binary_outliers, nfolds, costs, independent_variables)
river_model <- integratedEM(training_feat, training_discretized_outliers, phi_init, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta, pseudoc, gam_data$lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method, independent_variables)
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

