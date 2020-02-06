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
library(cowplot)
library(RColorBrewer)
source("watershed.R")
options(bitmapType = 'cairo', device = 'pdf')




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

make_posterior_predictions_object_vi <- function(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier_component, phi_outlier_component, number_of_dimensions, posterior_prob_test) {
 	prediction_output <- compute_all_exact_posterior_predictions_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier_component, phi_outlier_component, number_of_dimensions)

	predictions_list <- list()
	for (num_sample in 1:nrow(feat)) {
		sample_list <- list()
		for (combination_number in 1:nrow(prediction_output$combination)) {
			prob <- 1
			for (dimension in 1:number_of_dimensions) {
				if (prediction_output$combination[combination_number,dimension] == 1) {
					prob <- prob*posterior_prob_test[num_sample, dimension]
				} else {
					prob <- prob*(1.0-posterior_prob_test[num_sample, dimension])
				}
			}
			sample_list[[paste(prediction_output$combination[combination_number,]+1,collapse=" ")]] = prob

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


# Confusion matrix where rows are actual labels and columns are predicted labels
# Each row (actual label) is normalized by the sum across all columns (predicted label) in that row (actual label)
make_confusion_matrix <- function(predictions_object, binary_outliers_test2) {
	# Initialize confusion matrix to enries of zero and row names and row labels
	num_classes <- length(labels(predictions_object[[1]]))
	confusion_matrix <- matrix(0,num_classes, num_classes)
	rownames(confusion_matrix) = labels(predictions_object[[1]])
	colnames(confusion_matrix) = labels(predictions_object[[1]])
	# Loop through each test sample
	for (sample_num in 1:nrow(binary_outliers_test2)) {
		# Pseudo-gold standard label
		gold_standard <- paste(binary_outliers_test2[sample_num,], collapse=" ")
		gold_index <- which(labels(predictions_object[[sample_num]]) == gold_standard)
		# Get predicted label according to maximum posterior probability
		max_value = -1
		max_index = -1
		for (combo_num in 1:length(predictions_object[[sample_num]])) {
			val = predictions_object[[sample_num]][[combo_num]]
			if (val > max_value) {
				max_value = val
				max_index <- combo_num
			}
		}
		# Add count to confusion matrix
		confusion_matrix[gold_index, max_index] = confusion_matrix[gold_index,max_index] + 1
	}
	#confusion_matrix = confusion_matrix/rowSums(confusion_matrix)
	return(confusion_matrix)
}



# Confusion matrix where rows are actual labels and columns are predicted labels
generate_confusion_matrix <- function(feat_test, discrete_outliers_test1, binary_outliers_test2, watershed_model, inference_method, posterior_prob_test) {
	if (inference_method == "exact") {
		predictions_object <- make_posterior_predictions_object_exact_inference(feat_test, discrete_outliers_test1, watershed_model$theta_singleton, watershed_model$theta_pair, watershed_model$theta, watershed_model$phi$inlier_component, watershed_model$phi$outlier_component, watershed_model$number_of_dimensions)
	} else if (inference_method == "vi" | inference_method == "pseudolikelihood") {
		predictions_object <- make_posterior_predictions_object_vi(feat_test, discrete_outliers_test1, watershed_model$theta_singleton, watershed_model$theta_pair, watershed_model$theta, watershed_model$phi$inlier_component, watershed_model$phi$outlier_component, watershed_model$number_of_dimensions, posterior_prob_test)
	}
	confusion_matrix <- make_confusion_matrix(predictions_object, binary_outliers_test2+1)  	
	return(confusion_matrix)
}

# Helper function to remove NAs from vector
remove_na <- function(x) {
	return(x[!is.na(x)])
}

compute_confidence_interval_on_auprc_with_non_parametric_bootstrap <- function(predictions, labels, num_bootstrap_samples) {
	set.seed(2)
	# Number of samples we have
	num_samples <- length(predictions)
	# Initialize vector to keep track of auprc across bootstrap iterations
	auprc_arr <- c()
	# Compute aucprc across whole vector
	population_auprc <- pr.curve(scores.class0=predictions[labels==1], scores.class1=predictions[labels==0])$auc.integral
	# Loop through bootstrap samples
	for (bootstrap_sample in 1:num_bootstrap_samples) {
		bootstrap_indices = sample(1:num_samples,size=num_samples,replace=TRUE)
		bootstrap_predictions = predictions[bootstrap_indices]
		bootstrap_labels = labels[bootstrap_indices]
		bootstrap_auprc <- pr.curve(scores.class0=bootstrap_predictions[bootstrap_labels==1], scores.class1=bootstrap_predictions[bootstrap_labels==0])$auc.integral
		auprc_arr <- c(auprc_arr, bootstrap_auprc)
	}
	#print(quantile(auprc_arr,.025))
	#print(quantile(auprc_arr, .975))
	return(auprc_arr)
}

 #######################################
# Extract ROC curves and precision recall curves for test set (in each dimension seperately) using:
#### 1. Watershed
#### 2. GAM
#### 3. RNA-only
#######################################
compute_roc_across_dimensions <- function(number_of_dimensions, dimension_labels, posterior_prob_test, real_valued_outliers_test1, gam_posteriors, cadd_scores, binary_outliers_test2) {
	num_non_parametric_bootstrap_samples <- 20000
	roc_object_across_dimensions <- list()
	pos_list <- c()
	neg_list <- c()
  	# Loop through dimensions
  	for (dimension in 1:number_of_dimensions) {
  		# Name of dimension
  		dimension_name <- strsplit(dimension_labels[dimension],"_pval")[[1]][1]
  		# Pseudo gold standard
  		test_outlier_status <- binary_outliers_test2[,dimension]
  		# river predictions
  		# roc_obj <- roc(test_outlier_status, posterior_prob_test[,dimension])
  		roc_obj <- roc.curve(scores.class0 = remove_na(posterior_prob_test[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(posterior_prob_test[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
  		pr_obj <- pr.curve(scores.class0 = remove_na(posterior_prob_test[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(posterior_prob_test[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)

  		watershed_auc_prc_ci <- compute_confidence_interval_on_auprc_with_non_parametric_bootstrap(posterior_prob_test[,dimension], test_outlier_status, num_non_parametric_bootstrap_samples)

  		pos_list <- c(pos_list, remove_na(posterior_prob_test[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]))
  		neg_list <- c(neg_list, remove_na(posterior_prob_test[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]))
  		# Predictions with only RNA
  		#rna_only_roc_obj <- roc(test_outlier_status, real_valued_outliers_test1[,dimension])
  		rna_only_roc_obj <- roc.curve(scores.class0 = remove_na(real_valued_outliers_test1[,dimension][test_outlier_status==1]), scores.class1 = remove_na(real_valued_outliers_test1[,dimension][test_outlier_status==0]), curve = T)
  		rna_only_pr_obj <- pr.curve(scores.class0 = remove_na(real_valued_outliers_test1[,dimension][test_outlier_status==1]), scores.class1 = remove_na(real_valued_outliers_test1[,dimension][test_outlier_status==0]), curve = T)

  		# predictions with only genomic annotations
  		#gam_roc_obj <- roc(test_outlier_status, gam_posteriors[,dimension])
   		gam_roc_obj <- roc.curve(scores.class0 = remove_na(gam_posteriors[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(gam_posteriors[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
   		gam_pr_obj <- pr.curve(scores.class0 = remove_na(gam_posteriors[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(gam_posteriors[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)

   		gam_auc_prc_ci <- compute_confidence_interval_on_auprc_with_non_parametric_bootstrap(gam_posteriors[,dimension], test_outlier_status, num_non_parametric_bootstrap_samples)

  		# predictions with CADD
   		cadd_roc_obj <- roc.curve(scores.class0 = remove_na(cadd_scores[test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(cadd_scores[test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
   		cadd_pr_obj <- pr.curve(scores.class0 = remove_na(cadd_scores[test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(cadd_scores[test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)


		evaROC <-	
		 list(watershed_sens=roc_obj$curve[,2],
              watershed_spec=1-roc_obj$curve[,1],
         	  watershed_auc=roc_obj$auc,
         	  watershed_pr_auc=pr_obj$auc.integral,
         	  watershed_recall=pr_obj$curve[,1],
         	  watershed_precision=pr_obj$curve[,2],
         	  watershed_pr_auc_bootstraps=watershed_auc_prc_ci,
         	  GAM_sens=gam_roc_obj$curve[,2],
              GAM_spec=1-gam_roc_obj$curve[,1],
              GAM_auc=gam_roc_obj$auc,
         	  GAM_pr_auc=gam_pr_obj$auc.integral,
         	  GAM_recall=gam_pr_obj$curve[,1],
         	  GAM_precision=gam_pr_obj$curve[,2],
         	  GAM_pr_auc_bootstraps=gam_auc_prc_ci,
         	  CADD_sens=cadd_roc_obj$curve[,2],
              CADD_spec=1-cadd_roc_obj$curve[,1],
              CADD_auc=cadd_roc_obj$auc,
         	  CADD_pr_auc=cadd_pr_obj$auc.integral,
         	  CADD_recall=cadd_pr_obj$curve[,1],
         	  CADD_precision=cadd_pr_obj$curve[,2],
         	  rna_only_pr_auc=rna_only_pr_obj$auc.integral,
         	  rna_only_recall=rna_only_pr_obj$curve[,1],
         	  rna_only_precision=rna_only_pr_obj$curve[,2],
              rna_only_sens=rna_only_roc_obj$curve[,2],
              rna_only_spec=1-rna_only_roc_obj$curve[,1],
              rna_only_auc=rna_only_roc_obj$auc,
              num_positive_pairs=length(remove_na(posterior_prob_test[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])])),
              num_negative_pairs=length(remove_na(posterior_prob_test[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])])))

		 roc_object_across_dimensions[[dimension]] <- list(name=dimension_name, evaROC=evaROC)
	}

		pr_obj <- pr.curve(scores.class0=pos_list, scores.class1=neg_list, curve = T)

		evaROC <-	
		 list(watershed_pr_auc=pr_obj$auc.integral,
         	  watershed_recall=pr_obj$curve[,1],
         	  watershed_precision=pr_obj$curve[,2])

		roc_object_across_dimensions[[(number_of_dimensions + 1)]] <- list(name="joint", evaROC=evaROC)

	return(roc_object_across_dimensions)
}


#######################################
## Compute absolute risk for high watershed posterior variants
#######################################
absolute_risk_plot <- function(gam_predictions, watershed_predictions, data_input, pvalue, output_file) {
	#######################################
	# Load in all data (training and test)
	#######################################
	N2_pairs <- data_input$N2_pairs	

  	real_valued_outliers_test2 <- abs(rbind(data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),]))

  	watershed_thresholds <- c(.5,.6,.7,.8,.9)

  	thresh <- c()
  	absolute_risks <- c()
  	outlier_type <- c()
  	model_type <- c()

  	for (index in 1:length(watershed_thresholds)) {
  		watershed_threshold <- watershed_thresholds[index]

  		# splice
  		ws_variants <- watershed_predictions[,1] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,1] <= pvalue)/length(real_valued_outliers_test2[ws_variants,1])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "Splicing")
  		model_type <- c(model_type, "watershed")

  		# TE
  		ws_variants <- watershed_predictions[,2] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,2] <= pvalue)/length(real_valued_outliers_test2[ws_variants,2])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "TE")
  		model_type <- c(model_type, "watershed")

  		# ase
  		ws_variants <- watershed_predictions[,3] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,3] <= pvalue)/length(real_valued_outliers_test2[ws_variants,3])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "ASE")
  		model_type <- c(model_type, "watershed")

  		 # splice
  		ws_variants <- gam_predictions[,1] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,1] <= pvalue)/length(real_valued_outliers_test2[ws_variants,1])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "Splicing")
  		model_type <- c(model_type, "GAM")

  		# TE
  		ws_variants <- gam_predictions[,2] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,2] <= pvalue)/length(real_valued_outliers_test2[ws_variants,2])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "TE")
  		model_type <- c(model_type, "GAM")

  		# ase
  		ws_variants <- gam_predictions[,3] > watershed_threshold
  		risk <- sum(real_valued_outliers_test2[ws_variants,3] <= pvalue)/length(real_valued_outliers_test2[ws_variants,3])
  		thresh <- c(thresh, watershed_threshold)
  		absolute_risks <- c(absolute_risks, risk)
  		outlier_type <- c(outlier_type, "ASE")
  		model_type <- c(model_type, "GAM")


  	}
  	df <- data.frame(watershed_threshold=thresh, absolute_risk=absolute_risks, outlier_type=factor(outlier_type, levels=c("ASE","Splicing","TE")), model_type=factor(model_type, levels=c("watershed", "GAM")))


	p_watershed <- ggplot(data=df[as.character(df$model_type)=="watershed",], aes(x=watershed_threshold, y=absolute_risk, fill=outlier_type)) +
	geom_bar(stat="identity", color="black", position=position_dodge())+
  	gtex_v8_figure_theme() + 
  	labs(x="Watershed posterior", y="Absolute risk", fill="") + 
  	scale_fill_manual(values=c("#7F5A83", "#0D324D", "#BFCDE0"))

	p_gam <- ggplot(data=df[as.character(df$model_type)=="GAM",], aes(x=watershed_threshold, y=absolute_risk, fill=outlier_type)) +
	geom_bar(stat="identity", color="black", position=position_dodge())+
  	gtex_v8_figure_theme() + 
  	labs(x="GAM posterior", y="Absolute risk", fill="") + 
  	scale_fill_manual(values=c("#7F5A83", "#0D324D", "#BFCDE0"))
# Use custom colors
#p + scale_fill_manual(values=c('#999999','#E69F00'))
# Use brewer color palettes
#p + scale_fill_brewer(palette="Blues")

	p <- plot_grid(p_gam, p_watershed, ncol=1)

	ggsave(p, file=output_file,width = 7.2,height=6,units="in")
}

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}



roc_analysis <- function(data_input, number_of_dimensions, lambda_costs, pseudoc, inference_method, independent_variables, vi_step_size, vi_threshold, lambda_init) {
	#watershed_data <- readRDS(file_name)
	#watershed_model <- watershed_data$model_params
	#gam_data <- watershed_data$gam_model_params
	#posterior_prob_test <- watershed_data$watershed_predictions
	#gam_test_posteriors <- watershed_data$gam_predictions

	#######################################
	# Load in all data (training and test)
	#######################################
	feat_all <- data_input$feat
	discrete_outliers_all <- data_input$outliers_discrete
	binary_outliers_all <- data_input$outliers_binary
	fraction_binary_outliers_all <- data_input$fraction_outliers_binary
	N2_pairs <- data_input$N2_pairs

	#######################################
	# Extract training data
	#######################################
	feat_train <- feat_all[is.na(N2_pairs),]
	discrete_outliers_train <- discrete_outliers_all[is.na(N2_pairs),]
	binary_outliers_train <-  binary_outliers_all[is.na(N2_pairs),]

	#######################################
	# Extract Test data
	#######################################
  	feat_test <- rbind(feat_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], feat_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	discrete_outliers_test1 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	discrete_outliers_test2 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])
  	binary_outliers_test1 <- rbind(binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  	binary_outliers_test2 <- rbind(fraction_binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], fraction_binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])
  	# Absolute pvalues from test prediction data set (to be used for RNA-only analysis)
  	real_valued_outliers_test1 <- -log10(abs(rbind(data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])) + 1e-7)


	#######################################
	## Standardize Genomic Annotations (features)
	#######################################
	mean_feat <- apply(feat_all, 2, mean)
	sd_feat <- apply(feat_all, 2, sd)
 	feat_all <- scale(feat_all, center=mean_feat, scale=sd_feat)
 	feat_train <- scale(feat_train, center=mean_feat, scale=sd_feat)
 	feat_test <- scale(feat_test, center=mean_feat, scale=sd_feat)


 	#######################################
	## Fit Genomic Annotation Model (GAM)
	#######################################
	nfolds <- 5

	gam_data <- logistic_regression_genomic_annotation_model_cv(feat_train, binary_outliers_train, nfolds, lambda_costs, lambda_init)
	print(paste0(nfolds,"-fold cross validation on GAM yielded optimal lambda of ", gam_data$lambda))

	gam_posterior_test_obj <- update_independent_marginal_probabilities_exact_inference_cpp(feat_test, binary_outliers_test1, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta, matrix(0,2,2), matrix(0,2,2), number_of_dimensions, choose(number_of_dimensions, 2), FALSE)
	gam_test_posteriors <- gam_posterior_test_obj$probability


	#######################################
	### Initialize phi using genomic annotation model
	#######################################
	gam_posterior_train_obj <- update_independent_marginal_probabilities_exact_inference_cpp(feat_train, binary_outliers_test1, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta, matrix(0,2,2), matrix(0,2,2), number_of_dimensions, choose(number_of_dimensions, 2), FALSE)
	gam_train_posteriors <- gam_posterior_train_obj$probability
	phi_init <- map_phi_initialization(discrete_outliers_train, gam_train_posteriors, number_of_dimensions, pseudoc)

 	#######################################
	## Fit Watershed Model (using training data)
	#######################################
	lambda_singleton <- 0
  	lambda_pair <- gam_data$lambda
  	lambda <- gam_data$lambda

  	watershed_model <- integratedEM(feat_train, discrete_outliers_train, phi_init, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta, pseudoc, lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method, independent_variables, vi_step_size, vi_threshold)



 	#######################################
	## Get test data watershed posterior probabilities
	#######################################
 	posterior_info_test <- update_marginal_posterior_probabilities(feat_test, discrete_outliers_test1, watershed_model)
  	posterior_prob_test <- posterior_info_test$probability  # Marginal posteriors
  	posterior_pairwise_prob_test <- posterior_info_test$probability_pairwise  # Pairwise posteriors

   	#######################################
	## Compute confusion matrix from held out test data
	# Confusion matrix where rows are actual labels and columns are predicted labels
	# Each row (actual label) is normalized by the sum across all columns (predicted label) in that row (actual label)	#######################################
 	#######################################
	confusion_matrix <- generate_confusion_matrix(feat_test, discrete_outliers_test1, binary_outliers_test2, watershed_model, inference_method, posterior_prob_test)
 	#######################################
	# Extract ROC curves and precision recall curves for test set (in each dimension seperately) using:
	#### 1. Watershed
	#### 2. GAM
	#### 3. RNA-only
	#######################################
	cadd_scores <- as.vector(feat_test[,40])
	dimension_labels <- colnames(data_input$outliers_binary)
	roc_object_across_dimensions <- compute_roc_across_dimensions(number_of_dimensions, dimension_labels, posterior_prob_test, real_valued_outliers_test1, gam_test_posteriors, cadd_scores, binary_outliers_test2)

	return(list(roc=roc_object_across_dimensions, confusion=confusion_matrix, model_params=watershed_model, gam_model_params=gam_data, watershed_predictions=posterior_prob_test, held_out_labels=binary_outliers_test2, gam_predictions=gam_test_posteriors, cadd_predictions=cadd_scores))

}



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

get_min_num_discrete_outliers <- function(discrete) {
	arr <- c()
	arr <- c(arr, sum(discrete==1))
	arr <- c(arr, sum(discrete==2))
	arr <- c(arr, sum(discrete==3))
	return(min(arr))
}

load_watershed_data_for_gene_threshold_comparison <- function(input_file, standard_input_file, number_of_dimensions, pvalue_fraction, pvalue_threshold) {
	# LOAD IN STANDARD INPUT DATA (JUST TO GET SAME THRESHOLDS).. not to use data
	raw_data_standard <- read.table(standard_input_file, header=TRUE)
	standard_outlier_pvalues <- as.matrix(raw_data_standard[,(ncol(raw_data_standard)-number_of_dimensions):(ncol(raw_data_standard)-1)])

	standard_outliers_discrete <- get_discretized_outliers(standard_outlier_pvalues)

	standard_outliers_binary <- ifelse(abs(standard_outlier_pvalues)<=pvalue_threshold,1,0)

	standard_num_training_instances = sum(is.na(as.character(raw_data_standard[,"N2pair"])))
	standard_n2_pairs = as.character(raw_data_standard[,"N2pair"])
	#standard_num_training_outliers = sum(standard_outliers_binary[is.na(standard_n2_pairs),])
	# sample name as SubjectID:GeneName
	rownames(standard_outlier_pvalues) <- paste(raw_data_standard[,"SubjectID"], ":", raw_data_standard[,"GeneName"],sep="")

	# Load in actual data
	raw_data <- read.table(input_file, header=TRUE)
	# Get genomic features (first 2 columns are line identifiers and last (number_of_dimensions+1) columns are outlier status' and N2 pair
	feat <- raw_data[,3:(ncol(raw_data)-number_of_dimensions-1)]
	# sample name as SubjectID:GeneName
	rownames(feat) <- paste(raw_data[,"SubjectID"], ":", raw_data[,"GeneName"],sep="")
	# Pvalues of outlier status of a particular sample (rows) for a particular outlier type (columns)
	outlier_pvalues <- as.matrix(raw_data[,(ncol(raw_data)-number_of_dimensions):(ncol(raw_data)-1)])
	# sample name as SubjectID:GeneName
	rownames(outlier_pvalues) <- paste(raw_data[,"SubjectID"], ":", raw_data[,"GeneName"],sep="")
	# Convert outlier status into binary random variables
	fraction_outliers_binary <- ifelse(abs(outlier_pvalues)<=.1,1,0) # Strictly for initialization of binary output matrix
	N2_pairs=factor(raw_data[,"N2pair"], levels=unique(raw_data[,"N2pair"]))
	outliers_binary_temp <- ifelse(abs(outlier_pvalues)<=pvalue_threshold,1,0)

	# Convert outlier status into discretized random variables
	outliers_discrete <- get_discretized_outliers(outlier_pvalues)
	num_training_outliers = 0
	standard_num_training_outliers = 0
	pseudoc_arr <- c()
	for (dimension_num in 1:number_of_dimensions) {
		ordered <- sort(abs(standard_outlier_pvalues[,dimension_num]))
		max_val <- ordered[floor(length(ordered)*pvalue_fraction)]
   		print(max_val)
		fraction_outliers_binary[,dimension_num] <- ifelse(abs(outlier_pvalues[,dimension_num])<=max_val,1,0)


		standard_num_outliers_dimension <- get_min_num_discrete_outliers(standard_outliers_discrete[is.na(standard_n2_pairs),dimension_num])
		num_outliers_dimension <- get_min_num_discrete_outliers(outliers_discrete[is.na(as.character(N2_pairs)), dimension_num])
		frac = num_outliers_dimension/standard_num_outliers_dimension
		frac = sum(outliers_binary_temp[is.na(as.character(N2_pairs)),dimension_num])/sum(standard_outliers_binary[is.na(standard_n2_pairs),dimension_num])
		dimension_pseudoc = (frac*30.0)
		pseudoc_arr <- c(pseudoc_arr, dimension_pseudoc)
		num_training_outliers = num_outliers_dimension + num_training_outliers
		standard_num_training_outliers = standard_num_training_outliers + standard_num_outliers_dimension
	}
  	outliers_binary <- ifelse(abs(outlier_pvalues)<=pvalue_threshold,1,0)



  	#num_training_outliers = sum(outliers_binary_temp[is.na(as.character(N2_pairs)),])
	# Extract array of N2 pairs
	num_training_instances = sum(is.na(as.character(N2_pairs)))
	# Put all data into compact data structure
	#scaled_pseudoc = 30.0*num_training_instances/standard_num_training_instances
	# scaled_pseudoc = 30.0*num_training_outliers/standard_num_training_outliers
	scaled_pseudoc = 30.0

	data_input <- list(feat=as.matrix(feat), outlier_pvalues=as.matrix(outlier_pvalues),outliers_binary=as.matrix(outliers_binary), fraction_outliers_binary=as.matrix(fraction_outliers_binary),outliers_discrete=outliers_discrete, N2_pairs=N2_pairs, pseudoc=scaled_pseudoc)
	return(data_input)
}




#########################################
# Command line arguments
#########################################
input_file <- args[1]  # Watershed input file
output_stem <- args[2]  # Stem to save all output files
number_of_dimensions <- as.numeric(args[3])  # Dimensionality of space
n2_pair_pvalue_fraction <- as.numeric(args[4])  # For N2 Pair cross-validation, pick pvalue threshold for each outlier dimension such that this fraction of cases are positive examples
binary_pvalue_threshold <- as.numeric(args[5])  # Pvalue threshold to call binary outliers for genomic annotation model
standard_input_file <- args[6]




#####################
# Parameters
#####################
lambda_costs <- c(.1,.01,1e-3)
vi_step_size=.8
vi_threshold=1e-8
#lambda_init <- NA
lambda_init <- 0.001

set.seed(3)
print(input_file)
#######################################
## Load in data
#######################################
data_input <- load_watershed_data_for_gene_threshold_comparison(input_file, standard_input_file, number_of_dimensions, n2_pair_pvalue_fraction, binary_pvalue_threshold)
saveRDS(data_input, paste0(output_stem,"_data_input.rds"))
pseudoc = data_input$pseudoc
print(pseudoc)

#######################################
## Run models (RIVER and GAM) assuming edges (connections) between dimensions with exact inference
#######################################
independent_variables = "false"
inference_method = "exact"
output_root <- paste0(output_stem,"_inference_", inference_method, "_independent_", independent_variables)
roc_object_exact <- roc_analysis(data_input, number_of_dimensions, lambda_costs, pseudoc, inference_method, independent_variables, vi_step_size, vi_threshold, lambda_init)
saveRDS(roc_object_exact, paste0(output_root, "_roc_object3.rds"))
#roc_object_exact <- readRDS(paste0(output_root, "_roc_object.rds"))

print(roc_object_exact$roc[[1]]$evaROC$watershed_pr_auc)
print(roc_object_exact$roc[[1]]$evaROC$GAM_pr_auc)
print(roc_object_exact$roc[[2]]$evaROC$watershed_pr_auc)
print(roc_object_exact$roc[[2]]$evaROC$GAM_pr_auc)
print(roc_object_exact$roc[[3]]$evaROC$watershed_pr_auc)
print(roc_object_exact$roc[[3]]$evaROC$GAM_pr_auc)

#######################################
## Run models (RIVER and GAM) assuming no edges (connections) between dimensions
#######################################
independent_variables = "true"
inference_method = "exact"
output_root <- paste0(output_stem,"_inference_", inference_method, "_independent_", independent_variables)
roc_object_independent <- roc_analysis(data_input, number_of_dimensions, lambda_costs, pseudoc, inference_method, independent_variables, vi_step_size, vi_threshold, lambda_init)
saveRDS(roc_object_independent, paste0(output_root, "_roc_object3.rds"))
print(roc_object_independent$roc[[1]]$evaROC$watershed_pr_auc)
print(roc_object_independent$roc[[1]]$evaROC$GAM_pr_auc)
print(roc_object_independent$roc[[2]]$evaROC$watershed_pr_auc)
print(roc_object_independent$roc[[2]]$evaROC$GAM_pr_auc)
print(roc_object_independent$roc[[3]]$evaROC$watershed_pr_auc)
print(roc_object_independent$roc[[3]]$evaROC$GAM_pr_auc)








