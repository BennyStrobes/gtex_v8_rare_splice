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


get_discretized_outliers <- function(outlier_pvalues) {
	# initialize output
	outliers_discretized <- matrix(0,dim(outlier_pvalues)[1], dim(outlier_pvalues)[2])
	for (dimension in 1:ncol(outlier_pvalues)) {
		# Check if it is total expression
		if (min(outlier_pvalues[,dimension]) < 0) {
			under_expression = outlier_pvalues[,dimension] < 0
			log_pvalues = -log10(abs(outlier_pvalues[,dimension]) + 1e-6)
			log_pvalues[under_expression] = log_pvalues[under_expression]*-1
			#discretized <- cut(log_pvalues,breaks=c(-6.01,-4,-2,-1,1,2,4,6.01))
			discretized <- cut(log_pvalues, breaks=c(-6.01,-1,1,6.01))
		} else {
			log_pvalues = -log10(abs(outlier_pvalues[,dimension]) + 1e-6)
			# discretized <- cut(log_pvalues, 7)
			discretized <- cut(log_pvalues, breaks=c(-.01,1,4,6))
		}
		outliers_discretized[,dimension] = as.numeric(discretized)
	}
	colnames(outliers_discretized) = colnames(outlier_pvalues)
	return(outliers_discretized)
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

visualize_confusion_matrix <- function(confusion_matrix, output_file) {
	confusion_matrix = confusion_matrix/rowSums(confusion_matrix)
    melted_corr <- melt(confusion_matrix)

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1, levels=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"))
    melted_corr$X2 <- factor(melted_corr$X2, levels=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"))
    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) #+ scale_fill_gradient(low="grey",high="plum2")

    heatmap <- heatmap + scale_fill_distiller(palette = "Blues", direction=1, limits=c(0,1))

    heatmap <- heatmap + theme(text = element_text(size=12),axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11), axis.text.x = element_text(angle = 0, vjust=.5)) 
    heatmap <- heatmap + labs(x = "Observed Class", y = "Predicted Class",fill="")

    heatmap <- heatmap + scale_x_discrete(breaks=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"),labels=c("0 0 0", "1 0 0", "0 1 0", "0 0 1", "1 0 1","0 1 1", "1 1 0", "1 1 1"))
    heatmap <- heatmap + scale_y_discrete(breaks=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"),labels=c("0 0 0", "1 0 0", "0 1 0", "0 0 1", "1 0 1","0 1 1", "1 1 0", "1 1 1"))

    return(heatmap)
}

visualize_river_and_watershed_confusion_matrices <- function(watershed_confusion_matrix, river_confusion_matrix, output_file) {
	options(bitmapType = 'cairo', device = 'pdf')
	watershed_confusion_plot <- visualize_confusion_matrix(watershed_confusion_matrix)
	river_confusion_plot <- visualize_confusion_matrix(river_confusion_matrix)


	combined_confusion_plots <- plot_grid(river_confusion_plot, watershed_confusion_plot, ncol=1)
	ggsave(combined_confusion_plots, file=output_file,width = 24,height=24, units="cm")


}


# Confusion matrix where rows are actual labels and columns are predicted labels
generate_confusion_matrix <- function(feat_test, discrete_outliers_test1, binary_outliers_test2, watershed_model) {
	predictions_object <- make_posterior_predictions_object_exact_inference(feat_test, discrete_outliers_test1, watershed_model$theta_singleton, watershed_model$theta_pair, watershed_model$theta, watershed_model$phi$inlier_component, watershed_model$phi$outlier_component, watershed_model$number_of_dimensions)
	confusion_matrix <- make_confusion_matrix(predictions_object, binary_outliers_test2+1)  	
	return(confusion_matrix)
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



			lbfgs_output <- lbfgs(compute_exact_crf_likelihood_for_lbfgs, compute_exact_crf_gradient_for_lbfgs, gradient_variable_vec, feat=feat_train_fold, discrete_outliers=outliers_train_fold, posterior=outliers_train_fold, posterior_pairwise=pairwise_outliers_train_fold, phi=phi_placeholder, lambda=lambda, lambda_pair=0, lambda_singleton=0, independent_variables=independent_variables,invisible=0)
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
	lbfgs_output <- lbfgs(compute_exact_crf_likelihood_for_lbfgs, compute_exact_crf_gradient_for_lbfgs, gradient_variable_vec, feat=feat_train_shuff, discrete_outliers=binary_outliers_train_shuff, posterior=binary_outliers_train_shuff, posterior_pairwise=pairwise_binary_outliers_train_shuff, phi=phi_placeholder, lambda=best_lambda, lambda_pair=0, lambda_singleton=0, independent_variables=independent_variables,invisible=0)
	# Get optimized crf coefficients back into model_params format
	gam_parameters$theta_singleton <- lbfgs_output$par[1:number_of_dimensions]
	for (dimension in 1:number_of_dimensions) {
		gam_parameters$theta[,dimension] <- lbfgs_output$par[(number_of_dimensions + 1 + ncol(feat_train)*(dimension-1)):(number_of_dimensions + ncol(feat_train)*(dimension))]
	}
 	#gam_parameters$theta_pair[1,] <- lbfgs_output$par[(number_of_dimensions + (number_of_dimensions*ncol(feat_train)) + 1):length(lbfgs_output$par)]
	gam_parameters$theta_pair <- matrix(lbfgs_output$par[(number_of_dimensions + (number_of_dimensions*ncol(feat_train)) + 1):length(lbfgs_output$par)], ncol=choose(number_of_dimensions, 2),byrow=TRUE)

	return(list(lambda=best_lambda, gam_parameters=gam_parameters))
}


roc_analysis <- function(data_input, number_of_dimensions, phi_init, costs, pseudoc, inference_method, output_root, independent_variables) {
	#######################################
	# Load in all data (training and test)
	#######################################
	feat_all <- data_input$feat
	discrete_outliers_all <- data_input$outliers_discrete
	binary_outliers_all <- data_input$outliers_binary
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
  	binary_outliers_test2 <- rbind(binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])
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
	nfolds <- 2
	gam_data <- genomic_annotation_model_cv(feat_train, binary_outliers_train, nfolds, costs, independent_variables)
	print(paste0(nfolds,"-fold cross validation on GAM yielded optimal lambda of ", gam_data$lambda))
	# Predict on held out test data
	gam_posterior_test <- update_marginal_probabilities_exact_inference_cpp(feat_test, binary_outliers_test1, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta, phi_init$inlier_component, phi_init$outlier_component, number_of_dimensions, choose(number_of_dimensions, 2), FALSE)
	gam_posteriors <- gam_posterior_test$probability



 	#######################################
	## Fit Watershed Model (using training data)
	#######################################
	lambda_singleton <- 0
  	lambda_pair <- 0
  	lambda <- gam_data$lambda
  	watershed_model <- integratedEM(feat_train, discrete_outliers_train, phi_init, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta, pseudoc, lambda, lambda_singleton, lambda_pair, number_of_dimensions, inference_method, independent_variables, output_root)
	saveRDS(watershed_model, paste0(output_root, "_model_params.rds"))
	#watershed_model <- readRDS(paste0(output_root, "_model_params.rds"))


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
	confusion_matrix <- generate_confusion_matrix(feat_test, discrete_outliers_test1, binary_outliers_test2, watershed_model)


 	#######################################
	# Extract ROC curves and precision recall curves for test set (in each dimension seperately) using:
	#### 1. Watershed
	#### 2. GAM
	#### 3. RNA-only
	#######################################
	dimension_labels <- colnames(data_input$outliers_binary)
	roc_object_across_dimensions <- compute_roc_across_dimensions(number_of_dimensions, dimension_labels, posterior_prob_test, real_valued_outliers_test1, gam_posteriors, binary_outliers_test2)



  	return(list(roc=roc_object_across_dimensions, confusion=confusion_matrix))

}

# x is 1-specificity (false positive rate)
# y is sensitivity  (true positive rate)
plot_pr_all_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
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
plot_roc_all_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
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
plot_roc_gam_watershed_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
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

# x is 1-specificity (false positive rate)
# y is sensitivity  (true positive rate)
plot_roc_river_watershed_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
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
		tpr <- c(tpr, dimension_roc_object_ind$evaROC$watershed_sens)
		fpr <- c(fpr, 1-dimension_roc_object_ind$evaROC$watershed_sens)
		outlier_type <- c(outlier_type, rep(dimension_name, length(dimension_roc_object_ind$evaROC$watershed_sens)))
		prediction_type <- c(prediction_type, rep("river", length(dimension_roc_object_ind$evaROC$watershed_sens)))
	}
	df <- data.frame(tpr=tpr, fpr=fpr, outlier_type=factor(outlier_type), prediction_type=factor(prediction_type, levels=c("watershed","river")))
  

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

plot_pr_gam_watershed_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
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

plot_pr_river_watershed_comparison_curve <- function(roc_object, roc_object_ind, number_of_dimensions, output_file) {
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
output_stem <- args[3]  # Stem to save all output files
number_of_dimensions <- as.numeric(args[4])  # Dimensionality of space
inference_method <- args[5]  # Currently only open is "exact"
pseudoc <- as.numeric(args[6])

#####################
# Parameters
#####################
#costs= c(.01, 1e-3, 1e-4, 0)
costs= c(.01, 1e-3)
phi_init <- initialize_phi(3, number_of_dimensions) 



#######################################
## Load in data
#######################################
data_input <- load_watershed_data(input_file, number_of_dimensions, pvalue_threshold)




#######################################
## Run models (RIVER and GAM) assuming edges (connections) between dimensions
#######################################
independent_variables = "false_geno"
output_root <- paste0(output_stem, "_independent_", independent_variables)
roc_object_geno_edge <- roc_analysis(data_input, number_of_dimensions, phi_init, costs, pseudoc, inference_method, output_root, independent_variables)

#######################################
## Run models (RIVER and GAM) assuming edges (connections) between dimensions
#######################################
independent_variables = "false"
output_root <- paste0(output_stem, "_independent_", independent_variables)
roc_object <- roc_analysis(data_input, number_of_dimensions, phi_init, costs, pseudoc, inference_method, output_root, independent_variables)

#######################################
## Run models (RIVER and GAM) assuming no edges (connections) between dimensions
#######################################
independent_variables = "true"
output_root <- paste0(output_stem, "_independent_", independent_variables)
roc_object_ind <- roc_analysis(data_input, number_of_dimensions, phi_init, costs, pseudoc, inference_method, output_root, independent_variables)

saveRDS(roc_object_geno_edge, paste0(output_stem, "roc_object_geno_edge.rds"))
saveRDS(roc_object, paste0(output_stem, "roc_object.rds"))
saveRDS(roc_object_ind, paste0(output_stem, "roc_object_ind.rds"))

# roc_object <- readRDS("roc_object.rds")
# roc_object_ind <- readRDS("roc_object_ind.rds")
# roc and confusion


#######################################
## Visualize Confusion matrix for both RIVER and Watershed
#######################################
#visualize_river_and_watershed_confusion_matrices(roc_object$confusion, roc_object_ind$confusion, paste0(output_stem, "_confusion_heatmap_river_watershed_comparison.pdf"))


#######################################
## Visualize ROC curves for all models (RNA-only, GAM, RIVER, watershed) and all three outlier types (te, splice, ase)
#######################################
#plot_roc_all_comparison_curve(roc_object$roc, roc_object_ind$roc, number_of_dimensions, paste0(output_stem, "_all_comparison_roc.pdf"))

#######################################
## Visualize precision-recall curves for all models (RNA-only, GAM, RIVER, watershed) and all three outlier types (te, splice, ase)
#######################################
#plot_pr_all_comparison_curve(roc_object$roc, roc_object_ind$roc, number_of_dimensions, paste0(output_stem, "_all_comparison_pr.pdf"))



#######################################
## Visualize precision-recall curves for river-watershed comparison and all three outlier types (te, splice, ase)
#######################################
#plot_pr_river_watershed_comparison_curve(roc_object$roc, roc_object_ind$roc, number_of_dimensions, paste0(output_stem, "_watershed_river_comparison_pr.pdf"))


#######################################
## Visualize precision-recall curves for watershed-GAM comparison and all three outlier types (te, splice, ase)
#######################################
#plot_pr_gam_watershed_comparison_curve(roc_object$roc, roc_object_ind$roc, number_of_dimensions,  paste0(output_stem, "_watershed_gam_comparison_pr.pdf"))


#######################################
## Visualize ROC curves for watershed-GAM comparison and all three outlier types (te, splice, ase)
#######################################
#plot_roc_gam_watershed_comparison_curve(roc_object$roc, roc_object_ind$roc, number_of_dimensions, paste0(output_stem, "_watershed_gam_comparison_roc.pdf"))

