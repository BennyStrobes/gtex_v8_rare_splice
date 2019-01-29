#include <vector>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;


// NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(2);
NumericMatrix extract_all_binary_combinations(int n) {
	// Initialize output
	int nTemp = (int)pow(2, n) - 1;

	NumericMatrix combo_mat(nTemp + 1, n);

	for (int i = 0; i <= nTemp; i++)
	{
		for (int k = 0; k < n; k++)
		{
			if ((i >> k) & 0x1)
			{
				combo_mat(i,k) = 1;
			}
			else
			{
				combo_mat(i,k) = 0;
			}
		}
	}
	return combo_mat;
}

// NumericMatrix all_binary_combinations_ignoring_one_row = extract_all_binary_combinations_ignoring_one_row(4,3);
NumericMatrix extract_all_binary_combinations_ignoring_one_row(int n, int column_to_ignore) {
	NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(n-1);

	int N = all_binary_combinations_matrix.nrow();
	NumericMatrix combo_mat(N, n);
	for (int row_num = 0; row_num < N; row_num++) {
		int counter = 0;
		for (int column_num = 0; column_num < n; column_num++) {
			if (column_num == column_to_ignore) {
				combo_mat(row_num, column_num) = 1;
			} else {
				combo_mat(row_num, column_num) = all_binary_combinations_matrix(row_num, counter);
				counter += 1;
			}
		}
	}
	return combo_mat;
}

// NumericMatrix all_binary_combinations_ignoring_two_row = extract_all_binary_combinations_ignoring_two_row(2, 1, 0);
NumericMatrix extract_all_binary_combinations_ignoring_two_row(int n, int column_to_ignore1, int column_to_ignore2) {
	if (n == 2) {
		NumericMatrix combo_mat(1,2);
		combo_mat(0,0) = 1;
		combo_mat(0,1) = 1;
		return combo_mat;
	} else {
		NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(n-2);
		int nrow = all_binary_combinations_matrix.nrow();
		NumericMatrix combo_mat(nrow, n);
		for (int row_num = 0; row_num < nrow; row_num++) {
			int counter = 0;
			for (int column_num = 0; column_num < n; column_num++) {
				if (column_num == column_to_ignore1 || column_num == column_to_ignore2) {
					combo_mat(row_num, column_num) = 1;
				} else {
					combo_mat(row_num, column_num) = all_binary_combinations_matrix(row_num, counter);
					counter += 1;
				}
			}
		}
		return combo_mat;
	}
}

double un_normalized_crf_weight(NumericMatrix all_binary_combinations_matrix, int combination_number, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num) {
	// Initialize weight
	double weight = 0;
	for (int dimension=0; dimension < number_of_dimensions; dimension++) {
		weight += all_binary_combinations_matrix(combination_number, dimension)*theta_singleton(dimension);
		// Loop through features
		for (int d = 0; d < feat.ncol(); d++) {
			weight += all_binary_combinations_matrix(combination_number,dimension)*feat(sample_num,d)*theta(d,dimension);
		}
		for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
			if (dimension != dimension2) {
				weight += all_binary_combinations_matrix(combination_number,dimension)*all_binary_combinations_matrix(combination_number,dimension2)*theta_pair(dimension, dimension2);

			}
		}
		if (all_binary_combinations_matrix(combination_number, dimension) == 1) {
			weight += log(phi_outlier(dimension, discrete_outliers(sample_num, dimension) - 1));
		} else {
			weight += log(phi_inlier(dimension, discrete_outliers(sample_num, dimension) - 1));
		}
	}
	return weight;
}

double exact_normalization_constant(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num) {
	NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(number_of_dimensions);
	double val = 0;

	for (int combination_number = 0; combination_number < all_binary_combinations_matrix.nrow(); combination_number++) {
		double un_normalized_wight = un_normalized_crf_weight(all_binary_combinations_matrix, combination_number, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num);
		val += exp(un_normalized_wight);
	}

	//all_binary_combinations_matrix <- extract_all_binary_combinations(number_of_dimensions)
	//val <- 0
	//for (combination_number in 1:nrow(all_binary_combinations_matrix)) {
	//	weight <- 0
	//	zs <- all_binary_combinations_matrix[combination_number,]
	//	weight <- un_normalized_crf_posterior(zs, theta, theta_singleton, theta_pair, phi, feat_vec, discrete_outliers_vec, number_of_dimensions)
	//	val <- val + exp(weight)
	//}
	//return(as.numeric(log(val)))
	return log(val);
}

double exact_probability(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int combination_number,NumericMatrix all_binary_combinations_ignoring_one_row) {
	double prob = exp(un_normalized_crf_weight(all_binary_combinations_ignoring_one_row, combination_number, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num) - normalization_constant);
	return prob;
}

double exact_marginal_probability(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int dimension) {
	NumericMatrix all_binary_combinations_ignoring_one_row = extract_all_binary_combinations_ignoring_one_row(number_of_dimensions, dimension);
	double marginal_prob = 0;
	for (int combination_number = 0; combination_number < all_binary_combinations_ignoring_one_row.nrow(); combination_number++) {
		marginal_prob += exact_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, combination_number, all_binary_combinations_ignoring_one_row);
	}
	return marginal_prob;
}

// [[Rcpp::export]]
NumericMatrix update_marginal_probabilities_exact_inference_cpp(NumericMatrix probabilities, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions) {
	// Rcpp::Rcout << feat(1,1) << std::endl;  
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Compute normalization constant for this sample
		// normalization_constant <- exact_posterior_normalization_constant(feat[n,], discrete_outliers[n,], model_params$theta_singleton, model_params$theta_pair, model_params$theta, model_params$phi, model_params$number_of_dimensions)
		double normalization_constant = exact_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num);
		// Loop through dimensions
		int dimension = 0;
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			probabilities(sample_num, dimension) = exact_marginal_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension);
		}
	}
	return probabilities;
}