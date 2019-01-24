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


























#########################################
# Command line arguments
#########################################

pvalue_threshold <- as.numeric(args[1])  # Threshold for calling outliers
input_file <- args[2]  # Watershed input file
stem <- args[3]  # Used in output files as a unique identifier for this run
watershed_run_dir <- args[4]  # Output directory


