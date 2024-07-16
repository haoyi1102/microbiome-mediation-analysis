rm(list = ls())

library(MASS)
library(stats)
library(graphics)
library(energy)
library(matrixStats)
library(devtools)

source("modima.R")
source("MedOmniTest.R")
# install_github("chanw0/SparseMCMM")
library(SparseMCMM)
# install.packages("ccmm")
library(ccmm)
# install_github("yijuanhu/LDM")
library(LDM)
# install_github("mkoslovsky/MicroBVS")
library(MicroBVS)
# devtools::install_github("quranwu/MedZIM")
library(MedZIM)
library(microHIMA)

set.seed(1234)

# Define the directory where the results are stored
output_dir <- "simulation_data_4"

# List all .rds files in the directory
rds_files <- list.files(output_dir, pattern = "*.rds", full.names = TRUE)

# Load all .rds files into a list
simulation_results <- lapply(rds_files, readRDS)

# Function to preprocess data
preprocess_data <- function(result) {
  Y_vector <- as.numeric(result$Y)
  T_vector <- as.numeric(result$T)
  M_a_matrix <- as.matrix(t(result$absolute_M))
  colnames(M_a_matrix) <- paste0("taxon_", 1:ncol(M_a_matrix))
  
  M_matrix <- M_a_matrix / rowSums(M_a_matrix)
  M_nz_matrix <- apply(M_matrix, 1, function(row) {
    min_value <- min(row[row > 0])
    pseudo_count <- min_value / 100
    row[row == 0] <- pseudo_count
    return(row)
  })
  M_nz_matrix <- t(M_nz_matrix)
  
  list(Y_vector = Y_vector, T_vector = T_vector, M_nz_matrix = M_nz_matrix)
}

# Function to calculate metrics
calculate_metrics <- function(detected_indices, true_indices) {
  true_positives <- length(intersect(detected_indices, true_indices))
  false_positives <- length(setdiff(detected_indices, true_indices))
  false_negatives <- length(setdiff(true_indices, detected_indices))
  
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  list(precision = precision, recall = recall, f1_score = f1_score,
       true_positives = true_positives, false_positives = false_positives, false_negatives = false_negatives)
}

# Initialize metrics for CCMM, LDM, and HIMA methods
metrics_ccmm <- list(precision = 0, recall = 0, f1_score = 0, true_positives = 0, false_positives = 0, false_negatives = 0)
metrics_ldm <- list(precision = 0, recall = 0, f1_score = 0, true_positives = 0, false_positives = 0, false_negatives = 0)
metrics_hima <- list(precision = 0, recall = 0, f1_score = 0, true_positives = 0, false_positives = 0, false_negatives = 0)

# for (i in 1:length(simulation_results)) {
# Loop through each result and apply CCMM, LDM, and HIMA methods
for (i in 1:10) {
  # i = 3
  result <- simulation_results[[i]]
  preprocessed_data <- preprocess_data(result)
  
  Y_vector <- preprocessed_data$Y_vector
  T_vector <- preprocessed_data$T_vector
  M_nz_matrix <- preprocessed_data$M_nz_matrix
  
  spiked_features <- as.numeric(gsub("Feature", "", result$feature$feature_spiked))
  
  # Apply CCMM method
  result_ccmm <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
                      sig.level = 0.05, tol = 1e-06, max.iter = 5000)
  significant_indices_ccmm <- as.numeric(which(result_ccmm$IDE.CIs[1, ] > 0 | result_ccmm$IDE.CIs[2, ] < 0))
  
  # Apply LDM method
  data <- data.frame(Y = Y_vector, T = T_vector)
  result_ldm <- ldm(
    formula = M_nz_matrix ~ T + Y,
    data = data,
    seed = 1234,
    test.mediation = TRUE
  )
  detected_otus_ldm <- as.numeric(gsub("taxon_", "", result_ldm$med.detected.otu.omni))
  
  # Apply HIMA method
  result_hima <- tryCatch({
    mhima(exposure = T_vector, covariates = NULL, otu.com = M_nz_matrix, outcome = Y_vector)
  }, error = function(e) {
    NULL
  })
  
  detected_otus_hima <- if (!is.null(result_hima)) result_hima$ID else integer(0)
  
  # Calculate metrics for CCMM
  metrics_ccmm_i <- calculate_metrics(significant_indices_ccmm, spiked_features)
  metrics_ccmm$true_positives <- metrics_ccmm$true_positives + metrics_ccmm_i$true_positives
  metrics_ccmm$false_positives <- metrics_ccmm$false_positives + metrics_ccmm_i$false_positives
  metrics_ccmm$false_negatives <- metrics_ccmm$false_negatives + metrics_ccmm_i$false_negatives
  
  # Calculate metrics for LDM
  metrics_ldm_i <- calculate_metrics(detected_otus_ldm, spiked_features)
  metrics_ldm$true_positives <- metrics_ldm$true_positives + metrics_ldm_i$true_positives
  metrics_ldm$false_positives <- metrics_ldm$false_positives + metrics_ldm_i$false_positives
  metrics_ldm$false_negatives <- metrics_ldm$false_negatives + metrics_ldm_i$false_negatives
  
  # Calculate metrics for HIMA
  metrics_hima_i <- calculate_metrics(detected_otus_hima, spiked_features)
  metrics_hima$true_positives <- metrics_hima$true_positives + metrics_hima_i$true_positives
  metrics_hima$false_positives <- metrics_hima$false_positives + metrics_hima_i$false_positives
  metrics_hima$false_negatives <- metrics_hima$false_negatives + metrics_hima_i$false_negatives
}

# Finalize metrics calculations
metrics_ccmm$precision <- metrics_ccmm$true_positives / (metrics_ccmm$true_positives + metrics_ccmm$false_positives)
metrics_ccmm$recall <- metrics_ccmm$true_positives / (metrics_ccmm$true_positives + metrics_ccmm$false_negatives)
metrics_ccmm$f1_score <- 2 * (metrics_ccmm$precision * metrics_ccmm$recall) / (metrics_ccmm$precision + metrics_ccmm$recall)

metrics_ldm$precision <- metrics_ldm$true_positives / (metrics_ldm$true_positives + metrics_ldm$false_positives)
metrics_ldm$recall <- metrics_ldm$true_positives / (metrics_ldm$true_positives + metrics_ldm$false_negatives)
metrics_ldm$f1_score <- 2 * (metrics_ldm$precision * metrics_ldm$recall) / (metrics_ldm$precision + metrics_ldm$recall)

metrics_hima$precision <- metrics_hima$true_positives / (metrics_hima$true_positives + metrics_hima$false_positives)
metrics_hima$recall <- metrics_hima$true_positives / (metrics_hima$true_positives + metrics_hima$false_negatives)
metrics_hima$f1_score <- 2 * (metrics_hima$precision * metrics_hima$recall) / (metrics_hima$precision + metrics_hima$recall)


# Print the results for CCMM
cat("CCMM Results:\n")
cat("True Positives:", metrics_ccmm$true_positives, "\n")
cat("False Positives:", metrics_ccmm$false_positives, "\n")
cat("False Negatives:", metrics_ccmm$false_negatives, "\n")
cat("Precision:", metrics_ccmm$precision, "\n")
cat("Recall:", metrics_ccmm$recall, "\n")
cat("F1 Score:", metrics_ccmm$f1_score, "\n\n")

# Print the results for LDM
cat("LDM Results:\n")
cat("True Positives:", metrics_ldm$true_positives, "\n")
cat("False Positives:", metrics_ldm$false_positives, "\n")
cat("False Negatives:", metrics_ldm$false_negatives, "\n")
cat("Precision:", metrics_ldm$precision, "\n")
cat("Recall:", metrics_ldm$recall, "\n")
cat("F1 Score:", metrics_ldm$f1_score, "\n")

# Print the results for HIMA
cat("HIMA Results:\n")
cat("True Positives:", metrics_hima$true_positives, "\n")
cat("False Positives:", metrics_hima$false_positives, "\n")
cat("False Negatives:", metrics_hima$false_negatives, "\n")
cat("Precision:", metrics_hima$precision, "\n")
cat("Recall:", metrics_hima$recall, "\n")
cat("F1 Score:", metrics_hima$f1_score, "\n")
# 
# Print the results for CCMM
# > cat("CCMM Results:\n")
# CCMM Results:
#   > cat("True Positives:", metrics_ccmm$true_positives, "\n")
# True Positives: 14 
# > cat("False Positives:", metrics_ccmm$false_positives, "\n")
# False Positives: 7 
# > cat("False Negatives:", metrics_ccmm$false_negatives, "\n")
# False Negatives: 46 
# > cat("Precision:", metrics_ccmm$precision, "\n")
# Precision: 0.6666667 
# > cat("Recall:", metrics_ccmm$recall, "\n")
# Recall: 0.2333333 
# > cat("F1 Score:", metrics_ccmm$f1_score, "\n\n")
# F1 Score: 0.345679 
# 
# > 
#   > # Print the results for LDM
#   > cat("LDM Results:\n")
# LDM Results:
#   > cat("True Positives:", metrics_ldm$true_positives, "\n")
# True Positives: 3 
# > cat("False Positives:", metrics_ldm$false_positives, "\n")
# False Positives: 4 
# > cat("False Negatives:", metrics_ldm$false_negatives, "\n")
# False Negatives: 57 
# > cat("Precision:", metrics_ldm$precision, "\n")
# Precision: 0.4285714 
# > cat("Recall:", metrics_ldm$recall, "\n")
# Recall: 0.05 
# > cat("F1 Score:", metrics_ldm$f1_score, "\n")
# F1 Score: 0.08955224 
# > 
#   > # Print the results for HIMA
#   > cat("HIMA Results:\n")
# HIMA Results:
#   > cat("True Positives:", metrics_hima$true_positives, "\n")
# True Positives: 1 
# > cat("False Positives:", metrics_hima$false_positives, "\n")
# False Positives: 0 
# > cat("False Negatives:", metrics_hima$false_negatives, "\n")
# False Negatives: 59 
# > cat("Precision:", metrics_hima$precision, "\n")
# Precision: 1 
# > cat("Recall:", metrics_hima$recall, "\n")
# Recall: 0.01666667 
# > cat("F1 Score:", metrics_hima$f1_score, "\n")
# F1 Score: 0.03278689 
# > 