# Existing setup code
rm(list = ls())

library(MASS)
library(stats)
library(graphics)
library(energy)
library(matrixStats)
library(devtools)
library(gtools)

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
#install.packages("microHIMA_1.0.tar.gz", repos = NULL, type = "source")
library(microHIMA) #install.packages("ncvreg") # install.packages("hommel")

set.seed(1234)

# Define the directory where the results are stored
output_dir <- "simulation_data_4"

# List all .rds files in the directory
rds_files <- list.files(output_dir, pattern = "*.rds", full.names = TRUE)
rds_files <- mixedsort(rds_files)

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
  
  list(Y_vector = Y_vector, T_vector = T_vector, M_nz_matrix = M_nz_matrix, M_matrix = M_matrix)
}

# Function to calculate metrics
calculate_metrics <- function(detected_indices, true_indices) {
  true_positives <- length(intersect(detected_indices, true_indices))
  false_positives <- length(setdiff(detected_indices, true_indices))
  false_negatives <- length(setdiff(true_indices, detected_indices))
  
  precision <- ifelse(true_positives + false_positives > 0, true_positives / (true_positives + false_positives), 0)
  recall <- ifelse(true_positives + false_negatives > 0, true_positives / (true_positives + false_negatives), 0)
  f1_score <- ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), 0)
  
  list(precision = precision, recall = recall, f1_score = f1_score,
       true_positives = true_positives, false_positives = false_positives, false_negatives = false_negatives)
}

# Initialize results data frame in long format
results_long_df <- data.frame(method = character(), iteration = integer(), true_positives = integer(),
                              false_positives = integer(), false_negatives = integer(), precision = numeric(), 
                              recall = numeric(), f1_score = numeric())

# Loop through each result and apply methods i=2
# for (i in 1:length(simulation_results)) {
  for (i in 1:2) {
  result <- simulation_results[[i]]
  spiked_features <- as.numeric(gsub("Feature", "", result$feature$feature_spiked))
  
  preprocessed_data <- preprocess_data(result)
  Y_vector <- preprocessed_data$Y_vector
  T_vector <- preprocessed_data$T_vector
  M_nz_matrix <- preprocessed_data$M_nz_matrix
  M_matrix <- preprocessed_data$M_matrix
  M_absolute <- result$absolute_M
  
  # Apply CCMM method
  result_ccmm <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
                      sig.level = 0.05, tol = 1e-06, max.iter = 5000)
  significant_indices_ccmm <- as.numeric(which(result_ccmm$IDE.CIs[1, ] > 0 | result_ccmm$IDE.CIs[2, ] < 0))
  metrics_ccmm <- calculate_metrics(significant_indices_ccmm, spiked_features)
  results_long_df <- rbind(results_long_df, data.frame(method = "ccmm", iteration = i, true_positives = metrics_ccmm$true_positives,
                                                       false_positives = metrics_ccmm$false_positives, false_negatives = metrics_ccmm$false_negatives,
                                                       precision = metrics_ccmm$precision, recall = metrics_ccmm$recall, f1_score = metrics_ccmm$f1_score))
  
  # Apply LDM method
  data <- data.frame(Y = Y_vector, T = T_vector)
  result_ldm <- ldm(
    formula = M_nz_matrix ~ T + Y,
    data = data,
    seed = 1234,
    test.mediation = TRUE
  )
  detected_otus_ldm <- as.numeric(gsub("taxon_", "", result_ldm$med.detected.otu.omni))
  metrics_ldm <- calculate_metrics(detected_otus_ldm, spiked_features)
  results_long_df <- rbind(results_long_df, data.frame(method = "ldm", iteration = i, true_positives = metrics_ldm$true_positives,
                                                       false_positives = metrics_ldm$false_positives, false_negatives = metrics_ldm$false_negatives,
                                                       precision = metrics_ldm$precision, recall = metrics_ldm$recall, f1_score = metrics_ldm$f1_score))
  
  # Apply HIMA method
  result_hima <- tryCatch({
    mhima(exposure = T_vector, covariates = NULL, otu.com = M_nz_matrix, outcome = Y_vector)
  }, error = function(e) {
    NULL
  })
  detected_otus_hima <- if (!is.null(result_hima)) result_hima$ID else integer(0)
  metrics_hima <- calculate_metrics(detected_otus_hima, spiked_features)
  results_long_df <- rbind(results_long_df, data.frame(method = "hima", iteration = i, true_positives = metrics_hima$true_positives,
                                                       false_positives = metrics_hima$false_positives, false_negatives = metrics_hima$false_negatives,
                                                       precision = metrics_hima$precision, recall = metrics_hima$recall, f1_score = metrics_hima$f1_score))
  
  # Apply microbvs method
  model_real <- MCMC_Med(trt = T_vector, Y = Y_vector, Z = M_nz_matrix, taxa = 2, iterations = 3000)
  result_global <- Selection_Med2(model = model_real)
  detected_otus_microbvs <- which(result_global$selected)
  metrics_microbvs <- calculate_metrics(detected_otus_microbvs, spiked_features)
  results_long_df <- rbind(results_long_df, data.frame(method = "microbvs", iteration = i, true_positives = metrics_microbvs$true_positives,
                                                       false_positives = metrics_microbvs$false_positives, false_negatives = metrics_microbvs$false_negatives,
                                                       precision = metrics_microbvs$precision, recall = metrics_microbvs$recall, f1_score = metrics_microbvs$f1_score))
  
  # Apply MedZIM method
  taxon_name <- "taxon_"
  
  libsize <- colSums(M_absolute)
  dat <- data.frame(Y_vector, T_vector, libsize)
  dat <- cbind(dat, M_matrix)
  
  MedZIM_results <- MedZIM_func(
    dat = dat,
    xVar = "T_vector",
    yVar = "Y_vector",
    taxon_name = taxon_name,
    libSize_name = "libsize",
    obs_gt_0 = 2,
    obs_eq_0 = 2,
    inter_x_mg0 = TRUE,
    inter_x_m = FALSE,
    eval.max = 200,
    iter.max = 200,
    x_from = 0,
    x_to = 1,
    type1error = 0.05,
    paraJobs = 2
  )
  
  detected_otus_medzim <- c()
  
  for (taxon in names(MedZIM_results$fullList)) {
    pNIE_value <- MedZIM_results$fullList[taxon][[1]]
    pNIE_value <- as.numeric(pNIE_value[[1]]["pNIE"])
    if (!is.na(pNIE_value) && pNIE_value < 0.05) {
      detected_otus_medzim <- c(detected_otus_medzim, taxon)
    }
  }
  
  detected_otus_medzim <- as.numeric(gsub("taxon_", "", detected_otus_medzim))
  
  metrics_medzim <- calculate_metrics(detected_otus_medzim, spiked_features)
  results_long_df <- rbind(results_long_df, data.frame(method = "medzim", iteration = i, true_positives = metrics_medzim$true_positives,
                                                       false_positives = metrics_medzim$false_positives, false_negatives = metrics_medzim$false_negatives,
                                                       precision = metrics_medzim$precision, recall = metrics_medzim$recall, f1_score = metrics_medzim$f1_score))
}

# Transform results_long_df to wide format
library(tidyr)

results_wide_df <- results_long_df %>%
  pivot_wider(names_from = method, values_from = c(true_positives, false_positives, false_negatives, precision, recall, f1_score))

# Save results_wide_df to CSV
write.csv(results_wide_df, file = "results_wide_df.csv", row.names = FALSE)

