rm(list = ls())

library(MASS)
library(stats)
library(graphics)
library(energy)  
library(matrixStats)
library(devtools)
##
source("modima.R")
source("MedOmniTest.R")
# install_github("chanw0/SparseMCMM")
library(SparseMCMM)
# install.packages("ccmm")
library(ccmm)
# install_github("yijuanhu/LDM")
library(LDM)
# install_github( "mkoslovsky/MicroBVS")
library(MicroBVS)
# devtools::install_github("quranwu/MedZIM")
library(MedZIM)

set.seed(1234)

# Define the directory where the results are stored
output_dir <- "simulation_data_2"

# List all .rds files in the directory
rds_files <- list.files(output_dir, pattern = "*.rds", full.names = TRUE)

# Source the functions
source("split_results.R")
source("preprocess_data.R")

# Split the results into categories
split_results <- split_results(rds_files)

# Preprocess data for each category
preprocessed_results_by_feature <- preprocess_results(split_results$results_by_feature)
preprocessed_results_by_sample <- preprocess_results(split_results$results_by_sample)
preprocessed_results_by_metadata_effect_size <- preprocess_results(split_results$results_by_metadata_effect_size)
preprocessed_results_by_perc_feature_spiked_metadata <- preprocess_results(split_results$results_by_perc_feature_spiked_metadata)
preprocessed_results_by_median_read_depth <- preprocess_results(split_results$results_by_median_read_depth)
preprocessed_results_by_noise_sd <- preprocess_results(split_results$results_by_noise_sd)

print("Data has been loaded and preprocessed successfully.")

apply_modima_medtest <- function(preprocessed_data) {
  results <- list()
  
  for (feature in names(preprocessed_data)) {
    data <- preprocessed_data[[feature]]
    
    # MODIMA
    T_dist <- dist(data$T_vector)
    Y_dist <- dist(data$Y_vector)
    M_dist <- dist(data$M_a_matrix)
    
    modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
    modima_p_value <- modima_result$p.value
    
    # MEDTEST
    m_list <- list(euc = dist(data$M_a_matrix))
    medtest_result <- MedOmniTest(x = data$T_vector, y = data$Y_vector, m.list = m_list, z = NULL, nperm = 999)
    medtest_p_value <- medtest_result$permP
    
    # Store results
    results[[feature]] <- list(
      modima_p_value = modima_p_value,
      medtest_p_value = medtest_p_value
    )
  }
  
  return(results)
}

# Apply MODIMA and MedOmniTest on preprocessed_results_by_feature
modima_medtest_results <- apply_modima_medtest(preprocessed_results_by_feature)
