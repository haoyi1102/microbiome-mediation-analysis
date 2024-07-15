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

synthetic_data<-preprocessed_results_by_sample[["300"]]
synthetic_data<-synthetic_data$M_matrix
zero_count <- sum(synthetic_data == 0)
total_elements <- length(synthetic_data)
zero_proportion <- zero_count / total_elements
zero_proportion

print("Data has been loaded and preprocessed successfully.")

###############
# data <- preprocessed_results_by_feature[["10"]]
# Y_vector <- as.vector(data$Y_vector)  
# T_vector <- as.vector(data$T_vector)  
# M_a_matrix <- as.matrix(data$M_a_matrix)
# M_nz_matrix <- as.matrix(data$M_nz_matrix)
# ldm_data <- data.frame(Y = Y_vector, T = T_vector)
# 
# ldm_result <- ldm(
#   formula = M_nz_matrix ~ T + Y,
#   data = ldm_data,
#   seed = 1234,
#   test.mediation = TRUE
# )

###############


mediation_method <- function(preprocessed_data) {
  results <- list()
  
  for (feature in names(preprocessed_data)) {
   
    data <- preprocessed_data[[feature]]
    
    # data <- preprocessed_results_by_feature[["10"]]
    
    # Extract necessary vectors and matrices
    Y_vector <- as.vector(data$Y_vector)  
    T_vector <- as.vector(data$T_vector)  
    M_a_matrix <- as.matrix(data$M_a_matrix)
    M_nz_matrix <- as.matrix(data$M_nz_matrix)
    #print(M_nz_matrix)
    
    # MODIMA
    T_dist <- dist(T_vector)
    Y_dist <- dist(Y_vector)
    M_dist <- dist(M_a_matrix)
    
    modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
    modima_p_value <- modima_result$p.value
    
    # MEDTEST
    m_list <- list(euc = dist(M_a_matrix))
    medtest_result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)
    medtest_p_value <- medtest_result$permP
    
    # CCMM
    ccmm_result <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
                        sig.level = 0.05, tol = 1e-06, max.iter = 5000)
    ccmm_CI_lower <- as.numeric(ccmm_result$TIDE.CI["2.5%"])
    ccmm_CI_upper <- as.numeric(ccmm_result$TIDE.CI["97.5%"])
    
    # SparseMCMM
    # sparsemcmm_result <- SparseMCMM(T_vector, M_nz_matrix, Y_vector, n.split=1, num.per=200)
    # sparsemcmm_ome_p_value <- sparsemcmm_result$Test["OME"]
    sparsemcmm_ome_p_value <- 0
    
    # LDM-MED
    # ldm_data <- data.frame(Y = Y_vector, T = T_vector)
    # 
    # ldm_result <- ldm(
    #   formula = M_nz_matrix ~ T + Y,
    #   data = ldm_data,
    #   seed = 1234,
    #   test.mediation = TRUE
    # )
    
    #ldm_p_value <- ldm_result$med.p.global.omni
    ldm_p_value <- 0 
    
    # Store results in a data frame
    method_results <- data.frame(
      Method = c("MODIMA", "MedOmniTest", "CCMM", "SparseMCMM", "LDM-MED"),
      p_value = c(modima_p_value, medtest_p_value, NA, sparsemcmm_ome_p_value, ldm_p_value),
      CI_lower = c(NA, NA, ccmm_CI_lower, NA, NA),
      CI_upper = c(NA, NA, ccmm_CI_upper, NA, NA)
    )
    
    results[[feature]] <- method_results
  }
  
  return(results)
}

# Apply the updated function on preprocessed_results_by_feature
results_of_feature <- mediation_method(preprocessed_results_by_feature)
results_of_sample <- mediation_method(preprocessed_results_by_sample)
results_of_perc_feature_spiked_metadata <- mediation_method(preprocessed_results_by_perc_feature_spiked_metadata)
results_of_metadata_effect_size <- mediation_method(preprocessed_results_by_metadata_effect_size)
results_of_median_read_depth <- mediation_method(preprocessed_results_by_median_read_depth)
results_of_noise_sd  <- mediation_method(preprocessed_results_by_noise_sd )

# results_of_sample
# $`100`
# Method p_value CI_lower  CI_upper
# 1      MODIMA   0.160       NA        NA
# 2 MedOmniTest   0.121       NA        NA
# 3        CCMM      NA -1.46156 0.9859121
# 4  SparseMCMM   0.000       NA        NA
# 5     LDM-MED   0.000       NA        NA
# 
# $`300`
# Method p_value  CI_lower CI_upper
# 1      MODIMA   0.853        NA       NA
# 2 MedOmniTest   0.869        NA       NA
# 3        CCMM      NA -0.793345 2.756859
# 4  SparseMCMM   0.000        NA       NA
# 5     LDM-MED   0.000        NA       NA
# 
# $`50`
# Method p_value  CI_lower  CI_upper
# 1      MODIMA   0.426        NA        NA
# 2 MedOmniTest   0.736        NA        NA
# 3        CCMM      NA -2.190693 0.9168572
# 4  SparseMCMM   0.000        NA        NA
# 5     LDM-MED   0.000        NA        NA
# 
# > 
