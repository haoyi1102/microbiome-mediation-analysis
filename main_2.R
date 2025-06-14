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

# if(!require(mediation)) install.packages("mediation")
library(mediation)

# 加载并行计算包
library(parallel)
library(doParallel)

set.seed(1234)

# Define the directory where the results are stored
output_dir <- "simulation_data_3"

# List all .rds files in the directory
rds_files <- list.files(output_dir, pattern = "*.rds", full.names = TRUE)

# Load all .rds files into a list
simulation_results <- lapply(rds_files, readRDS)

# Initialize counters for scores
modima_score <- 0
medtest_score <- 0
simple_mediation_score <- 0
total_tests <- length(simulation_results)

# 设置并行计算
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Loop through each result and apply MODIMA and MedOmniTest

for (i in 1:total_tests) {
  result <- simulation_results[[i]]
  
  # Extract necessary data
  Y_vector <- as.numeric(result$Y)
  T_vector <- as.numeric(result$T)
  M_a_matrix <- as.matrix(t(result$absolute_M))
  colnames(M_a_matrix) <- paste0("taxon_", 1:ncol(M_a_matrix))
  
  # MODIMA
  T_dist <- dist(T_vector)
  Y_dist <- dist(Y_vector)
  M_dist <- dist(M_a_matrix)
  
  modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
  modima_p_value <- modima_result$p.value
  
  if (modima_p_value < 0.05) {
    modima_score <- modima_score + 1
  }
  
  # MEDTEST
  m_list <- list(euc = dist(M_a_matrix))
  medtest_result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)
  medtest_p_value <- medtest_result$permP
  
  if (medtest_p_value < 0.05) {
    medtest_score <- medtest_score + 1
  }
  
  
  # # Simple mediation analysis for each mediator using parallel computing
  # mediation_p_values <- foreach(j = 1:ncol(M_a_matrix), .combine = c, .packages = "mediation") %dopar% {
  #   M_vector <- M_a_matrix[, j]
  #   
  #   model.M <- lm(M_vector ~ T_vector)
  #   model.Y <- lm(Y_vector ~ T_vector + M_vector)
  #   
  #   med.out <- mediate(model.M, model.Y, treat = "T_vector", mediator = "M_vector", boot = TRUE, sims = 100)
  #   return(med.out$d.avg.p)
  # }
  # 
  # 
  # # Combine the p-values from multiple mediation analyses
  # combined_p_value <- min(mediation_p_values) 
  # 
  # if (combined_p_value < 0.5) {
  #   simple_mediation_score <- simple_mediation_score + 1
  # }
  # 
  # print(mediation_p_values)
  print(i)
}

stopCluster(cl)

# Calculate accuracies
modima_accuracy <- modima_score / total_tests
medtest_accuracy <- medtest_score / total_tests
# simple_mediation_accuracy <- simple_mediation_score / total_tests

# Print results
print(paste("MODIMA accuracy:", modima_accuracy))
print(paste("MedOmniTest accuracy:", medtest_accuracy))
# print(paste("Simple mediation accuracy:", simple_mediation_accuracy))

# #> # Print results
# > print(paste("MODIMA accuracy:", modima_accuracy))
# [1] "MODIMA accuracy: 0.56"
# > print(paste("MedOmniTest accuracy:", medtest_accuracy))
# [1] "MedOmniTest accuracy: 0.38"
