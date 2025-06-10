rm(list = ls())

library(MASS)
library(stats)
library(graphics)
library(energy)  
library(matrixStats)
library(devtools)
library(foreach)
library(parallel)
library(doParallel)
library(gtools)
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
# install.packages("rrBLUP")
# install.packages("CompQuadForm")
source("highmed2019.r")
source('fromSKAT.R')
set.seed(1234)

#mediation effect

# read the results
output_dir <- "simulation_data_sample_new"

# List all .rds files in the directory
rds_files <- list.files(output_dir, pattern = "*.rds", full.names = TRUE)
rds_files <- mixedsort(rds_files)

# Function to preprocess a single result
preprocess_single_result <- function(result) {
  Y_vector <- as.numeric(result$Y)
  T_vector <- as.numeric(result$T)
  M_a_matrix <- as.matrix(t(result$absolute_M))
  colnames(M_a_matrix) <- paste0("taxon_", 1:ncol(M_a_matrix))
  
  # Extract M_matrix and M_nz_matrix using a separate function
  matrices <- compute_matrices(M_a_matrix)
  
  list(Y_vector = Y_vector, T_vector = T_vector, M_a_matrix = M_a_matrix, M_matrix = matrices$M_matrix, M_nz_matrix = matrices$M_nz_matrix, feature = result$feature)
}

# Function to compute M_matrix and M_nz_matrix
compute_matrices <- function(M_a_matrix) {
  M_matrix <- M_a_matrix / rowSums(M_a_matrix)
  
  temp <- M_a_matrix
  temp[temp == 0] <- 0.5
  M_nz_matrix <- temp / rowSums(temp)
  
  list(M_matrix = M_matrix, M_nz_matrix = M_nz_matrix)
}

# Function to preprocess all results in a category
preprocess_results <- function(results) {
  lapply(results, preprocess_single_result)
}

# Split the results into categories
split_results <- function(files) {
  results <- list()
  for (file in files) {
    result <- readRDS(file)
    # Assuming the file naming convention includes the variable name and value
    name_parts <- strsplit(basename(file), "_")[[1]]
    
    if (length(name_parts) < 4) {
      next
    }
    
    variable <- name_parts[2]
    value <- name_parts[3]
    
    if (!is.null(results[[variable]])) {
      if (!is.null(results[[variable]][[value]])) {
        results[[variable]][[value]] <- c(results[[variable]][[value]], list(result))
      } else {
        results[[variable]][[value]] <- list(result)
      }
    } else {
      results[[variable]] <- list()
      results[[variable]][[value]] <- list(result)
    }
  }
  return(results)
}

# Function to set zero percentage in M_matrix
set_zero_percentage <- function(M_matrix, percentage) {
  num_samples <- nrow(M_matrix)
  num_features <- ncol(M_matrix)
  
  for (i in 1:num_samples) {
    zero_count <- sum(M_matrix[i, ] == 0)
    total_elements <- num_features
    current_zero_percentage <- (zero_count / total_elements) * 100
    
    while (current_zero_percentage < percentage) {
      non_zero_values <- M_matrix[i, M_matrix[i, ] != 0]
      if (length(non_zero_values) == 0) break  # Avoid infinite loop
      min_value_index <- which(M_matrix[i, ] == min(non_zero_values))
      M_matrix[i, min_value_index] <- 0
      zero_count <- sum(M_matrix[i, ] == 0)
      current_zero_percentage <- (zero_count / total_elements) * 100
    }
  }
  
  return(M_matrix)
}


# 输入矩阵处理：函数接收一个计数矩阵，并确定需要子采样的深度和目标零值比例。
# 
# 子采样步骤：对每个样本进行子采样，以减少每个样本中的总读取数，同时保持样本内变量的相对比例。
# 
# 生成零值：在子采样后的矩阵中，根据计数值的大小生成一定比例的零值，计数值越小的元素越容易被置为零。
# 
# 输出矩阵：最终输出的是一个经过子采样和零值生成的稀疏矩阵。



# Split the results into categories
split_results <- split_results(rds_files)# Preprocess data for each category
preprocessed_results_by_sample <- lapply(split_results$n, preprocess_results)
preprocessed_data <- preprocessed_results_by_sample[["sample"]]

# Zero percentage settings
zero_percentages <- c(30)

for (percentage in zero_percentages) { 
  results <- list()
  
  for (i in 1:length(preprocessed_data)) { 
    loop_start_time <- Sys.time()  
    
    data_list <- preprocessed_data[i]
    data = data_list[[1]]
    
    # Extract necessary vectors and matrices
    Y_vector <- as.vector(data$Y_vector)
    T_vector <- as.vector(data$T_vector)
    M_matrix <- as.matrix(data$M_matrix)
    M_a_matrix <- as.matrix(data$M_a_matrix)
    
    # Set zero percentage for M_matrix
    M_a_matrix <- set_zero_percentage(M_a_matrix, percentage)
    M_matrix <- M_a_matrix / rowSums(M_a_matrix)
    
    # MODIMA with M_matrix
    modima_time_matrix <- system.time({
      T_dist <- dist(T_vector)
      Y_dist <- dist(Y_vector)
      M_dist_matrix <- dist(M_matrix)
      modima_result_matrix <- modima(T_dist, M_dist_matrix, Y_dist, nrep = 999)
      modima_p_value_matrix <- modima_result_matrix$p.value
    })
    print(paste("MODIMA (M_matrix) completed in:", modima_time_matrix[3], "seconds"))
    
    # MODIMA with M_a_matrix
    modima_time_amatrix <- system.time({
      M_dist_amatrix <- dist(M_a_matrix)
      modima_result_amatrix <- modima(T_dist, M_dist_amatrix, Y_dist, nrep = 999)
      modima_p_value_amatrix <- modima_result_amatrix$p.value
    })
    print(paste("MODIMA (M_a_matrix) completed in:", modima_time_amatrix[3], "seconds"))
    
    # MEDTEST with M_matrix
    medtest_time_matrix <- system.time({
      m_list_matrix <- list(euc = dist(M_matrix))
      medtest_result_matrix <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list_matrix, z = NULL, nperm = 999)
      medtest_p_value_matrix <- medtest_result_matrix$permP
    })
    print(paste("MEDTEST (M_matrix) completed in:", medtest_time_matrix[3], "seconds"))
    
    # MEDTEST with M_a_matrix
    medtest_time_amatrix <- system.time({
      m_list_amatrix <- list(euc = dist(M_a_matrix))
      medtest_result_amatrix <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list_amatrix, z = NULL, nperm = 999)
      medtest_p_value_amatrix <- medtest_result_amatrix$permP
    })
    print(paste("MEDTEST (M_a_matrix) completed in:", medtest_time_amatrix[3], "seconds"))
    
    loop_end_time <- Sys.time()  
    print(paste("Loop", i, "completed in:", loop_end_time - loop_start_time, "mins"))
    
    # Store results
    results[[i]] <- list(
      modima_p_value_matrix = modima_p_value_matrix,
      modima_p_value_amatrix = modima_p_value_amatrix,
      medtest_p_value_matrix = medtest_p_value_matrix,
      medtest_p_value_amatrix = medtest_p_value_amatrix
    )
  }
  
  # Convert results to a data frame
  results_df <- do.call(rbind, lapply(results, function(x) {
    c(modima_p_value_matrix = x$modima_p_value_matrix,
      modima_p_value_amatrix = x$modima_p_value_amatrix,
      medtest_p_value_matrix = x$medtest_p_value_matrix,
      medtest_p_value_amatrix = x$medtest_p_value_amatrix)
  }))
  
  # Convert the list of results into a data frame
  results_df <- as.data.frame(results_df)
  print(results_df)
  
  # Save the data frame to a CSV file
  write.csv(results_df, file = paste0("sample_", percentage, "_modima_medtest.csv"), row.names = FALSE)
  
  print(paste("Results for", percentage, "% zero values have been saved"))
}

print("All results have been saved")