rm(list = ls())
setwd("/home/haoyi/microbiome mediation analysis") 
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
set_zero_percentage_new <- function(count_matrix, target_depth = NULL) {
  # 如果没有指定 target_depth，使用每个样本中最小的总读取数
  if (is.null(target_depth)) {
    target_depth <- min(rowSums(count_matrix))
  }
  
  # 初始化一个矩阵用于存储子采样结果
  subsampled_matrix <- matrix(0, nrow = nrow(count_matrix), ncol = ncol(count_matrix))
  colnames(subsampled_matrix) <- colnames(count_matrix)
  rownames(subsampled_matrix) <- rownames(count_matrix)
  
  # 对每个样本（行）进行子采样
  for (i in 1:nrow(count_matrix)) {
    original_counts <- count_matrix[i, ]
    total_reads <- sum(original_counts)
    
    if (total_reads > target_depth) {
      # 按比例计算每个变量的抽样概率
      prob <- original_counts / total_reads
      
      # 根据概率分布随机抽取 target_depth 次读取
      selected_indices <- sample(
        rep(seq_along(original_counts), original_counts), 
        size = target_depth, 
        replace = FALSE
      )
      
      # 计算每个变量被抽中的次数
      subsampled_matrix[i, ] <- tabulate(selected_indices, nbins = length(original_counts))
    } else {
      # 如果总读取数小于或等于目标深度，直接保留原始计数
      subsampled_matrix[i, ] <- original_counts
    }
  }
  
  return(subsampled_matrix)
}
# 
# set_zero_percentage_new <- function(count_matrix, target_depth = NULL, target_zero_proportion) {
#   # Convert the target_zero_proportion from percentage to a proportion
#   target_zero_proportion <- target_zero_proportion / 100
#   
#   # If target_depth is not specified, use the minimum total read count across all samples
#   if (is.null(target_depth)) {
#     target_depth <- min(rowSums(count_matrix))
#   }
#   
#   # Initialize a matrix to store the subsampled results
#   subsampled_matrix <- matrix(0, nrow = nrow(count_matrix), ncol = ncol(count_matrix))
#   colnames(subsampled_matrix) <- colnames(count_matrix)
#   rownames(subsampled_matrix) <- rownames(count_matrix)
#   
#   # Perform subsampling for each sample (row)
#   for (i in 1:nrow(count_matrix)) {
#     original_counts <- count_matrix[i, ]
#     total_reads <- sum(original_counts)
#     
#     if (total_reads > target_depth) {
#       # Calculate the sampling probability for each variable based on its proportion of the total reads
#       prob <- original_counts / total_reads
#       
#       # Randomly select target_depth reads according to the sampling probability distribution
#       selected_indices <- sample(
#         rep(seq_along(original_counts), original_counts), 
#         size = target_depth, 
#         replace = FALSE
#       )
#       
#       # Calculate the number of times each variable was selected
#       subsampled_matrix[i, ] <- tabulate(selected_indices, nbins = length(original_counts))
#     } else {
#       # If the total read count is less than or equal to the target depth, keep the original counts
#       subsampled_matrix[i, ] <- original_counts
#     }
#     
#     # Generate zero values
#     non_zero_indices <- which(subsampled_matrix[i, ] != 0)
#     current_zero_count <- sum(subsampled_matrix[i, ] == 0)
#     target_zero_count <- ceiling(target_zero_proportion * length(subsampled_matrix[i, ]))
#     
#     additional_zeros_needed <- target_zero_count - current_zero_count
#     
#     if (additional_zeros_needed > 0) {
#       # Calculate the probability of generating zeros, inversely proportional to the count value
#       probabilities <- 1 / (subsampled_matrix[i, non_zero_indices] + 1)
#       probabilities <- probabilities / sum(probabilities)  # Normalize probabilities
#       
#       # Select elements to be set to zero based on the calculated probabilities
#       zero_indices <- sample(non_zero_indices, size = min(additional_zeros_needed, length(non_zero_indices)), 
#                              prob = probabilities, replace = FALSE)
#       
#       subsampled_matrix[i, zero_indices] <- 0
#     }
#   }
#   
#   return(subsampled_matrix)
# }
# 输入矩阵处理：函数接收一个计数矩阵，并确定需要子采样的深度和目标零值比例。
# 
# 子采样步骤：对每个样本进行子采样，以减少每个样本中的总读取数，同时保持样本内变量的相对比例。
# 
# 生成零值：在子采样后的矩阵中，根据计数值的大小生成一定比例的零值，计数值越小的元素越容易被置为零。
# 
# 输出矩阵：最终输出的是一个经过子采样和零值生成的稀疏矩阵。



# Split the results into categories
split_results <- split_results(rds_files)

# Preprocess data for each category
preprocessed_results_by_sample <- lapply(split_results$n, preprocess_results)
preprocessed_data <- preprocessed_results_by_sample[["sample"]]

#i=2 i=3

# for (i in 1:length(preprocessed_data)){
#   data_list <- preprocessed_data[i]
#   data = data_list[[1]]
#   # data$feature
#   # Extract necessary vectors and matrices
#   Y_vector <- as.vector(data$Y_vector)
#   T_vector <- as.vector(data$T_vector)
#   M_matrix <- as.matrix(data$M_matrix)
#   M_a_matrix <- as.matrix(data$M_a_matrix)
#   M_nz_matrix <- as.matrix(data$M_nz_matrix)
#   zero_count <- sum(M_matrix == 0)
#   total_elements <- length(M_matrix)
#   zero_proportion <- (zero_count / total_elements) * 100
#   print(i)
#   print(zero_proportion)
#   print(data$feature)
# 
# }

# Zero percentage settings
zero_percentages <- c(1)

for (percentage in zero_percentages) { #percentage =30
  results <- list()
  
  for (i in 1:length(preprocessed_data)) { #i=3
    loop_start_time <- Sys.time()  # Start timing for the entire loop i=1
    
    data_list <- preprocessed_data[i]
    data = data_list[[1]]
    
    # Extract necessary vectors and matrices
    Y_vector <- as.vector(data$Y_vector)
    T_vector <- as.vector(data$T_vector)
    M_matrix <- as.matrix(data$M_matrix)
    M_a_matrix <- as.matrix(data$M_a_matrix)
    
    # Set zero percentage
    # sum(M_matrix==0)
    #M_matrix <- set_zero_percentage(M_matrix, percentage)
    #M_matrix= set_zero_percentage_new(count_matrix=M_matrix, target_depth = NULL, target_zero_proportion= percentage) 
    M_a_matrix= set_zero_percentage_new(count_matrix=M_a_matrix, target_depth = 50) 
    M_matrix <- M_a_matrix / rowSums(M_a_matrix)
    
    # Compute M_nz_matrix after setting zero percentage
    temp <- M_matrix
    temp[temp == 0] <- 0.5
    M_nz_matrix <- temp / rowSums(temp)
    
    # MODIMA
    modima_time <- system.time({
      T_dist <- dist(T_vector)
      Y_dist <- dist(Y_vector)
      M_dist <- dist(M_matrix)
      modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
      modima_p_value <- modima_result$p.value
    })
    print(paste("MODIMA completed in:", modima_time[3], "seconds"))
    
    # MEDTEST
    medtest_time <- system.time({
      m_list <- list(euc = dist(M_matrix))
      medtest_result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)
      medtest_p_value <- medtest_result$permP
    })
    print(paste("MEDTEST completed in:", medtest_time[3], "seconds"))
    
    # CCMM
    ccmm_time <- system.time({
      
      tryCatch({
        ccmm_result <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
                            sig.level = 0.05, tol = 1e-06, max.iter = 5000)
        ccmm_CI_lower <- as.numeric(ccmm_result$TIDE.CI["2.5%"])
        ccmm_CI_upper <- as.numeric(ccmm_result$TIDE.CI["97.5%"])
      }, error = function(e) {
        
        ccmm_CI_lower <- -1
        ccmm_CI_upper <- 1
        
        print(paste("Error in CCMM:", e$message))
      })
    })
    print(paste("CCMM completed in:", ccmm_time[3], "seconds")) 
    
    # LDM-MED
    ldm_time <- system.time({
      data_ldm <- data.frame(Y = Y_vector, T = T_vector, M_nz_matrix = M_nz_matrix)
      ldm_p_value <- tryCatch({
        ldm_result <- ldm(
          formula = M_nz_matrix ~ T + Y,
          data = data_ldm,
          seed = 1234,
          test.mediation = TRUE
        )
        ldm_result$med.p.global.omni
      }, error = function(e) {
        NULL  
      })
    })
    print(paste("LDM-MED completed in:", ldm_time[3], "seconds"))
    
    # PERMANOVA-MED
    permanova_time <- system.time({
      permanova_p_value <- tryCatch({
        data_ldm <- data.frame(Y = Y_vector, T = T_vector)
        formula <- as.formula("M_nz_matrix ~ T + Y")
        res_perm <- permanovaFL(
          formula = formula,
          data = data_ldm,
          dist.method = "bray",  
          seed = 1234,
          n.perm.max = 1000,
          n.cores = 4,
          test.mediation = TRUE,
          verbose = TRUE
        )
        res_perm$med.p.permanova
      }, error = function(e) {
        NULL  
      })
    })
    print(paste("PERMANOVA-MED completed in:", permanova_time[3], "seconds"))
    
    loop_end_time <- Sys.time()  # End timing for the entire loop
    print(paste("Loop", i, "completed in:", loop_end_time - loop_start_time, "mins"))
    
    # Store results
    results[[i]] <- list(
      modima_p_value = modima_p_value,
      medtest_p_value = medtest_p_value,
      ccmm_CI_lower = ccmm_CI_lower,
      ccmm_CI_upper = ccmm_CI_upper,
      ldm_p_value = ldm_p_value,
      permanova_p_value = permanova_p_value
    )
  }
  
  # Convert results to a data frame
  results_df <- do.call(rbind, lapply(results, function(x) {
    c(modima_p_value = x$modima_p_value,
      medtest_p_value = x$medtest_p_value,
      ccmm_CI_lower = x$ccmm_CI_lower,
      ccmm_CI_upper = x$ccmm_CI_upper,
      ldm_p_value = x$ldm_p_value,
      permanova_p_value = x$permanova_p_value)
  }))
  
  # Convert the list of results into a data frame
  results_df <- as.data.frame(results_df)
  print(results_df)
  # Save the data frame to a CSV file
  write.csv(results_df, file = paste0("sample_subsample_", percentage, ".csv"), row.names = FALSE)
  
  print(paste("Results for", percentage, "% zero values have been saved"))
}

print("All results have been saved")