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
output_dir <- "simulation_data_feature"

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
preprocessed_data <- preprocessed_results_by_sample[["feature"]]

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
zero_percentages <- c(10, 30, 50, 70)

for (percentage in zero_percentages) { #percentage =30
  results <- list()
  
  for (i in 1:length(preprocessed_data)) { #i=1
    loop_start_time <- Sys.time()  # Start timing for the entire loop i=1
    
    data_list <- preprocessed_data[i]
    data = data_list[[1]]
    
    # Extract necessary vectors and matrices
    Y_vector <- as.vector(data$Y_vector)
    T_vector <- as.vector(data$T_vector)
    M_matrix <- as.matrix(data$M_matrix)
    M_a_matrix <- as.matrix(data$M_a_matrix)
    
    # Set zero percentage
    M_matrix <- set_zero_percentage(M_matrix, percentage)
    #M_matrix= set_zero_percentage_new(count_matrix=M_matrix, target_depth = NULL, target_zero_proportion= percentage) 
    
    
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
  write.csv(results_df, file = paste0("feature_", percentage, ".csv"), row.names = FALSE)
  
  print(paste("Results for", percentage, "% zero values have been saved"))
}

print("All results have been saved")

############## end mediation effect


# read the results
output_dir <- "simulation_data_feature_no_mediation"

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

# Split the results into categories
split_results <- split_results(rds_files)

# Preprocess data for each category
preprocessed_results_by_sample <- lapply(split_results$n, preprocess_results)
preprocessed_data <- preprocessed_results_by_sample[["feature"]]

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
zero_percentages <- c(10, 30, 50, 70)

for (percentage in zero_percentages) {
  results <- list()
  
  for (i in 1:length(preprocessed_data)) {
    loop_start_time <- Sys.time()  # Start timing for the entire loop i=1 percentage = 30
    
    data_list <- preprocessed_data[i]
    data = data_list[[1]]
    
    # Extract necessary vectors and matrices
    Y_vector <- as.vector(data$Y_vector)
    T_vector <- as.vector(data$T_vector)
    M_matrix <- as.matrix(data$M_matrix)
    M_a_matrix <- as.matrix(data$M_a_matrix)
    
    # Set zero percentage
    M_matrix <- set_zero_percentage(M_matrix, percentage)
    
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
      ccmm_result <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
                          sig.level = 0.05, tol = 1e-06, max.iter = 5000)
      ccmm_CI_lower <- as.numeric(ccmm_result$TIDE.CI["2.5%"])
      ccmm_CI_upper <- as.numeric(ccmm_result$TIDE.CI["97.5%"])
    })
    print(paste("CCMM completed in:", ccmm_time[3], "seconds"))
    
    # LDM-MED
    ldm_time <- system.time({
      data_ldm <- data.frame(Y = Y_vector, T = T_vector, M_nz_matrix = M_matrix)
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
  write.csv(results_df, file = paste0("feature_noMe_", percentage, ".csv"), row.names = FALSE)
  
  print(paste("Results for", percentage, "% zero values have been saved"))
}

print("All results have been saved")










#####################################################################################3




# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Function to generate power results from a CSV file
generate_power_results <- function(file_name, zero_percentage) {
  results_df <- read.csv(file_name)
  
  # Add sample size information
  results_df$sample_size <- rep(c(10, 20, 40), each = 100)
  
  # Replace missing values with 1 (acceptance)
  results_df[is.na(results_df)] <- 1
  
  # Rename methods for better readability
  results_long <- melt(results_df, id.vars = "sample_size", 
                       measure.vars = c("modima_p_value", "medtest_p_value", "ldm_p_value", "permanova_p_value"),
                       variable.name = "Method", value.name = "p_value")
  
  results_long$Method <- recode(results_long$Method,
                                "modima_p_value" = "MODIMA",
                                "medtest_p_value" = "MedTest",
                                "ldm_p_value" = "LDM",
                                "permanova_p_value" = "PERMANOVA")
  
  # Calculate power (p-value < 0.05)
  power_results <- results_long %>%
    group_by(sample_size, Method) %>%
    summarise(power = mean(p_value < 0.05, na.rm = TRUE))
  
  # Calculate rejection rate based on confidence intervals (reject if not containing 0)
  results_df <- results_df %>%
    mutate(ccmm_reject = ifelse(ccmm_CI_lower > 0 | ccmm_CI_upper < 0, 1, 0))
  
  # Add ccmm_reject to power_results
  ccmm_power_results <- results_df %>%
    group_by(sample_size) %>%
    summarise(Method = "CCMM", power = mean(ccmm_reject, na.rm = TRUE))
  
  # Combine all methods' power_results
  power_results <- bind_rows(power_results, ccmm_power_results)
  
  # Replace zero power values with a small value
  power_results$power[power_results$power == 0] <- 0.001
  
  # Add zero percentage information
  power_results$Zero_Percentage <- zero_percentage
  
  return(power_results)
}

# Generate power results for each sample
power_results_10 <- generate_power_results("feature_10.csv", "10% Zero Values")
power_results_30 <- generate_power_results("feature_30.csv", "30% Zero Values")
power_results_50 <- generate_power_results("feature_50.csv", "50% Zero Values")
power_results_70 <- generate_power_results("feature_70.csv", "70% Zero Values")

# Combine all power results
all_power_results <- bind_rows(power_results_10, power_results_30, power_results_50, power_results_70)

# Plot power analysis with facet_grid
ggplot(all_power_results, aes(x = Method, y = power, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Power Analysis by Zero Value Percentage",
       x = "Method",
       y = "Power") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "feature Size") +
  coord_cartesian(ylim = c(0, 1)) +  # Fixed y-axis range
  facet_grid(. ~ Zero_Percentage)

######################################################################333333


############### type-1 error
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Function to generate Type-1 error results from a CSV file
generate_type1_error_results <- function(file_name, zero_percentage) {
  results_df <- read.csv(file_name)
  
  # Add sample size information
  results_df$sample_size <- rep(c(10, 20, 40), each = 100)
  
  # Replace missing values with 1 (acceptance)
  results_df[is.na(results_df)] <- 1
  
  # Rename methods for better readability
  results_long <- melt(results_df, id.vars = "sample_size", 
                       measure.vars = c("modima_p_value", "medtest_p_value"
                                        # "sparsemcmm_ome_p_value.OME", 
                                        , "ldm_p_value", "permanova_p_value"),
                       variable.name = "Method", value.name = "p_value")
  
  results_long$Method <- recode(results_long$Method,
                                "modima_p_value" = "MODIMA",
                                "medtest_p_value" = "MedTest",
                                # "sparsemcmm_ome_p_value.OME" = "SparseMCMM",
                                "ldm_p_value" = "LDM",
                                "permanova_p_value" = "PERMANOVA")
  
  # Calculate Type-1 error (p-value < 0.05)
  type1_error_results <- results_long %>%
    group_by(sample_size, Method) %>%
    summarise(type1_error = mean(p_value < 0.05, na.rm = TRUE))
  
  # Calculate rejection rate based on confidence intervals (reject if not containing 0)
  results_df <- results_df %>%
    mutate(ccmm_reject = ifelse(ccmm_CI_lower > 0 | ccmm_CI_upper < 0, 1, 0)
           # ,
           # microbvs_reject = ifelse(`microbvs_result.2.5.` > 0 | `microbvs_result.97.5.` < 0, 1, 0)
           
    )
  
  # Add ccmm_reject to type1_error_results
  ccmm_type1_error_results <- results_df %>%
    group_by(sample_size) %>%
    summarise(Method = "CCMM", type1_error = mean(ccmm_reject, na.rm = TRUE))
  
  # microbvs_type1_error_results <- results_no_mediation_df %>%
  #   group_by(sample_size) %>%
  #   summarise(Method = "MicroBVS", type1_error = mean(microbvs_reject, na.rm = TRUE))
  
  # Combine all methods' type1_error_results
  type1_error_results <- bind_rows(type1_error_results, ccmm_type1_error_results
                                   # , microbvs_type1_error_results
  )
  
  # Replace zero type1_error values with a small value
  type1_error_results$type1_error[type1_error_results$type1_error == 0] <- 0.001
  
  # Add zero percentage information
  type1_error_results$Zero_Percentage <- zero_percentage
  
  return(type1_error_results)
}

# Generate Type-1 error results for each sample
type1_error_results_10 <- generate_type1_error_results("feature_noMe_10.csv", "10% Zero Values")
type1_error_results_30 <- generate_type1_error_results("feature_noMe_30.csv", "30% Zero Values")
type1_error_results_50 <- generate_type1_error_results("feature_noMe_50.csv", "50% Zero Values")
type1_error_results_70 <- generate_type1_error_results("feature_noMe_70.csv", "70% Zero Values")

# Combine all Type-1 error results
all_type1_error_results <- bind_rows(type1_error_results_10, type1_error_results_30, type1_error_results_50, type1_error_results_70)

# Plot Type-1 error analysis with facet_grid
ggplot(all_type1_error_results, aes(x = Method, y = type1_error, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Type-1 Error Analysis by Zero Value Percentage",
       x = "Method",
       y = "Type-1 Error") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "feature Size") +
  coord_cartesian(ylim = c(0, 1)) +  # Fixed y-axis range
  facet_grid(. ~ Zero_Percentage)