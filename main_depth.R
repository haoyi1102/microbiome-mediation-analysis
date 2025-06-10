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
output_dir <- "10-10-new_method_simulation_data_depth_new"

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



# Split the results into categories
split_results <- split_results(rds_files)

# Preprocess data for each category
preprocessed_results_by_sample <- lapply(split_results$median, preprocess_results)
preprocessed_data <- preprocessed_results_by_sample[["read"]]



# Zero percentage settings
zero_percentages <- c(1,1000,100)
# 

###################################33￥￥￥￥￥￥￥￥￥￥￥￥￥￥444

# library(ggplot2)
# 
# # 初始化一个空向量用于存储每个数据集中零值的百分比
# zero_percentages <- numeric(length(preprocessed_data))
# 
# # 遍历 preprocessed_data，计算每个 M_matrix 中的零值百分比
# for (i in 1:length(preprocessed_data)) {
#   data <- preprocessed_data[[i]]  # 获取第i个数据集
#   M_matrix <- as.matrix(data$M_matrix)  # 转换为矩阵
# 
#   # 计算零值百分比
#   zero_percentages[i] <- sum(M_matrix == 0) / length(M_matrix) * 100  # 转换为百分比
# }
# 
# # 计算总的平均零值百分比
# average_zero_percentage <- mean(zero_percentages)
# 
# # 打印平均零值百分比
# cat("总共平均的零值百分比为: ", average_zero_percentage, "%\n")
# 
# # 将零值百分比转化为数据框，用于绘制频率图
# zero_percent_df <- data.frame(
#   Zero_Percentage = zero_percentages  # 每个数据集的零值百分比
# )

# # 绘制零值百分比的频率图（直方图）
# ggplot(zero_percent_df, aes(x = Zero_Percentage)) +
#   geom_histogram(binwidth = 5, fill = "steelblue", color = "black", alpha = 0.7) +
#   theme_minimal() +
#   labs(title = "零值百分比的频率分布",
#        x = "零值百分比 (%)",
#        y = "频率") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


###################################￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥444444

for (percentage in zero_percentages) { #percentage = 1
  results <- list()
  
  for (i in 1:length(preprocessed_data)) { #i=52
    loop_start_time <- Sys.time()  # Start timing for the entire loop i=90
    
    data_list <- preprocessed_data[i]
    data = data_list[[1]]
    #percentage =1
    # Extract necessary vectors and matrices
    Y_vector <- as.vector(data$Y_vector)
    T_vector <- as.vector(data$T_vector)
    M_matrix <- as.matrix(data$M_matrix)
    M_a_matrix <- as.matrix(data$M_a_matrix)
    
    # Set zero percentage
    # 
    
    
    
    if (percentage != 1) {
      M_a_matrix <- set_zero_percentage_new(count_matrix = M_a_matrix, target_depth = percentage)
    }
    M_matrix <- M_a_matrix / rowSums(M_a_matrix)
    print(sum(M_matrix==0)/length(M_matrix))
    # M_matrix[is.na(M_matrix)] <- 0
    # sum(M_matrix==0) sum(M_matrix!=0)
    
    # print("%:")
    
    # Compute M_nz_matrix after setting zero percentage
    temp <- M_matrix
    temp[temp == 0] <- 0.5
    M_nz_matrix <- temp / rowSums(temp)
    
    # MODIMA
    modima_time <- system.time({
      modima_p_value <- tryCatch({
        T_dist <- dist(T_vector)
        Y_dist <- dist(Y_vector)
        M_dist <- dist(M_matrix)
        modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
        modima_result$p.value
      }, error = function(e) {
        print(paste("Error in MODIMA:", e$message))
        
      })
    })
    print(paste("MODIMA completed in:", modima_time[3], "seconds"))
    
    
    
    # mediationStat = function(dE, dM, dR){
    #   EM <- energy::bcdcor(dE, dM)
    #   MRE <- energy::pdcor(dM, dR, dE)
    #   res <- EM*MRE
    #   return(res)
    # }
    # #MODIMA p-vals
    # energy_p = function(E, M, R, nrep=999){
    #   dE <- stats::dist(E)
    #   dM <- stats::dist(M)
    #   dR <- stats::dist(R)
    #   p1 = energy::dcor.test(dE, dM, R = nrep)$p.value
    #   p2 = energy::dcor.test(dM, dR, R = nrep)$p.value
    #   p3 = spdcov.test(dM, dR, dE, R = nrep)$p.value
    #   p4 = spdcov.test(dR, dM, dE, R = nrep)$p.value
    #   p5 = energy::pdcov.test(dM, dR, dE, R = nrep)$p.value
    #   c(mediationStat(dE, dM, dR), p1, p2, p3, p4, p5, max(p1, p2, p3, p4, p5))
    # }
    # energy.p = energy_p(E=T_vector, M=M_matrix, R=Y_vector)
    # 
    # # MEDTEST
    
    # MEDTEST
    medtest_time <- system.time({
      medtest_p_value <- tryCatch({
        m_list <- list(euc = dist(M_matrix))
        medtest_result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)
        medtest_result$permP
      }, error = function(e) {
        print(paste("Error in MEDTEST:", e$message))
        1  # 返回1表示错误发生
      })
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
  
  write.csv(results_df, file = paste0("10_10_depth_", percentage, ".csv"), row.names = FALSE)
  
  print(paste("Results for", percentage, "% zero values have been saved"))
}

print("All results have been saved")







##################

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Function to generate power results from a CSV file
generate_power_results <- function(file_name, zero_percentage) {
  results_df <- read.csv(file_name)
  
  # Add sample size information
  results_df$sample_size <- rep(c(1000, 10000, 100000), each = 50)
  
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
  power_results$Zero_Percentage <- factor(zero_percentage, levels = c("low", "medium", "high"))
  
  
  return(power_results)
}

# Generate power results for each sample
power_results_30 <- generate_power_results("10_10_depth_1.csv", "low")
power_results_50 <- generate_power_results("10_10_depth_1000.csv", "medium")
power_results_70 <- generate_power_results("10_10_depth_100.csv", "high")

# Combine all power results
all_power_results <- bind_rows( power_results_30, power_results_50, power_results_70)

# Plot power analysis with facet_grid
ggplot(all_power_results, aes(x = Method, y = power, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Power Analysis by Zero Value Percentage",
       x = "Method",
       y = "Power") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Depth") +
  coord_cartesian(ylim = c(0, 1)) +  # Fixed y-axis range
  facet_grid(. ~ Zero_Percentage)


ggsave(filename = "10-10-pic/sample_Prevalence_1_hist.png", plot = last_plot(), width = 10, height = 8)

