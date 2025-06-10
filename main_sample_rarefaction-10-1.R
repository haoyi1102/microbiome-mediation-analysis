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

##正确的

#mediation effect

# read the results
output_dir <- "10-10-new_method_simulation_data_sample_new"
#output_dir <- "test"
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
preprocessed_results_by_sample <- lapply(split_results$n, preprocess_results)
preprocessed_data <- preprocessed_results_by_sample[["sample"]]



# Zero percentage settings
#zero_percentages <- c(1,1000,100)
zero_percentages <- c(1)
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
  
  for (i in 1:length(preprocessed_data)) { # i=1
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
    
    # # MODIMA
    # modima_time <- system.time({
    #   modima_p_value <- tryCatch({
    #     T_dist <- dist(T_vector)
    #     Y_dist <- dist(Y_vector)
    #     M_dist <- dist(M_matrix)
    #     modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
    #     modima_result$p.value
    #   }, error = function(e) {
    #     print(paste("Error in MODIMA:", e$message))
    #     
    #   })
    # })
    # print(paste("MODIMA completed in:", modima_time[3], "seconds"))
    
    
    
   

    # MEDTEST
    # medtest_time <- system.time({
    #   medtest_p_value <- tryCatch({
    #     m_list <- list(euc = dist(M_matrix))
    #     medtest_result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)
    #     medtest_result$permP
    #   }, error = function(e) {
    #     print(paste("Error in MEDTEST:", e$message))
    #     1  # 返回1表示错误发生
    #   })
    # })
    # print(paste("MEDTEST completed in:", medtest_time[3], "seconds"))
    
    #CCMM
    # ccmm_time <- system.time({
    # 
    #   tryCatch({
    #     ccmm_result <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
    #                         sig.level = 0.05, tol = 1e-06, max.iter = 5000)
    #     ccmm_CI_lower <- as.numeric(ccmm_result$TIDE.CI["2.5%"])
    #     ccmm_CI_upper <- as.numeric(ccmm_result$TIDE.CI["97.5%"])
    #   }, error = function(e) {
    # 
    #     ccmm_CI_lower <- -1
    #     ccmm_CI_upper <- 1
    # 
    #     print(paste("Error in CCMM:", e$message))
    #   })
    # })
    # print(paste("CCMM completed in:", ccmm_time[3], "seconds"))
    
    # LDM-MED
    # ldm_time <- system.time({
    #   data_ldm <- data.frame(Y = Y_vector, T = T_vector, M_nz_matrix = M_nz_matrix)
    #   ldm_p_value <- tryCatch({
    #     ldm_result <- ldm(
    #       formula = M_nz_matrix ~ T + Y,
    #       data = data_ldm,
    #       seed = 1234,
    #       test.mediation = TRUE
    #     )
    #     ldm_result$med.p.global.omni
    #   }, error = function(e) {
    #     NULL  
    #   })
    # })
    # print(paste("LDM-MED completed in:", ldm_time[3], "seconds"))
    
    # PERMANOVA-MED
    # permanova_time <- system.time({
    #   permanova_p_value <- tryCatch({
    #     data_ldm <- data.frame(Y = Y_vector, T = T_vector)
    #     formula <- as.formula("M_nz_matrix ~ T + Y")
    #     res_perm <- permanovaFL(
    #       formula = formula,
    #       data = data_ldm,
    #       dist.method = "bray",
    #       seed = 1234,
    #       n.perm.max = 1000,
    #       n.cores = 4,
    #       test.mediation = TRUE,
    #       verbose = TRUE
    #     )
    #     res_perm$med.p.permanova
    #   }, error = function(e) {
    #     NULL
    #   })
    # })
    # print(paste("PERMANOVA-MED completed in:", permanova_time[3], "seconds"))
    
    
    # MicroBVS
    microbvs_time <- system.time({
      tryCatch({
        # Run the MCMC model for MicroBVS
        model_real <- MCMC_Med(trt = T_vector, Y = Y_vector, Z = M_matrix, taxa = 2)
        
        # Perform mediation analysis and extract global quantile results
        result_global <- Selection_Med2(model = model_real)
        microbvs_global_lower <- result_global$global_quantile[1]  # Lower bound of the credible interval
        microbvs_global_upper <- result_global$global_quantile[2]  # Upper bound of the credible interval
        
      }, error = function(e) {
        # Handle errors gracefully
        microbvs_global_lower <- NA
        microbvs_global_upper <- NA
        
        print(paste("Error in MicroBVS:", e$message))
      })
    })
    print(paste("MicroBVS completed in:", microbvs_time[3], "seconds"))
    
    # SparseMCMM
    sparsemcmm_time <- system.time({
      tryCatch({
        # Run SparseMCMM
        sparsemcmm_result <-SparseMCMM(T_vector, M_nz_matrix, Y_vector, n.split = 1, num.per = 10)
        
        # Extract the p-value for OME
        sparsemcmm_ome_p_value <- sparsemcmm_result$Test["OME"]
        
      }, error = function(e) {
        # Handle errors gracefully
        sparsemcmm_ome_p_value <- NA
        print(paste("Error in SparseMCMM:", e$message))
      })
    })
    
    print(paste("SparseMCMM completed in:", sparsemcmm_time[3], "seconds"))
    
     
    #result_hima <- mhima(exposure = T_vector, covariates = NULL, otu.com = M_nz_matrix, outcome = Y_vector)

    # # taxon_name <- "taxon_"
    # # 
    # # libsize <- colSums(M_a_matrix)
    # # dat <- data.frame(Y_vector, T_vector, libsize)
    # # dat <- cbind(dat, M_matrix)
    # # 
    # # MedZIM_results <- MedZIM_func(
    # #   dat = dat,
    # #   xVar = "T_vector",
    # #   yVar = "Y_vector",
    # #   taxon_name = taxon_name,
    # #   libSize_name = "libsize",
    # #   obs_gt_0 = 2,
    # #   obs_eq_0 = 2,
    # #   inter_x_mg0 = TRUE,
    # #   inter_x_m = FALSE,
    # #   eval.max = 200,
    # #   iter.max = 200,
    # #   x_from = 0,
    # #   x_to = 1,
    # #   type1error = 0.05,
    # #   paraJobs = 2
    # # )
    # 
    # X_npem <- as.data.frame(T_vector)
    # colnames(X_npem) <- "Treatment"
    # M_npem <- apply(M_matrix, c(1, 2), function(s) log(s + 1))
    # Y_npem <- as.factor(Y_vector)
    # 
    # # Run NPEM method with Univariate Mediator
    # result_npem <- NPEM(X_npem, M_npem, Y_npem, method = "UV")
    # 
    # mcmc_end_time <- Sys.time()
    # mcmc_time_taken <- as.numeric(difftime(mcmc_end_time, mcmc_start_time, units = "secs"))
    # print(paste("MCMC_Med completed in:", mcmc_time_taken, "seconds"))
    # ###
    
    
    loop_end_time <- Sys.time()  # End timing for the entire loop
    print(paste("Loop", i, "completed in:", loop_end_time - loop_start_time, "mins"))
    
    # Store results
    # results[[i]] <- list(
    #   modima_p_value = modima_p_value,
    #   medtest_p_value = medtest_p_value,
    #   ccmm_CI_lower = ccmm_CI_lower,
    #   ccmm_CI_upper = ccmm_CI_upper,
    #   ldm_p_value = ldm_p_value,
    #   permanova_p_value = permanova_p_value
    # )
     results[[i]] <- list(
       
       microbvs_CI_lower =microbvs_global_lower,
       microbvs_CI_upper = microbvs_global_upper,
       sparsemcmm_p_value = sparsemcmm_ome_p_value
     )
    
  }
  
  # Convert results to a data frame
  results_df <- do.call(rbind, lapply(results, function(x) {
    c(
      microbvs_CI_lower = x$microbvs_CI_lower,
     
      microbvs_CI_upper = x$microbvs_CI_upper,
      sparsemcmm_p_value = x$sparsemcmm_p_value)
  }))
  
  # Convert the list of results into a data frame
  results_df <- as.data.frame(results_df)
  print(results_df)
  # Save the data frame to a CSV file
 
  write.csv(results_df, file = paste0("10-10-different-gen-1.1_sample_test_A", percentage, ".csv"), row.names = FALSE)
  
  print(paste("Results for", percentage, "% zero values have been saved"))
 }

print("All results have been saved")




#################


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
output_dir <- "simulation_data_sample_new_no_mediation"

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

# set_zero_percentage_new <- function(count_matrix, scaling_factor = 1) {
#   # 计算概率矩阵，概率与计数值成反比
#   prob_matrix <- scaling_factor / (count_matrix + scaling_factor)
#   
#   # 初始化一个矩阵用于存储结果
#   zero_matrix <- matrix(0, nrow = nrow(count_matrix), ncol = ncol(count_matrix))
#   colnames(zero_matrix) <- colnames(count_matrix)
#   rownames(zero_matrix) <- rownames(count_matrix)
#   
#   # 使用伯努利试验生成零值
#   for (i in 1:nrow(count_matrix)) {
#     for (j in 1:ncol(count_matrix)) {
#       zero_matrix[i, j] <- ifelse(runif(1) < prob_matrix[i, j], 0, count_matrix[i, j])
#     }
#   }
#   return(zero_matrix)
# }
# 
# set_zero_percentage_new <- function(count_matrix, target_zero_proportion) {
#   # 初始的 scaling_factor 值
#   scaling_factor <- 1
#   current_zero_proportion <- 0
#   
#   # 第一步：使用伯努利随机变量初步生成零值
#   while (abs(current_zero_proportion - target_zero_proportion) > 1) {
#     # 计算概率矩阵
#     prob_matrix <- scaling_factor / (count_matrix + scaling_factor)
#     
#     # 初始化一个矩阵用于存储结果
#     zero_matrix <- matrix(0, nrow = nrow(count_matrix), ncol = ncol(count_matrix))
#     
#     # 使用伯努利试验生成零值
#     for (i in 1:nrow(count_matrix)) {
#       for (j in 1:ncol(count_matrix)) {
#         zero_matrix[i, j] <- ifelse(runif(1) < prob_matrix[i, j], 0, count_matrix[i, j])
#       }
#     }
#     
#     # 计算当前的零值比例
#     zero_count <- sum(zero_matrix == 0)
#     total_elements <- length(zero_matrix)
#     current_zero_proportion <- (zero_count / total_elements) * 100
#     
#     # 调整 scaling_factor 以逼近目标比例
#     if (current_zero_proportion < target_zero_proportion) {
#       scaling_factor <- scaling_factor * 0.9  # 增加零值的概率
#     } else if (current_zero_proportion > target_zero_proportion) {
#       scaling_factor <- scaling_factor * 1.1  # 减少零值的概率
#     }
#   }
#   
#   # 第二步：精确调整零值比例
#   zero_count <- sum(zero_matrix == 0)
#   target_zero_count <- ceiling(target_zero_proportion / 100 * total_elements)
#   
#   if (zero_count < target_zero_count) {
#     # 增加零值：找到非零值并将部分置为零，直到达到目标比例
#     non_zero_indices <- which(zero_matrix != 0, arr.ind = TRUE)
#     additional_zeros_needed <- target_zero_count - zero_count
#     selected_indices <- non_zero_indices[sample(nrow(non_zero_indices), additional_zeros_needed), ]
#     zero_matrix[selected_indices] <- 0
#   } else if (zero_count > target_zero_count) {
#     # 减少零值：找到零值并将部分恢复为原始值，直到达到目标比例
#     zero_indices <- which(zero_matrix == 0, arr.ind = TRUE)
#     excessive_zeros <- zero_count - target_zero_count
#     selected_indices <- zero_indices[sample(nrow(zero_indices), excessive_zeros), ]
#     zero_matrix[selected_indices] <- count_matrix[selected_indices]
#   }
#   
#   return(zero_matrix)
# }

########## this one
# set_zero_percentage_new <- function(count_matrix, target_zero_proportion) {
#   total_elements <- length(count_matrix)
#   
#   # 调整 scaling_factor，减少零值生成
#   scaling_factor <- mean(count_matrix, na.rm = TRUE) / ((target_zero_proportion / 100) + 1)
#   
#   # 计算概率矩阵，并确保概率在 [0, 1] 范围内
#   prob_matrix <- scaling_factor / (count_matrix + scaling_factor)
#   prob_matrix <- pmin(pmax(prob_matrix, 0), 1)
#   prob_matrix[is.na(prob_matrix)] <- 0
#   
#   # 使用伯努利随机变量生成零值
#   zero_matrix <- ifelse(runif(total_elements) < prob_matrix, 0, count_matrix)
#   
#   # 确保返回的矩阵没有 NA 值
#   zero_matrix[is.na(zero_matrix)] <- 0
#   
#   return(matrix(zero_matrix, nrow = nrow(count_matrix), ncol = ncol(count_matrix)))
# }


# set_zero_percentage_new <- function(count_matrix, target_zero_proportion) {
#   total_elements <- length(count_matrix)
#   target_zero_count <- ceiling(target_zero_proportion / 100 * total_elements)
#   
#   # 初步估计 scaling_factor
#   estimated_scaling_factor <- mean(count_matrix) / (target_zero_proportion / 100)
#   
#   # 计算概率矩阵
#   prob_matrix <- estimated_scaling_factor / (count_matrix + estimated_scaling_factor)
#   
#   # 使用向量化操作生成零值
#   zero_matrix <- ifelse(runif(total_elements) < prob_matrix, 0, count_matrix)
#   
#   # 精确调整零值比例
#   zero_count <- sum(zero_matrix == 0)
#   
#   # 检查并调整零值数量
#   if (zero_count < target_zero_count) {
#     # 增加零值：找到非零值并将部分置为零，直到达到目标比例
#     non_zero_indices <- which(zero_matrix != 0)
#     additional_zeros_needed <- target_zero_count - zero_count
#     selected_indices <- sample(non_zero_indices, additional_zeros_needed)
#     zero_matrix[selected_indices] <- 0
#   } else if (zero_count > target_zero_count) {
#     # 减少零值：找到零值并将部分恢复为原始值，直到达到目标比例
#     zero_indices <- which(zero_matrix == 0)
#     excessive_zeros <- zero_count - target_zero_count
#     selected_indices <- sample(zero_indices, excessive_zeros)
#     zero_matrix[selected_indices] <- count_matrix[selected_indices]
#   }
#   
#   return(matrix(zero_matrix, nrow = nrow(count_matrix), ncol = ncol(count_matrix)))
# }

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
zero_percentages <- c(1000,100,30)

for (percentage in zero_percentages) { #percentage = 1000
  results <- list()
  
  for (i in 1:length(preprocessed_data)) { #i=5
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
    
    
    
    M_a_matrix <- set_zero_percentage_new(count_matrix=M_a_matrix, target_depth=percentage)
    M_matrix <- M_a_matrix / rowSums(M_a_matrix)
    M_matrix[is.na(M_matrix)] <- 0
    # sum(M_matrix==0)
    
    print("%:")
    print(sum(M_matrix==0)/length(M_matrix))
    # Compute M_nz_matrix after setting zero percentage
    temp <- M_matrix
    temp[temp == 0] <- 0.01
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
    
    # MEDTEST
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
  write.csv(results_df, file = paste0("sample_rarefraction_NoMe_", percentage, ".csv"), row.names = FALSE)
  
  print(paste("Results for", percentage, "% zero values have been saved"))
}

print("All results have been saved")


###########33
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Function to generate power results from a CSV file
generate_power_results <- function(file_name, zero_percentage) {
  results_df <- read.csv(file_name)
  
  # Add sample size information
  results_df$sample_size <- rep(c(50, 100, 300), each = 100)
  
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

# Generate power results for each sample 是Prevalence :power_results_30 <- generate_power_results("10_10_pre__new_sample_test_1.csv", "low")
power_results_30 <- generate_power_results("10_1_new_sample_test_1.csv", "low")
power_results_50 <- generate_power_results("10_1_new_sample_test_1000.csv", "medium")
power_results_70 <- generate_power_results("10_1_new_sample_test_100.csv", "high")

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
  scale_fill_discrete(name = "Sample Size") +
  coord_cartesian(ylim = c(0, 1)) +  # Fixed y-axis range
  facet_grid(. ~ Zero_Percentage)


ggsave(filename = "10-10-pic/sample_Prevalence_1_hist_1.1.png", plot = last_plot(), width = 10, height = 8)
################333


# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Function to generate power results from a CSV file
generate_power_results <- function(file_name, zero_percentage) {
  results_df <- read.csv(file_name)
  
  # Add sample size information
  results_df$sample_size <- rep(c(50, 100, 300), each = 100)
  
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
  
  # Add zero percentage information
  power_results$Zero_Percentage <- factor(zero_percentage, levels = c("low", "medium", "high"))
  
  return(power_results)
}

# Generate power results for each sample
power_results_low <- generate_power_results("10_1_new_sample_test_1.csv", "low")
power_results_medium <- generate_power_results("10_1_new_sample_test_1000.csv", "medium")
power_results_high <- generate_power_results("10_1_new_sample_test_100.csv", "high")

# Combine all power results
all_power_results <- bind_rows(power_results_low, power_results_medium, power_results_high)

# Plot power analysis as a line plot with different shapes for each Zero_Percentage level
ggplot(all_power_results, aes(x = sample_size, y = power, group = interaction(Method, Zero_Percentage), color = Method)) +
  geom_line(aes(linetype = Zero_Percentage), size = 1) +
  geom_point(aes(shape = Zero_Percentage), size = 3) +
  theme_minimal() +
  labs(title = "Power Analysis by Method and Zero Value Percentage",
       x = "Sample Size",
       y = "Power") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_discrete(name = "Method") +
  scale_linetype_manual(name = "Zero Percentage", values = c("solid", "dashed", "dotted"), labels = c("Low", "Medium", "High")) +
  scale_shape_manual(name = "Zero Percentage", values = c(16, 17, 18), labels = c("Low", "Medium", "High")) +
  coord_cartesian(ylim = c(0, 1))  # Fixed y-axis range




# 保存图片到新文件夹
ggsave(filename = "10-10-pic/sample_1_line.png", plot = last_plot(), width = 10, height = 8)
##########333
###############3
################

############### type-1 error
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Function to generate Type-1 error results from a CSV file
generate_type1_error_results <- function(file_name, zero_percentage) {
  results_df <- read.csv(file_name)
  
  # Add sample size information
  results_df$sample_size <- rep(c(50, 100, 300), each = 80)
  
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
  type1_error_results$Zero_Percentage <- factor(zero_percentage, levels = c("Low", "Medium", "High"))
  
  return(type1_error_results)
}

# Generate Type-1 error results for each sample

type1_error_results_30 <- generate_type1_error_results("sample_rarefraction_NoMe_1000.csv", "Low")
type1_error_results_50 <- generate_type1_error_results("sample_rarefraction_NoMe_100.csv", "Medium")
type1_error_results_70 <- generate_type1_error_results("sample_rarefraction_NoMe_30.csv", "High")

# Combine all Type-1 error results
all_type1_error_results <- bind_rows( type1_error_results_30, type1_error_results_50, type1_error_results_70)

# Plot Type-1 error analysis with facet_grid
ggplot(all_type1_error_results, aes(x = Method, y = type1_error, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Type-1 Error Analysis by Zero Value Percentage",
       x = "Method",
       y = "Type-1 Error") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample Size") +
  coord_cartesian(ylim = c(0, 1)) +  # Fixed y-axis range
  facet_grid(. ~ Zero_Percentage)


