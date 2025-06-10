# 设置工作目录和加载所需包
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

source("modima.R")
source("MedOmniTest.R")
library(SparseMCMM)
library(ccmm)
library(LDM)
library(MicroBVS)
library(MedZIM)
source("highmed2019.r")
source('fromSKAT.R')

set.seed(1234)

# 读取生成的数据
output_dir <- "10-10-new_method_simulation_data_sample_new"
repeat_count <- 100

# 列出目录中的所有 .rds 文件
rds_files <- list.files(output_dir, pattern = "*.rds", full.names = TRUE)
rds_files <- mixedsort(rds_files)

# 定义处理每个结果的函数
preprocess_single_result <- function(result) {
  Y_vector <- as.numeric(result$Y)
  T_vector <- as.numeric(result$T)
  M_a_matrix <- as.matrix(t(result$absolute_M))
  colnames(M_a_matrix) <- paste0("taxon_", 1:ncol(M_a_matrix))
  
  matrices <- compute_matrices(M_a_matrix)
  
  list(Y_vector = Y_vector, T_vector = T_vector, M_a_matrix = M_a_matrix, M_matrix = matrices$M_matrix, M_nz_matrix = matrices$M_nz_matrix, feature = result$feature)
}

compute_matrices <- function(M_a_matrix) {
  M_matrix <- M_a_matrix / rowSums(M_a_matrix)
  
  temp <- M_a_matrix
  temp[temp == 0] <- 0.5
  M_nz_matrix <- temp / rowSums(temp)
  
  list(M_matrix = M_matrix, M_nz_matrix = M_nz_matrix)
}

preprocess_results <- function(results) {
  lapply(results, preprocess_single_result)
}

split_results <- function(files) {
  results <- list()
  for (file in files) {
    result <- readRDS(file)
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

# 分割结果文件
split_results <- split_results(rds_files)

# 对每个样本数据进行预处理
preprocessed_results_by_sample <- lapply(split_results$n, preprocess_results)
preprocessed_data <- preprocessed_results_by_sample[["sample"]]

results_low <- list()
results_medium <- list()
results_high <- list()

for (i in 1:length(preprocessed_data)) {
  loop_start_time <- Sys.time()  # 开始计时
  
  data_list <- preprocessed_data[i]
  data <- data_list[[1]]
  
  # 提取向量和矩阵
  Y_vector <- as.vector(data$Y_vector)
  T_vector <- as.vector(data$T_vector)
  M_matrix <- as.matrix(data$M_matrix)
  M_a_matrix <- as.matrix(data$M_a_matrix)
  
  M_matrix <- M_a_matrix / rowSums(M_a_matrix)
  print(sum(M_matrix == 0) / length(M_matrix))
  
  temp <- M_matrix
  temp[temp == 0] <- 0.5
  M_nz_matrix <- temp / rowSums(temp)
  
  # MODIMA 分析
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
  
  # MEDTEST 分析
  medtest_time <- system.time({
    medtest_p_value <- tryCatch({
      m_list <- list(euc = dist(M_matrix))
      medtest_result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)
      medtest_result$permP
    }, error = function(e) {
      print(paste("Error in MEDTEST:", e$message))
      1
    })
  })
  print(paste("MEDTEST completed in:", medtest_time[3], "seconds"))
  
  # CCMM 分析
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
  
  # LDM-MED 分析
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
  
  # PERMANOVA-MED 分析
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
  
  loop_end_time <- Sys.time()  # 结束计时
  print(paste("Loop", i, "completed in:", loop_end_time - loop_start_time, "mins"))
  
 

# 将结果保存为 CSV 文件
if (i <= repeat_count || (i > 3 * repeat_count && i <= 4 * repeat_count) || (i > 6 * repeat_count && i <= 7 * repeat_count)) {
  results_low[[i]] <- list(
    modima_p_value = modima_p_value,
    medtest_p_value = medtest_p_value,
    ccmm_CI_lower = ccmm_CI_lower,
    ccmm_CI_upper = ccmm_CI_upper,
    ldm_p_value = ldm_p_value,
    permanova_p_value = permanova_p_value
  )
} else if ((i > repeat_count && i <= 2 * repeat_count) || (i > 4 * repeat_count && i <= 5 * repeat_count) || (i > 7 * repeat_count && i <= 8 * repeat_count)) {
  results_medium[[i]] <- list(
    modima_p_value = modima_p_value,
    medtest_p_value = medtest_p_value,
    ccmm_CI_lower = ccmm_CI_lower,
    ccmm_CI_upper = ccmm_CI_upper,
    ldm_p_value = ldm_p_value,
    permanova_p_value = permanova_p_value
  )
} else if ((i > 2 * repeat_count && i <= 3 * repeat_count) || (i > 5 * repeat_count && i <= 6 * repeat_count) || (i > 8 * repeat_count && i <= 9 * repeat_count)) {
  results_high[[i]] <- list(
    modima_p_value = modima_p_value,
    medtest_p_value = medtest_p_value,
    ccmm_CI_lower = ccmm_CI_lower,
    ccmm_CI_upper = ccmm_CI_upper,
    ldm_p_value = ldm_p_value,
    permanova_p_value = permanova_p_value
  )
}}

write.csv(do.call(rbind, lapply(results_low, as.data.frame)), file = "low_results.csv", row.names = FALSE)
write.csv(do.call(rbind, lapply(results_medium, as.data.frame)), file = "medium_results.csv", row.names = FALSE)
write.csv(do.call(rbind, lapply(results_high, as.data.frame)), file = "high_results.csv", row.names = FALSE)
print("All results have been saved")


#######################################################
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

# Generate power results for each sample
power_results_30 <- generate_power_results("low_results.csv", "low")
power_results_50 <- generate_power_results("medium_results.csv", "medium")
power_results_70 <- generate_power_results("high_results.csv", "high")

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

ggsave(filename = "10-10-pic/sample_2_hist.png", plot = last_plot(), width = 10, height = 8)







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
power_results_low <- generate_power_results("low_results.csv", "low")
power_results_medium <- generate_power_results("medium_results.csv", "medium")
power_results_high <- generate_power_results("high_results.csv", "high")

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
ggsave(filename = "10-10-pic/sample_2_line.png", plot = last_plot(), width = 10, height = 8)

