# 清理工作空间并加载必要的库
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

# 加载自定义函数和所需的源文件
source("modima.R")
source("MedOmniTest.R")
library(SparseMCMM)
library(ccmm)
library(LDM)
library(MicroBVS)
library(MedZIM)

# Function to preprocess a single result
preprocess_single_result <- function(result) {
  Y_vector <- as.numeric(result$Y)
  T_vector <- as.numeric(result$T)
  M_a_matrix <- as.matrix(t(result$absolute_M))
  colnames(M_a_matrix) <- paste0("taxon_", 1:ncol(M_a_matrix))
  
  # Extract M_matrix and M_nz_matrix using a separate function
  matrices <- compute_matrices(M_a_matrix)
  
  list(Y_vector = Y_vector, T_vector = T_vector, M_a_matrix = M_a_matrix, 
       M_matrix = matrices$M_matrix, M_nz_matrix = matrices$M_nz_matrix, feature = result$feature)
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

# Split the results into categories based on filenames
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

# 设置结果保存的主目录
main_output_dir <- "simulation_data_alpha_beta_results"
if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir)
}

# 读取每个alpha_beta文件夹中的数据文件
alpha_beta_folders <- list.dirs(path = "simulation_data_alpha_beta", full.names = TRUE, recursive = FALSE)

# 遍历每个alpha_beta文件夹进行分析
for (folder in alpha_beta_folders) { # folder = "simulation_data_alpha_beta/alpha_0_beta_0"  
  
  # 获取当前文件夹的名称并为该文件夹创建结果子目录
  folder_name <- basename(folder)
  result_subfolder <- file.path(main_output_dir, folder_name)
  
  # 如果结果文件夹已经存在，并且包含结果文件，则跳过这个文件夹
  if (dir.exists(result_subfolder) && length(list.files(result_subfolder, pattern = "*.csv", full.names = TRUE)) > 0) {
    print(paste("Skipping folder", folder_name, "as results already exist."))
    next
  }
  
  if (!dir.exists(result_subfolder)) {
    dir.create(result_subfolder)
  }
  
  # 列出当前文件夹中的所有 .rds 文件
  rds_files <- list.files(folder, pattern = "*.rds", full.names = TRUE)
  rds_files <- mixedsort(rds_files)
  
  # 使用 split_results 将文件组织成不同类别
  split_results_list <- split_results(rds_files)
  
  # 初始化一个列表用于存储所有结果
  results <- list()
  
  # 遍历所有的 .rds 文件并进行分析
  for (i in seq_along(rds_files)) { # i =90
    loop_start_time <- Sys.time()  # 记录循环开始时间
    
    result <- readRDS(rds_files[i])
    preprocessed_data <- preprocess_single_result(result)  # 使用 preprocess_single_result 进行预处理
    
    # 提取必要的向量和矩阵
    Y_vector <- as.vector(preprocessed_data$Y_vector)
    T_vector <- as.vector(preprocessed_data$T_vector)
    M_matrix <- as.matrix(preprocessed_data$M_matrix)
    M_nz_matrix <- as.matrix(preprocessed_data$M_nz_matrix)
    
    # MODIMA 分析
    modima_p_value <- tryCatch({
      T_dist <- dist(T_vector)
      Y_dist <- dist(Y_vector)
      M_dist <- dist(M_matrix)
      modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
      modima_result$p.value
    }, error = function(e) {
      print(paste("Error in MODIMA:", e$message))
      NA
    })
    
    # MEDTEST 分析
    medtest_p_value <- tryCatch({
      m_list <- list(euc = dist(M_matrix))
      medtest_result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)
      medtest_result$permP
    }, error = function(e) {
      print(paste("Error in MEDTEST:", e$message))
      NA
    })
    
    # CCMM 分析
    ccmm_CI_lower <- ccmm_CI_upper <- NA
    tryCatch({
      ccmm_result <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
                          sig.level = 0.05, tol = 1e-06, max.iter = 5000)
      ccmm_CI_lower <- as.numeric(ccmm_result$TIDE.CI["2.5%"])
      ccmm_CI_upper <- as.numeric(ccmm_result$TIDE.CI["97.5%"])
    }, error = function(e) {
      print(paste("Error in CCMM:", e$message))
    })
    
    # LDM-MED 分析
    ldm_p_value <- tryCatch({
      data_ldm <- data.frame(Y = Y_vector, T = T_vector, M_nz_matrix = M_nz_matrix)
      ldm_result <- ldm(formula = M_nz_matrix ~ T + Y, data = data_ldm, seed = 1234, test.mediation = TRUE)
      ldm_result$med.p.global.omni
    }, error = function(e) {
      print(paste("Error in LDM-MED:", e$message))
      NA
    })
    
    # PERMANOVA-MED 分析
    permanova_p_value <- tryCatch({
      data_ldm <- data.frame(Y = Y_vector, T = T_vector)
      formula <- as.formula("M_nz_matrix ~ T + Y")
      res_perm <- permanovaFL(formula = formula, data = data_ldm, dist.method = "bray", seed = 1234,
                              n.perm.max = 1000, n.cores = 4, test.mediation = TRUE, verbose = TRUE)
      res_perm$med.p.permanova
    }, error = function(e) {
      print(paste("Error in PERMANOVA-MED:", e$message))
      NA
    })
    
    loop_end_time <- Sys.time()  # 记录循环结束时间
    loop_duration <- loop_end_time - loop_start_time
    print(paste("Loop", i, "completed in:", loop_duration, "seconds"))
    
    # 将分析结果存储到列表中
    results[[i]] <- list(
      modima_p_value = modima_p_value,
      medtest_p_value = medtest_p_value,
      ccmm_CI_lower = ccmm_CI_lower,
      ccmm_CI_upper = ccmm_CI_upper,
      ldm_p_value = ldm_p_value,
      permanova_p_value = permanova_p_value
    )
  }
  
  # 将结果列表转换为数据框并保存为CSV文件
  results_df <- do.call(rbind, lapply(results, function(x) {
    c(modima_p_value = x$modima_p_value,
      medtest_p_value = x$medtest_p_value,
      ccmm_CI_lower = x$ccmm_CI_lower,
      ccmm_CI_upper = x$ccmm_CI_upper,
      ldm_p_value = x$ldm_p_value,
      permanova_p_value = x$permanova_p_value)
  }))
  
  results_df <- as.data.frame(results_df)
  output_file <- file.path(result_subfolder, paste0(folder_name, "_results.csv"))
  write.csv(results_df, file = output_file, row.names = FALSE)
  
  print(paste("Results for folder", folder_name, "have been saved successfully."))
}

print("All analyses completed and results saved.")














##########################33333

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the results_long_df from the CSV file
results_long_df <- read.csv("results_long_df.csv")

# Create a vector of effect sizes
effect_sizes <- rep(c(1, seq(5, 20, by = 5)), each = 100)# each controls sample

# Map iterations to effect sizes
results_long_df <- results_long_df %>%
  mutate(effect_size = effect_sizes[iteration])

# Remove unused columns
# results_long_df <- results_long_df %>%
#   select(method, iteration, true_positives, false_positives, effect_size)

results_long_df <- dplyr::select(results_long_df, method, iteration, true_positives, false_positives, effect_size)


# Group by effect_size and method, then summarize the metrics
summary_df <- results_long_df %>%
  group_by(effect_size, method) %>%
  summarise(
    true_positives = sum(true_positives, na.rm = TRUE),
    false_positives = sum(false_positives, na.rm = TRUE),
  )

# Calculate false_negatives based on total_positives and true_positives
summary_df <- summary_df %>%
  mutate(false_negatives = 0.2*20*100 - true_positives)

# # Calculate precision, recall, and f1_score
# summary_df <- summary_df %>%
#   mutate(
#     precision = ifelse(true_positives + false_positives > 0, true_positives / (true_positives + false_positives), 0),
#     recall = ifelse(true_positives + false_negatives > 0, true_positives / (true_positives + false_negatives), 0),
#     f1_score = ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), 0)
#   )
# 
# results_long_df_long <- summary_df %>%
#   pivot_longer(cols = c(precision, recall, f1_score), names_to = "metric", values_to = "value")
# 
# ggplot(results_long_df_long, aes(x = effect_size, y = value, color = method, linetype = metric)) +
#   geom_line(size = 0.8) +  # Adjust line size to be thinner
#   geom_point(size = 2) +
#   labs(title = "Performance Metrics vs Effect Size",
#        x = "Effect Size",
#        y = "Metric Value",
#        color = "Method",
#        linetype = "Metric") +
#   theme_minimal() +
#   theme(legend.position = "right",  # Move legend to the right
#         legend.title = element_text(size = 12),  # Increase legend title size
#         legend.text = element_text(size = 10),  # Increase legend text size
#         plot.title = element_text(hjust = 0.5, size = 16),  # Center and increase title size
#         axis.title = element_text(size = 14),  # Increase axis title size
#         axis.text = element_text(size = 12))  # Increase axis text size

# Calculate precision, recall, and f1_score
summary_df <- summary_df %>%
  mutate(
    precision = ifelse(true_positives + false_positives > 0, true_positives / (true_positives + false_positives), 0),
    recall = ifelse(true_positives + false_negatives > 0, true_positives / (true_positives + false_negatives), 0),
    f1_score = ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), 0)
  )

# Transform data to long format for plotting
results_long_df_long <- summary_df %>%
  pivot_longer(cols = c(precision, recall, f1_score), names_to = "metric", values_to = "value")


library(openxlsx)
write.xlsx(summary_df, file = "results_summary.xlsx")


# Plotting precision
ggplot(results_long_df_long %>% filter(metric == "precision"), aes(x = effect_size, y = value, color = method)) +
  geom_line(size = 0.8) +  # Adjust line size to be thinner
  geom_point(size = 2) +
  labs(title = "Precision vs Effect Size",
       x = "Effect Size",
       y = "Precision",
       color = "Method") +
  theme_minimal() +
  theme(legend.position = "right",  # Move legend to the right
        legend.title = element_text(size = 12),  # Increase legend title size
        legend.text = element_text(size = 10),  # Increase legend text size
        plot.title = element_text(hjust = 0.5, size = 16),  # Center and increase title size
        axis.title = element_text(size = 14),  # Increase axis title size
        axis.text = element_text(size = 12))  # Increase axis text size

# Plotting recall
ggplot(results_long_df_long %>% filter(metric == "recall"), aes(x = effect_size, y = value, color = method)) +
  geom_line(size = 0.8) +  # Adjust line size to be thinner
  geom_point(size = 2) +
  labs(title = "Recall vs Effect Size",
       x = "Effect Size",
       y = "Recall",
       color = "Method") +
  theme_minimal() +
  theme(legend.position = "right",  # Move legend to the right
        legend.title = element_text(size = 12),  # Increase legend title size
        legend.text = element_text(size = 10),  # Increase legend text size
        plot.title = element_text(hjust = 0.5, size = 16),  # Center and increase title size
        axis.title = element_text(size = 14),  # Increase axis title size
        axis.text = element_text(size = 12))  # Increase axis text size

# Plotting f1_score
ggplot(results_long_df_long %>% filter(metric == "f1_score"), aes(x = effect_size, y = value, color = method)) +
  geom_line(size = 0.8) +  # Adjust line size to be thinner
  geom_point(size = 2) +
  labs(title = "F1 Score vs Effect Size",
       x = "Effect Size",
       y = "F1 Score",
       color = "Method") +
  theme_minimal() +
  theme(legend.position = "right",  # Move legend to the right
        legend.title = element_text(size = 12),  # Increase legend title size
        legend.text = element_text(size = 10),  # Increase legend text size
        plot.title = element_text(hjust = 0.5, size = 16),  # Center and increase title size
        axis.title = element_text(size = 14),  # Increase axis title size
        axis.text = element_text(size = 12))  # Increase axis text size











# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Function to generate power results from a CSV file
generate_power_results <- function(file_name, zero_percentage) {
  results_df <- read.csv(file_name)
  
  # Add noise size information
  results_df$noise_size <- rep(c(0, 0.25,0.5,0.75,1) , each = 10)
  
  # Replace missing values with 1 (acceptance)
  results_df[is.na(results_df)] <- 1
  
  # Rename methods for better readability
  results_long <- melt(results_df, id.vars = "noise_size", 
                       measure.vars = c("modima_p_value", "medtest_p_value", "ldm_p_value", "permanova_p_value"),
                       variable.name = "Method", value.name = "p_value")
  
  results_long$Method <- recode(results_long$Method,
                                "modima_p_value" = "MODIMA",
                                "medtest_p_value" = "MedTest",
                                "ldm_p_value" = "LDM",
                                "permanova_p_value" = "PERMANOVA")
  
  # Calculate power (p-value < 0.05)
  power_results <- results_long %>%
    group_by(noise_size, Method) %>%
    summarise(power = mean(p_value < 0.05, na.rm = TRUE))
  
  # Calculate rejection rate based on confidence intervals (reject if not containing 0)
  results_df <- results_df %>%
    mutate(ccmm_reject = ifelse(ccmm_CI_lower > 0 | ccmm_CI_upper < 0, 1, 0))
  
  # Add ccmm_reject to power_results
  ccmm_power_results <- results_df %>%
    group_by(noise_size) %>%
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
power_results_30 <- generate_power_results("effect_test_.csv", "30% Zero Values")

# Combine all power results
all_power_results <- power_results_30

# Plot power analysis with line chart
ggplot(all_power_results, aes(x = noise_size, y = power, color = Method, group = Method)) +
  geom_line(size = 0.8) +  # Make lines thinner
  geom_point(size = 1.5) +  # Smaller points
  theme_minimal() +
  labs(title = "Power Analysis by Noise Size",
       x = "Noise Size",
       y = "Power") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",  # Move legend to bottom
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 10),  # Adjust legend text size
        plot.title = element_text(hjust = 0.5)) +  # Center the plot title
  scale_color_brewer(palette = "Set1") +  # Use a color palette for better distinction
  coord_cartesian(ylim = c(0, 1)) +  # Fixed y-axis range
  facet_grid(. ~ Zero_Percentage)