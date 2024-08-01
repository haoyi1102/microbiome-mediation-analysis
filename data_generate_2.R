rm(list = ls())
library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)

# Load the generate_simulation_data function from data_generate_1.R
source("data_generate_1.R")


library(ggplot2)

# 定义生成零值占比的函数
# i=0
# generate_zero_proportion <- function() {
#   result <- generate_simulation_data(
#     n_sample = 100, template = "Vaginal", n_feature = 20,
#     metadata_effect_size = 0, perc_feature_spiked_metadata = 1,
#     median_read_depth = 10000, alpha_0 = 1, alpha_T = 0.5,
#     alpha_M_value = 0.2, noise_sd = 0.5
#   )
#   zero_count <- sum(result$absolute_M == 0)
#   total_elements <- length(result$absolute_M)
#   zero_proportion <- zero_count / total_elements
#   i=i+1
#   print(i)
#   return(zero_proportion)
# }
# 
# # 生成1000次模拟数据并收集零值占比
# 
# zero_proportions <- replicate(1000, generate_zero_proportion())
# 
# # 绘制零值占比的分布柱状图
# zero_proportion_df <- data.frame(zero_proportion = zero_proportions)
# ggplot(zero_proportion_df, aes(x = zero_proportion)) +
#   geom_histogram(binwidth = 0.01, fill = "red", color = "black", alpha = 0.7) +
#   labs(title = "Distribution of Zero Proportion in Simulated Data (vaginal)",
#        x = "Zero Proportion",
#        y = "Frequency") +
#   theme_minimal()

# 直接拼数据可行 替换feature的时候要注意 是不是起作用的feature
# 没有哪一个变量直接改变0 metadata_effect_size和perc_feature_spiked_metadata都不行
# # 
# # # Print the results
# # print(result)

# Create a new directory to save the results


# Function to save individual results to .rds file
save_result <- function(result, variable, value, index) {
  rds_file <- file.path(output_dir, paste0("result_", variable, "_", value, "_", index, ".rds"))
  saveRDS(result, file = rds_file)
}
load("fitted_model_1.RData")
# Fixed parameters
template <- fitted
alpha_0 <- 1
alpha_T <- 0.5
alpha_M_value <- 0.2

# Parameter values to iterate over
n_features <- c(10, 20, 40)
n_samples <- c(50, 100, 300)
metadata_effect_sizes <- c(0.5, 1, 10)
perc_feature_spiked_metadatas <- c(0.1, 0.2, 0.5)
median_read_depths <- c(100, 1000, 10000)
noise_sds <- c(0.5, 1, 3)

# # Generate 100 datasets for each sample size and save the results
# output_dir <- "simulation_data_sample"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir)
# }
# for (n_sample in n_samples) {
#   for (i in 1:20) {
#     # Parameters defined within the loop
#     n_feature <- 20
#     metadata_effect_size <- 1
#     perc_feature_spiked_metadata <- 0.1
#     median_read_depth <- 10000
#     noise_sd <- 0.5
#     
#     result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = n_feature, 
#                                        metadata_effect_size = metadata_effect_size, perc_feature_spiked_metadata = perc_feature_spiked_metadata, 
#                                        median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T, 
#                                        alpha_M_value = alpha_M_value, noise_sd = noise_sd)
#     save_result(result, "n_sample", n_sample, i)
#   }
# }
# 
# print("All results have been saved successfully.")
###########33333
output_dir <- "simulation_data_sample_new"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
for (n_sample in n_samples) {
  for (i in 1:20) {
    # Parameters defined within the loop
    n_feature <- 20
    metadata_effect_size <- 8
    perc_feature_spiked_metadata <- 0.2
    median_read_depth <- 100000
    noise_sd <- 0.5
    
    result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = n_feature, 
                                       metadata_effect_size = metadata_effect_size, perc_feature_spiked_metadata = perc_feature_spiked_metadata, 
                                       median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                       alpha_M_value = alpha_M_value, noise_sd = noise_sd)
    save_result(result, "n_sample", n_sample, i)
  }
}

print("All results have been saved successfully.")

# no_mediation_effect
output_dir <- "simulation_data_sample_new_no_mediation"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (n_sample in n_samples) {
  for (i in 1:10) {
    # Parameters defined within the loop i=1 n_sample=100
    n_feature <- 20
    metadata_effect_size <- 1
    perc_feature_spiked_metadata <- 0.1
    median_read_depth <- 10000
    noise_sd <- 0.5
    
    result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = n_feature, 
                                       metadata_effect_size = metadata_effect_size, perc_feature_spiked_metadata = perc_feature_spiked_metadata, 
                                       median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                       alpha_M_value = alpha_M_value, noise_sd = noise_sd)
    
    # Generate no mediation effect dataset by permuting the mediator matrix
    result_no_mediation <- result
    col_names <- colnames(result$absolute_M)  # Save column names
    result_no_mediation$absolute_M <- result$absolute_M[, sample(ncol(result$absolute_M))]
    colnames(result_no_mediation$absolute_M) <- col_names  # Restore column names
    save_result(result_no_mediation, "n_sample", n_sample, i)  }
}

print("All results have been saved successfully.")
############
output_dir <- "simulation_data_sample_no_mediation"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (n_sample in n_samples) {
  for (i in 1:20) {
    # Parameters defined within the loop i=1 n_sample=100
    n_feature <- 20
    metadata_effect_size <- 1
    perc_feature_spiked_metadata <- 0.1
    median_read_depth <- 10000
    noise_sd <- 0.5
    
    result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = n_feature, 
                                       metadata_effect_size = metadata_effect_size, perc_feature_spiked_metadata = perc_feature_spiked_metadata, 
                                       median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                       alpha_M_value = alpha_M_value, noise_sd = noise_sd)
    
    # Generate no mediation effect dataset by permuting the mediator matrix
    result_no_mediation <- result
    col_names <- colnames(result$absolute_M)  # Save column names
    result_no_mediation$absolute_M <- result$absolute_M[, sample(ncol(result$absolute_M))]
    colnames(result_no_mediation$absolute_M) <- col_names  # Restore column names
    save_result(result_no_mediation, "n_sample", n_sample, i)  }
}

print("All results have been saved successfully.")

# Generate 100 datasets for each feature size and save the results
output_dir <- "simulation_data_feature"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
for (n_feature in n_features) {
  for (i in 1:20) {
    # Parameters defined within the loop
    n_sample <- 100
    metadata_effect_size <- 1
    perc_feature_spiked_metadata <- 0.1
    median_read_depth <- 10000
    noise_sd <- 0.5
    
    result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = n_feature, 
                                       metadata_effect_size = metadata_effect_size, perc_feature_spiked_metadata = perc_feature_spiked_metadata, 
                                       median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                       alpha_M_value = alpha_M_value, noise_sd = noise_sd)
    save_result(result, "n_feature", n_feature, i)
  }
}

print("All results have been saved successfully.")

# no_mediation_effect
output_dir <- "simulation_data_feature_no_mediation"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (n_feature in n_features) {
  for (i in 1:20) {
    # Parameters defined within the loop i=1 n_sample=100
    n_sample <- 100
    metadata_effect_size <- 1
    perc_feature_spiked_metadata <- 0.1
    median_read_depth <- 10000
    noise_sd <- 0.5
    
    result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = n_feature, 
                                       metadata_effect_size = metadata_effect_size, perc_feature_spiked_metadata = perc_feature_spiked_metadata, 
                                       median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                       alpha_M_value = alpha_M_value, noise_sd = noise_sd)
    
    # Generate no mediation effect dataset by permuting the mediator matrix
    result_no_mediation <- result
    col_names <- colnames(result$absolute_M)  # Save column names
    result_no_mediation$absolute_M <- result$absolute_M[, sample(ncol(result$absolute_M))]
    colnames(result_no_mediation$absolute_M) <- col_names  # Restore column names
    save_result(result_no_mediation, "n_feature", n_feature, i)  }
}

print("All results have been saved successfully.")





















