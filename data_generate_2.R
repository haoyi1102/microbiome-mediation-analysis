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

# Fixed parameters
template <-  readRDS("sparsedossa_baseline.rds")
template_medium <-  readRDS("sparsedossa_baseline_medium.rds")
template_high <-  readRDS("sparsedossa_baseline_high.rds")
alpha_0 <- 1
alpha_T <- 0.5
alpha_M_value <- 0.2

# Parameter values to iterate over
n_features <- c(10, 20, 40)
n_samples <- c(50, 100, 300)
metadata_effect_sizes <- c(0.5, 1, 10)
perc_feature_spiked_metadatas <- c(0.1, 0.2, 0.5)
median_read_depths <- c(1000, 10000, 100000)
noise_sds <- c(0.5, 1, 1.5,2,2.5,3,3.5,4,4.5,5)

########################33
output_dir <- "11_18_simulation_data_per_feature_valid_new"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop to define parameters and run simulations
for (num_spiked_features in seq(1, 10)) {  # Spiked features from 1 to 50% of total features (assuming 20 total)
  for (i in 1:100) {  # Number of replicates
    # Define parameters inside the loop
    n_sample <- 100  # Fixed sample size
    n_features <- 50  # Total feature size
    perc_feature_spiked_metadata <- num_spiked_features / n_features  # Percentage of spiked features
    metadata_effect_size <- 8  # Metadata effect size
    median_read_depth <- 100000  # Median read depth
    alpha_0 <- 1  # Fixed alpha_0
    alpha_T <- 0.5  # Fixed alpha_T
    alpha_M_value <- 0.2  # Fixed alpha_M_value
    noise_sd <- 0.5  # Noise standard deviation
    
    # Generate simulation data
    result <- generate_simulation_data(
      n_sample = n_sample,
      template = template,
      n_feature = n_features,
      metadata_effect_size = metadata_effect_size,
      perc_feature_spiked_metadata = perc_feature_spiked_metadata,
      median_read_depth = median_read_depth,
      alpha_0 = alpha_0,
      alpha_T = alpha_T,
      alpha_M_value = alpha_M_value,
      noise_sd = noise_sd
    )
    
    save_result(result, "spiked_features", num_spiked_features, i)
  }
}

print("All results have been saved successfully.")
#####################


# #######################333 10-10
# output_dir <- "10-10-new_method_simulation_data_sample_new"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir)
# }
# 
# 
# templates <- list(template, template_medium, template_high)
# 
# 
# for (n_sample in n_samples) {
#   iteration_count <- 1 
#   
#   for (template_index in 1:length(templates)) {
#     current_template <- templates[[template_index]]
#     
#     for (i in 1:100) {  
#       # 模拟参数设置
#       n_feature <- 20
#       metadata_effect_size <- 8
#       perc_feature_spiked_metadata <- 0.2
#       median_read_depth <- 100000
#       noise_sd <- 0.5
#       
#       # 生成模拟数据
#       result <- generate_simulation_data(
#         n_sample = n_sample,
#         template = current_template,
#         n_feature = n_feature,
#         metadata_effect_size = metadata_effect_size,
#         perc_feature_spiked_metadata = perc_feature_spiked_metadata,
#         median_read_depth = median_read_depth,
#         alpha_0 = alpha_0,
#         alpha_T = alpha_T,
#         alpha_M_value = alpha_M_value,
#         noise_sd = noise_sd
#       )
#       
#       # 保存结果，使用编号进行命名
#       save_result(result, "n_sample", n_sample, iteration_count)
#       
#       # 增加编号计数器
#       iteration_count <- iteration_count + 1
#     }
#   }
# }
# 
# print("All results have been saved successfully.")
# 
# 

######################3

###########33333
output_dir <- "10-10-pervenlence-new method simulation_data_sample_new"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
for (n_sample in n_samples) {
  for (i in 1:100) {
    # Parameters defined within the loop n_sample= 100 set.seed(1253) i = 1
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
  
}}

print("All results have been saved successfully.")

# no_mediation_effect
output_dir <- "simulation_data_sample_new_no_mediation"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (n_sample in n_samples) {
  for (i in 1:80) {
    # Parameters defined within the loop i=1 n_sample=100
    n_feature <- 20
    metadata_effect_size <- 8
    perc_feature_spiked_metadata <- 0.2
    median_read_depth <- 100000
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
    metadata_effect_size <- 8
    perc_feature_spiked_metadata <- 0.2
    median_read_depth <- 100000
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
output_dir <- "new_simulation_data_feature"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
for (n_feature in n_features) {
  for (i in 1:100) {
    # Parameters defined within the loop
    n_sample <- 100
    metadata_effect_size <- 8
    perc_feature_spiked_metadata <- 0.2
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
  for (i in 1:100) {
    # Parameters defined within the loop i=1 n_sample=100
    n_sample <- 100
    metadata_effect_size <- 8
    perc_feature_spiked_metadata <- 0.2
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






output_dir <- "new method simulation_data_noise"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
for (noise_sd in noise_sds) {
  for (i in 1:100) {
    # Parameters defined within the loop
    n_sample <- 100
    n_feature <- 20
    metadata_effect_size <- 8
    perc_feature_spiked_metadata <- 0.2
    median_read_depth <- 10000
    
    
    result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = n_feature,
                                       metadata_effect_size = metadata_effect_size, perc_feature_spiked_metadata = perc_feature_spiked_metadata,
                                       median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T,
                                       alpha_M_value = alpha_M_value, noise_sd = noise_sd)
    save_result(result, "noise_sd", noise_sd, i)
  }
}

print("All results have been saved successfully.")




output_dir <- "simulation_data_noise_no_mediation"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (noise_sd in noise_sds) {
  for (i in 1:100) {
    # Parameters defined within the loop
    n_sample <- 100
    n_feature <- 20
    metadata_effect_size <- 8
    perc_feature_spiked_metadata <- 0.2
    median_read_depth <- 10000
    
    result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = n_feature, 
                                       metadata_effect_size = metadata_effect_size, perc_feature_spiked_metadata = perc_feature_spiked_metadata, 
                                       median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                       alpha_M_value = alpha_M_value, noise_sd = noise_sd)
    
    # Generate no mediation effect dataset by permuting the mediator matrix
    result_no_mediation <- result
    col_names <- colnames(result$absolute_M)  # Save column names
    result_no_mediation$absolute_M <- result$absolute_M[, sample(ncol(result$absolute_M))]
    colnames(result_no_mediation$absolute_M) <- col_names  # Restore column names
    save_result(result, "noise_sd", noise_sd, i)
  }
}

print("All results have been saved successfully.")




# 创建一个新的输出目录
output_dir <- "new_method_simulation_data_effect_size"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 设定不同的效应大小
effect_sizes <- c(0, 5)  # 这里可以根据需要设置五个不同的效应大小

# 循环每个效应大小并生成模拟数据
for (effect_size in effect_sizes) {
  for (i in 1:10) {
    # 固定参数设置
    n_sample <- 100
    n_feature <- 20
    perc_feature_spiked_metadata <- 0.2
    median_read_depth <- 10000
    
    # 调用生成模拟数据的函数
    result <- generate_simulation_data(
      n_sample = n_sample,
      template = template,
      n_feature = n_feature,
      metadata_effect_size = effect_size,  # 使用不同的效应大小
      perc_feature_spiked_metadata = perc_feature_spiked_metadata,
      median_read_depth = median_read_depth,
      alpha_0 = alpha_0,
      alpha_T = alpha_T,
      alpha_M_value = alpha_M_value,
      noise_sd = 0.5  # 噪声水平固定为0.5
    )
    
    # 保存结果
    save_result(result, "effect_size", effect_size, i)
  }
}

print("所有结果已成功保存。")

#################################3 test

# 总文件夹路径
main_output_dir <- "simulation_data_alpha_beta"

# 如果总文件夹不存在，则创建
if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir)
}

# 设定参数
effect_sizes <- c(0, 2, 5, 10)  # Alpha值
alpha_M_values <- c(0, 0.25, 0.5, 1)  # Beta值
sample_sizes <- c(20, 50, 100, 200)  # 样本大小
num_datasets_per_combination <- 50  # 每种组合生成的数据集数量

# 生成模拟数据并保存到相应文件夹
for (effect_size in effect_sizes) {
  for (alpha_M_value in alpha_M_values) {
    
    # 创建以 alpha 和 beta 值命名的子文件夹
    subfolder_name <- paste0("alpha_", effect_size, "_beta_", alpha_M_value)
    subfolder_path <- file.path(main_output_dir, subfolder_name)
    if (!dir.exists(subfolder_path)) {
      dir.create(subfolder_path)
    }
    
    for (sample_size in sample_sizes) {
      for (dataset_index in 1:num_datasets_per_combination) {
        
        # 生成模拟数据
        result <- generate_simulation_data(
          n_sample = sample_size,
          template = template,
          n_feature = 20,
          metadata_effect_size = effect_size,
          perc_feature_spiked_metadata = 0.2,
          median_read_depth = 10000,
          # median_read_depth = 100000,
          alpha_0 = 1,
          alpha_T = 0.5,
          alpha_M_value = alpha_M_value,
          noise_sd = 0.5
        )
        
        # 保存生成的模拟数据到文件，以文件名标识样本大小和数据集编号
        file_name <- paste0("sample_", sample_size, "_dataset_", dataset_index, ".rds")
        saveRDS(result, file = file.path(subfolder_path, file_name))
      }
    }
  }
}

print("模拟数据生成完毕，并已保存在各自的文件夹中。")

############### median


output_dir <- "10-10-new_method_simulation_data_depth_new"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (median_read_depth in median_read_depths) {
  for (i in 1:50) {
    # Parameters defined within the loop
    n_sample <- 100  # Keeping the sample size fixed as 100
    n_feature <- 20
    metadata_effect_size <- 8
    perc_feature_spiked_metadata <- 0.2
    noise_sd <- 0.5
    
    result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = n_feature,
                                       metadata_effect_size = metadata_effect_size, perc_feature_spiked_metadata = perc_feature_spiked_metadata,
                                       median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T,
                                       alpha_M_value = alpha_M_value, noise_sd = noise_sd)
    save_result(result, "median_read_depth", median_read_depth, i)
  }
}
