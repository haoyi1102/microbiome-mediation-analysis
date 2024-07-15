rm(list = ls())
library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)

# Load the generate_simulation_data function from data_generate_1.R
source("data_generate_1.R")


# # Use the generate_simulation_data function with desired parameters
# result <- generate_simulation_data(n_sample = 500, template = "Stool", n_feature = 20,
#                                    metadata_effect_size = 1, perc_feature_spiked_metadata = 0.1,
#                                    median_read_depth = 10000, alpha_0 = 1, alpha_T = 0.5,
#                                    alpha_M_value = 0.2, noise_sd = 0.5)
# # 
# # # Print the results
# # print(result)

# Create a new directory to save the results
output_dir <- "simulation_data_1"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Function to save individual results to .rds file
save_result <- function(result, variable, value) {
  rds_file <- file.path(output_dir, paste0("result_", variable, "_", value, ".rds"))
  saveRDS(result, file = rds_file)
}

# Fixed parameters
template <- "Stool"
alpha_0 <- 1
alpha_T <- 0.5
alpha_M_value <- 0.2

# Parameter values to iterate over
n_features <- c(10, 20, 50)
n_samples <- c(50, 100, 300)
metadata_effect_sizes <- c(0.5, 1, 10)
perc_feature_spiked_metadatas <- c(0.1, 0.2, 0.5)
median_read_depths <- c(100, 1000, 10000)
noise_sds <- c(0.5, 1, 3)

# Generate datasets and save the results, changing one parameter at a time
for (n_feature in n_features) {
  result <- generate_simulation_data(n_sample = 100, template = template, n_feature = n_feature, 
                                     metadata_effect_size = 1, perc_feature_spiked_metadata = 0.1, 
                                     median_read_depth = 10000, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                     alpha_M_value = alpha_M_value, noise_sd = 0.5)
  save_result(result, "n_feature", n_feature)
}

for (n_sample in n_samples) {
  result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = 20, 
                                     metadata_effect_size = 1, perc_feature_spiked_metadata = 0.1, 
                                     median_read_depth = 10000, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                     alpha_M_value = alpha_M_value, noise_sd = 0.5)
  save_result(result, "n_sample", n_sample)
}

for (metadata_effect_size in metadata_effect_sizes) {
  result <- generate_simulation_data(n_sample = 100, template = template, n_feature = 20, 
                                     metadata_effect_size = metadata_effect_size, perc_feature_spiked_metadata = 0.1, 
                                     median_read_depth = 10000, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                     alpha_M_value = alpha_M_value, noise_sd = 0.5)
  save_result(result, "metadata_effect_size", metadata_effect_size)
}

for (perc_feature_spiked_metadata in perc_feature_spiked_metadatas) {
  result <- generate_simulation_data(n_sample = 100, template = template, n_feature = 20, 
                                     metadata_effect_size = 1, perc_feature_spiked_metadata = perc_feature_spiked_metadata, 
                                     median_read_depth = 10000, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                     alpha_M_value = alpha_M_value, noise_sd = 0.5)
  save_result(result, "perc_feature_spiked_metadata", perc_feature_spiked_metadata)
}

for (median_read_depth in median_read_depths) {
  result <- generate_simulation_data(n_sample = 100, template = template, n_feature = 20, 
                                     metadata_effect_size = 1, perc_feature_spiked_metadata = 0.1, 
                                     median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                     alpha_M_value = alpha_M_value, noise_sd = 0.5)
  save_result(result, "median_read_depth", median_read_depth)
}

for (noise_sd in noise_sds) {
  result <- generate_simulation_data(n_sample = 100, template = template, n_feature = 20, 
                                     metadata_effect_size = 1, perc_feature_spiked_metadata = 0.1, 
                                     median_read_depth = 10000, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                     alpha_M_value = alpha_M_value, noise_sd = noise_sd)
  save_result(result, "noise_sd", noise_sd)
}

print("All results have been saved successfully.")