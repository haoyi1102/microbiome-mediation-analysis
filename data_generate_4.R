rm(list = ls())
library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)

# Load the generate_simulation_data function from data_generate_1.R
source("data_generate_1.R")

# Create a new directory to save the results
output_dir <- "simulation_data_4"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Function to save individual results to .rds file
save_result <- function(result, index) {
  rds_file <- file.path(output_dir, paste0("result_", index, ".rds"))
  saveRDS(result, file = rds_file)
}

# Fixed parameters for the simulation
n_sample <- 100

n_feature <- 20
perc_feature_spiked_metadata <- 0.2
median_read_depth <- 100000
alpha_0 <- 1
alpha_T <- 0.5
alpha_M_value <- 1
noise_sd <- 0.5

# Run the simulation 110 times, changing metadata_effect_size every 10 iterations
index <- 1
effect_sizes <- c(1, seq(5, 20, by = 5))
template <-  readRDS("sparsedossa_baseline.rds")
for (effect_size in effect_sizes) {
  for (i in 1:100) {
    result <- generate_simulation_data(n_sample = n_sample, template = template, n_feature = n_feature, 
                                       metadata_effect_size = effect_size, perc_feature_spiked_metadata = perc_feature_spiked_metadata, 
                                       median_read_depth = median_read_depth, alpha_0 = alpha_0, alpha_T = alpha_T, 
                                       alpha_M_value = alpha_M_value, noise_sd = noise_sd)
    save_result(result, index)
    index <- index + 1
  }
}

print("All results have been saved successfully.")
