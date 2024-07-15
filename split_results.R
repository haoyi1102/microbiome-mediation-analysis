# Function to split the results into different categories
split_results <- function(files) {
  results_by_feature <- list()
  results_by_sample <- list()
  results_by_metadata_effect_size <- list()
  results_by_perc_feature_spiked_metadata <- list()
  results_by_median_read_depth <- list()
  results_by_noise_sd <- list()
  
  for (file in files) {
    result <- readRDS(file)
    if (grepl("n_feature_", file)) {
      n_feature <- sub(".*n_feature_(\\d+).*", "\\1", file)
      results_by_feature[[n_feature]] <- result
    } else if (grepl("n_sample_", file)) {
      n_sample <- sub(".*n_sample_(\\d+).*", "\\1", file)
      results_by_sample[[n_sample]] <- result
    } else if (grepl("metadata_effect_size_", file)) {
      metadata_effect_size <- sub(".*metadata_effect_size_(\\d+\\.?\\d*).*", "\\1", file)
      results_by_metadata_effect_size[[metadata_effect_size]] <- result
    } else if (grepl("perc_feature_spiked_metadata_", file)) {
      perc_feature_spiked_metadata <- sub(".*perc_feature_spiked_metadata_(\\d+\\.?\\d*).*", "\\1", file)
      results_by_perc_feature_spiked_metadata[[perc_feature_spiked_metadata]] <- result
    } else if (grepl("median_read_depth_", file)) {
      median_read_depth <- sub(".*median_read_depth_(\\d+).*", "\\1", file)
      results_by_median_read_depth[[median_read_depth]] <- result
    } else if (grepl("noise_sd_", file)) {
      noise_sd <- sub(".*noise_sd_(\\d+\\.?\\d*).*", "\\1", file)
      results_by_noise_sd[[noise_sd]] <- result
    }
  }
  
  list(
    results_by_feature = results_by_feature,
    results_by_sample = results_by_sample,
    results_by_metadata_effect_size = results_by_metadata_effect_size,
    results_by_perc_feature_spiked_metadata = results_by_perc_feature_spiked_metadata,
    results_by_median_read_depth = results_by_median_read_depth,
    results_by_noise_sd = results_by_noise_sd
  )
}