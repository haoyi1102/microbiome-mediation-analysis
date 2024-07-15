


# data("Stool_subset")
# Stool_subset[3,2]<-1678541
# Stool_subset[4,2]<-2056565
# Stool_subset[4,4]<-125656
# IBD <- fit_SparseDOSSA2(lambda = 1,data = Stool_subset)
# ?fit_SparseDOSSA2
# ?control_fit


generate_simulation_data <- function(n_sample = 100, template = "Stool", n_feature = 20, 
                                     metadata_effect_size = 1, perc_feature_spiked_metadata = 0.1, 
                                     median_read_depth = 10000, alpha_0 = 1, alpha_T = 0.5, 
                                     alpha_M_value = 0.2, noise_sd = 0.1) {


temp = n_sample/2
metadata <- data.frame(treatment = rep(c("Control", "Treatment"), each = temp))
metadata$treatment_binary <- ifelse(metadata$treatment == "Control", 0, 1)
metadata_matrix <- as.matrix(metadata$treatment_binary)

# Define a function to generate data; some seeds may not work
generate_data <- function(seed) {
  set.seed(seed)
  tryCatch({
    sim_data <- SparseDOSSA2(
      template = template, 
      n_sample = n_sample, 
      n_feature = n_feature,
      spike_metadata = "abundance", 
      metadata_effect_size = metadata_effect_size,
      perc_feature_spiked_metadata = perc_feature_spiked_metadata,
      metadata_matrix = metadata_matrix,
      median_read_depth = median_read_depth,  
      verbose = TRUE
    )
    return(sim_data)
  }, error = function(e) {
    return(NULL)
  })
}

# Run in a loop until successful
sim_data <- NULL
while(is.null(sim_data)) {
  seed <- sample(1:10000, 1)
  cat("Using seed:", seed, "\n")
  sim_data <- generate_data(seed = seed)
}

# Extract generated count data
synthetic_data <- sim_data$simulated_data
generated_metadata <- sim_data$spike_metadata$metadata_matrix
feature_metadata_spike_df <- sim_data$spike_metadata$feature_metadata_spike_df

zero_count <- sum(synthetic_data == 0)
total_elements <- length(synthetic_data)
zero_proportion <- zero_count / total_elements


#
# compositional can lead to co-linearity problems that make traditional variance calculations inapplicable.
# # Function to convert count data to compositional data
# to_compositional <- function(data) {
#   t(apply(data, 2, function(x) x / sum(x)))
# }
# 
# # Convert synthetic_data to compositional data
# compositional_data <- to_compositional(synthetic_data)
# 
# # Calculate variance for each sample
# sample_variances <- apply(compositional_data, 2, var)
# print("Sample Variances:")
# print(sample_variances)

# Aitchison #Function to convert count data to compositional data with pseudo-count
to_compositional <- function(data, pseudo_count = 1e-6) {
  data <- data + pseudo_count
  t(apply(data, 2, function(x) x / sum(x)))
}

# Convert synthetic_data to compositional data
compositional_data <- to_compositional(synthetic_data)

# Calculate Aitchison variance for each sample
aitchison_variance <- function(comp_data) {
  log_data <- log(comp_data)
  clr_data <- t(apply(log_data, 1, function(x) x - mean(x)))
  return(apply(clr_data, 2, var))
}

# Apply the function to compositional data
sample_variances <- aitchison_variance(compositional_data)

mean_variance <- mean(sample_variances)



# Convert treatment group variable to numeric: Control -> 0, Treatment -> 1
metadata$treatment_numeric <- as.numeric(metadata$treatment == "Treatment")

# Initialize alpha_M with the same length as n_feature - 1
alpha_M <- rep(alpha_M_value, n_feature - 1)

# **Adjust alpha_M coefficients based on feature_metadata_spike_df**
for (i in 1:nrow(feature_metadata_spike_df)) {
  feature_index <- as.numeric(gsub("Feature", "", feature_metadata_spike_df$feature_spiked[i]))
  if (feature_metadata_spike_df$effect_size[i] < 0) {
    alpha_M[feature_index] <- -alpha_M_value
  }
}



# Initialize outcome Y
outcome_Y <- rep(alpha_0, n_sample)

# Add treatment variable term
outcome_Y <- outcome_Y + metadata$treatment_numeric * alpha_T

# Compute log contrasts and add to outcome Y
for (i in 1:n_sample) {
  log_contrasts <- numeric(n_feature - 1)
  for (j in 1:(n_feature - 1)) {
    ratio <- (synthetic_data[j, i] + 1) / (synthetic_data[n_feature, i] + 1)
    if (ratio == 0) {
      log_contrasts[j] <- 0
    } else {
      log_contrasts[j] <- log(ratio)
    }
  }
  outcome_Y[i] <- outcome_Y[i] + sum(alpha_M * log_contrasts)
}

# Add noise term
outcome_Y <- outcome_Y + rnorm(n_sample, mean = 0, sd = noise_sd)
# Return the results
list(T = metadata$treatment_numeric, Y = outcome_Y, absolute_M = synthetic_data, feature = feature_metadata_spike_df,
     zero_proportion = zero_proportion, mean_variance = mean_variance)
}


