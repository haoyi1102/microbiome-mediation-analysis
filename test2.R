# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

# Function to process data and calculate power
process_data <- function(file_path, sample_sizes, effect_size,n) {
  results_df <- read.csv(file_path)
  results_df$sample_size <- rep(sample_sizes, each = n)
  results_df[is.na(results_df)] <- 1
  
  results_long <- melt(results_df, id.vars = "sample_size", 
                       measure.vars = c("modima_p_value", "medtest_p_value", 
                                        "sparsemcmm_ome_p_value.OME", "ldm_p_value", 
                                        "permanova_p_value"),
                       variable.name = "Method", value.name = "p_value")
  
  results_long$Method <- recode(results_long$Method,
                                "modima_p_value" = "MODIMA",
                                "medtest_p_value" = "MedTest",
                                "sparsemcmm_ome_p_value.OME" = "SparseMCMM",
                                "ldm_p_value" = "LDM",
                                "permanova_p_value" = "PERMANOVA")
  
  power_results <- results_long %>%
    group_by(sample_size, Method) %>%
    summarise(power = mean(p_value < 0.05, na.rm = TRUE))
  
  results_df <- results_df %>%
    mutate(ccmm_reject = ifelse(ccmm_CI_lower > 0 | ccmm_CI_upper < 0, 1, 0),
           microbvs_reject = ifelse(`microbvs_result.2.5.` > 0 | `microbvs_result.97.5.` < 0, 1, 0))
  
  ccmm_power_results <- results_df %>%
    group_by(sample_size) %>%
    summarise(Method = "CCMM", power = mean(ccmm_reject, na.rm = TRUE))
  
  microbvs_power_results <- results_df %>%
    group_by(sample_size) %>%
    summarise(Method = "MicroBVS", power = mean(microbvs_reject, na.rm = TRUE))
  
  power_results <- bind_rows(power_results, ccmm_power_results, microbvs_power_results)
  power_results$power[power_results$power == 0] <- 0.001
  
  return(power_results)
}

# Function to process data and calculate type-1 error
process_data_error <- function(file_path, sample_sizes,n) {
  results_df <- read.csv(file_path)
  results_df$sample_size <- rep(sample_sizes, each = n)
  results_df[is.na(results_df)] <- 1
  
  results_long <- melt(results_df, id.vars = "sample_size", 
                       measure.vars = c("modima_p_value", "medtest_p_value", 
                                        "sparsemcmm_ome_p_value.OME", "ldm_p_value", 
                                        "permanova_p_value"),
                       variable.name = "Method", value.name = "p_value")
  
  results_long$Method <- recode(results_long$Method,
                                "modima_p_value" = "MODIMA",
                                "medtest_p_value" = "MedTest",
                                "sparsemcmm_ome_p_value.OME" = "SparseMCMM",
                                "ldm_p_value" = "LDM",
                                "permanova_p_value" = "PERMANOVA")
  
  type1_error_results <- results_long %>%
    group_by(sample_size, Method) %>%
    summarise(type1_error = mean(p_value < 0.05, na.rm = TRUE))
  
  results_df <- results_df %>%
    mutate(ccmm_reject = ifelse(ccmm_CI_lower > 0 | ccmm_CI_upper < 0, 1, 0),
           microbvs_reject = ifelse(`microbvs_result.2.5.` > 0 | `microbvs_result.97.5.` < 0, 1, 0))
  
  ccmm_type1_error_results <- results_df %>%
    group_by(sample_size) %>%
    summarise(Method = "CCMM", type1_error = mean(ccmm_reject, na.rm = TRUE))
  
  microbvs_type1_error_results <- results_df %>%
    group_by(sample_size) %>%
    summarise(Method = "MicroBVS", type1_error = mean(microbvs_reject, na.rm = TRUE))
  
  type1_error_results <- bind_rows(type1_error_results, ccmm_type1_error_results, microbvs_type1_error_results)
  type1_error_results$type1_error[type1_error_results$type1_error == 0] <- 0.001
  
  return(type1_error_results)
}

# Process data
power_results_20 <- process_data("mediation_analysis_results_sample_new.csv", c(50, 100, 300), 20,10)
power_results_1 <- process_data("mediation_analysis_results_sample.csv", c(50, 100, 300), 1,20)
type1_error_results_20 <- process_data_error("mediation_analysis_results_new_no_mediation_effect_sample.csv", c(50, 100, 300),10)
type1_error_results_1 <- process_data_error("mediation_analysis_results_no_mediation_effect__sample.csv", c(50, 100, 300),20)

# Plot power and type-1 error results
power_plot_20 <- ggplot(power_results_20, aes(x = Method, y = power, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Power Analysis (Effect Size = 20)",
       x = "Method",
       y = "Power") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample Size")

power_plot_1 <- ggplot(power_results_1, aes(x = Method, y = power, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Power Analysis (Effect Size = 1)",
       x = "Method",
       y = "Power") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample Size")

type1_error_plot_20 <- ggplot(type1_error_results_20, aes(x = Method, y = type1_error, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Type-1 Error Analysis (Effect Size = 20)",
       x = "Method",
       y = "Type-1 Error") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample Size")

type1_error_plot_1 <- ggplot(type1_error_results_1, aes(x = Method, y = type1_error, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Type-1 Error Analysis (Effect Size = 1)",
       x = "Method",
       y = "Type-1 Error") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample Size")

# Combine all plots into a single layout
grid.arrange(power_plot_20, power_plot_1, type1_error_plot_20, type1_error_plot_1, ncol = 2)