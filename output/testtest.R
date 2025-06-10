
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
  type1_error_results$Zero_Percentage <- zero_percentage
  
  return(type1_error_results)
}

# Generate Type-1 error results for each sample
type1_error_results_10 <- generate_type1_error_results("sample_noMe_10.csv", "10% Zero Values")
type1_error_results_30 <- generate_type1_error_results("sample_noMe_30.csv", "30% Zero Values")
type1_error_results_50 <- generate_type1_error_results("sample_noMe_50.csv", "50% Zero Values")
type1_error_results_70 <- generate_type1_error_results("sample_noMe_70.csv", "70% Zero Values")

# Combine all Type-1 error results
all_type1_error_results <- bind_rows(type1_error_results_10, type1_error_results_30, type1_error_results_50, type1_error_results_70)

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

