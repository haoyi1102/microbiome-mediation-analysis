library(ggplot2)
library(reshape2)
library(dplyr)

# Function to generate power results from a CSV file
generate_power_results <- function(file_name, zero_percentage) {
  results_df <- read.csv(file_name)
  #"11-21_sample_test_1.csv"
  # Add sample size information
  results_df$sample_size <- rep(c( 100), each = 100)
  
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
  power_results$Zero_Percentage <- factor(zero_percentage, levels = c("low"))
  
  
  return(power_results)
}

# Generate power results for each sample
power_results_30 <- generate_power_results("11-21_sample_test_1.csv", "low")


# Combine all power results
all_power_results <- bind_rows( power_results_30)

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
















library(ggplot2)
library(reshape2)
library(dplyr)
library(knitr)

# Function to generate power results from a CSV file
generate_power_results <- function(file_name, zero_percentage) {
  results_df <- read.csv(file_name)
  
  # Add sample size information
  results_df$sample_size <- rep(c(100), each = 100)
  
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
  power_results$Zero_Percentage <- factor(zero_percentage, levels = c("low"))
  
  return(power_results)
}

# Generate power results for each sample
power_results_30 <- generate_power_results("11-21_sample_test_1.csv", "low")

# Combine all power results
all_power_results <- bind_rows(power_results_30)

# Create a table summarizing the results
power_table <- all_power_results %>%
  arrange(Zero_Percentage, Method) %>%
  select(Zero_Percentage, sample_size, Method, power)

# Display the table
kable(power_table, col.names = c("Zero Percentage", "Sample Size", "Method", "Power"), 
      caption = "Power Analysis Results by Zero Value Percentage and Method", format = "markdown")
