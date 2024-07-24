# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid) 

# Function to process data and calculate power
process_data <- function(file_path, sample_sizes, effect_size, n) {
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
  power_results$power[power_results$power == 0] <- 0.005
  
  return(power_results)
}

# Function to process data and calculate type-1 error
process_data_error <- function(file_path, sample_sizes, n) {
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
power_results_20 <- process_data("mediation_analysis_results_sample_new.csv", c(50, 100, 300), 20, 10)
power_results_1 <- process_data("mediation_analysis_results_sample.csv", c(50, 100, 300), 1, 20)
type1_error_results_20 <- process_data_error("mediation_analysis_results_new_no_mediation_effect_sample.csv", c(50, 100, 300), 10)
type1_error_results_1 <- process_data_error("mediation_analysis_results_no_mediation_effect__sample.csv", c(50, 100, 300), 20)

# Plot power and type-1 error results
simplify_plot <- function(plot, y_label) {
  plot +
    theme_minimal() +
    labs(x = "Method", y = y_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.box.margin = margin(t = -10)) +
    scale_fill_discrete(name = "Sample Size")
}

power_plot_1 <- simplify_plot(
  ggplot(power_results_1, aes(x = Method, y = power, fill = as.factor(sample_size))) +
    geom_bar(stat = "identity", position = "dodge"),
  "Power"
)

power_plot_20 <- simplify_plot(
  ggplot(power_results_20, aes(x = Method, y = power, fill = as.factor(sample_size))) +
    geom_bar(stat = "identity", position = "dodge"),
  "Power"
)

type1_error_plot_1 <- simplify_plot(
  ggplot(type1_error_results_1, aes(x = Method, y = type1_error, fill = as.factor(sample_size))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red"),
  "Type-1 Error"
)

type1_error_plot_20 <- simplify_plot(
  ggplot(type1_error_results_20, aes(x = Method, y = type1_error, fill = as.factor(sample_size))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red"),
  "Type-1 Error"
)

# Extract legend
extract_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

legend <- extract_legend(power_plot_1)

# Remove legend from individual plots
power_plot_1 <- power_plot_1 + theme(legend.position = "none")
power_plot_20 <- power_plot_20 + theme(legend.position = "none")
type1_error_plot_1 <- type1_error_plot_1 + theme(legend.position = "none")
type1_error_plot_20 <- type1_error_plot_20 + theme(legend.position = "none")

# Create labels for the plots
label_grobs <- list(
  textGrob("(a)", x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center", gp = gpar(fontsize = 10)),
  textGrob("(b)", x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center", gp = gpar(fontsize = 10)),
  textGrob("(c)", x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center", gp = gpar(fontsize = 10)),
  textGrob("(d)", x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center", gp = gpar(fontsize = 10))
)

# Combine plots with labels into a single layout
combined_plot <- arrangeGrob(
  grobs = list(
    arrangeGrob(power_plot_1, bottom = label_grobs[[1]]),
    arrangeGrob(power_plot_20, bottom = label_grobs[[2]]),
    arrangeGrob(type1_error_plot_1, bottom = label_grobs[[3]]),
    arrangeGrob(type1_error_plot_20, bottom = label_grobs[[4]])
  ),
  ncol = 2
)

final_plot <- grid.arrange(combined_plot, legend, ncol = 1, heights = c(10, 1))

# Display the final plot
print(final_plot)

