rm(list = ls())
setwd("/home/haoyi/microbiome mediation analysis") 

# Load required libraries
library(MASS)
library(stats)
library(graphics)
library(energy)  
library(matrixStats)
library(devtools)
library(foreach)
library(parallel)
library(doParallel)
library(gtools)

# Load custom functions and dependencies
source("modima.R")
source("MedOmniTest.R")
source("highmed2019.r")
source("fromSKAT.R")

# Load additional mediation analysis packages
library(SparseMCMM)
library(ccmm)
library(LDM)
library(MicroBVS)
library(MedZIM)

# Set seed for reproducibility
set.seed(1234)

output_dir <- "11_18_simulation_data_per_feature_valid_new"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List all .rds files in the directory
rds_files <- list.files(output_dir, pattern = "*.rds", full.names = TRUE)
rds_files <- mixedsort(rds_files)

# Function to preprocess a single result
preprocess_single_result <- function(result) {
  Y_vector <- as.numeric(result$Y)
  T_vector <- as.numeric(result$T)
  M_a_matrix <- as.matrix(t(result$absolute_M))
  colnames(M_a_matrix) <- paste0("taxon_", 1:ncol(M_a_matrix))
  
  # Extract M_matrix and M_nz_matrix using a separate function
  matrices <- compute_matrices(M_a_matrix)
  
  list(Y_vector = Y_vector, T_vector = T_vector, M_a_matrix = M_a_matrix, M_matrix = matrices$M_matrix, M_nz_matrix = matrices$M_nz_matrix, feature = result$feature)
}

# Function to compute M_matrix and M_nz_matrix
compute_matrices <- function(M_a_matrix) {
  M_matrix <- M_a_matrix / rowSums(M_a_matrix)
  
  temp <- M_a_matrix
  temp[temp == 0] <- 0.5
  M_nz_matrix <- temp / rowSums(temp)
  
  list(M_matrix = M_matrix, M_nz_matrix = M_nz_matrix)
}

# Function to preprocess all results in a category
preprocess_results <- function(results) {
  lapply(results, preprocess_single_result)
}

# Split the results into categories
split_results <- function(files) {
  results <- list()
  for (file in files) {
    result <- readRDS(file)
    # Assuming the file naming convention includes the variable name and value
    name_parts <- strsplit(basename(file), "_")[[1]]
    
    if (length(name_parts) < 4) {
      next
    }
    
    variable <- name_parts[2]
    value <- name_parts[3]
    
    if (!is.null(results[[variable]])) {
      if (!is.null(results[[variable]][[value]])) {
        results[[variable]][[value]] <- c(results[[variable]][[value]], list(result))
      } else {
        results[[variable]][[value]] <- list(result)
      }
    } else {
      results[[variable]] <- list()
      results[[variable]][[value]] <- list(result)
    }
  }
  return(results)
}

# Function to set zero percentage in M_matrix
set_zero_percentage <- function(M_matrix, percentage) {
  num_samples <- nrow(M_matrix)
  num_features <- ncol(M_matrix)
  
  for (i in 1:num_samples) {
    zero_count <- sum(M_matrix[i, ] == 0)
    total_elements <- num_features
    current_zero_percentage <- (zero_count / total_elements) * 100
    
    while (current_zero_percentage < percentage) {
      non_zero_values <- M_matrix[i, M_matrix[i, ] != 0]
      if (length(non_zero_values) == 0) break  # Avoid infinite loop
      min_value_index <- which(M_matrix[i, ] == min(non_zero_values))
      M_matrix[i, min_value_index] <- 0
      zero_count <- sum(M_matrix[i, ] == 0)
      current_zero_percentage <- (zero_count / total_elements) * 100
    }
  }
  
  return(M_matrix)
}





# Split the results into categories
split_results <- split_results(rds_files)

# Preprocess data for each category
preprocessed_results_by_sample <- lapply(split_results$spiked, preprocess_results)
preprocessed_data <- preprocessed_results_by_sample[["features"]]

#i=2 i=3

# for (i in 1:length(preprocessed_data)){
#   data_list <- preprocessed_data[i]
#   data = data_list[[1]]
#   # data$feature
#   # Extract necessary vectors and matrices
#   Y_vector <- as.vector(data$Y_vector)
#   T_vector <- as.vector(data$T_vector)
#   M_matrix <- as.matrix(data$M_matrix)
#   M_a_matrix <- as.matrix(data$M_a_matrix)
#   M_nz_matrix <- as.matrix(data$M_nz_matrix)
#   zero_count <- sum(M_matrix == 0)
#   total_elements <- length(M_matrix)
#   zero_proportion <- (zero_count / total_elements) * 100
#   print(i)
#   print(zero_proportion)
#   print(data$feature)
# 
# }

# Zero percentage settings
# zero_percentages <- c( 30)
# 
# for (percentage in zero_percentages) { #percentage =30
results <- list()

for (i in 1:length(preprocessed_data)) { # i=850
  loop_start_time <- Sys.time()  # Start timing for the entire loop i=1
  
  data_list <- preprocessed_data[i]
  data = data_list[[1]]
  
  # Extract necessary vectors and matrices
  Y_vector <- as.vector(data$Y_vector)
  T_vector <- as.vector(data$T_vector)
  M_matrix <- as.matrix(data$M_matrix)
  M_a_matrix <- as.matrix(data$M_a_matrix)
  sum(M_matrix==0)/length(M_matrix)
  # Set zero percentage
  # M_matrix <- set_zero_percentage(M_matrix, percentage)
  #M_matrix= set_zero_percentage_new(count_matrix=M_matrix, target_depth = NULL, target_zero_proportion= percentage) 
  
  
  # Compute M_nz_matrix after setting zero percentage
  temp <- M_matrix
  temp[temp == 0] <- 0.5
  M_nz_matrix <- temp / rowSums(temp)
  
  # MODIMA
  modima_time <- system.time({
    modima_p_value <- tryCatch({
      T_dist <- dist(T_vector)
      Y_dist <- dist(Y_vector)
      M_dist <- dist(M_matrix)
      modima_result <- modima(T_dist, M_dist, Y_dist, nrep = 999)
      modima_result$p.value
    }, error = function(e) {
      message(paste("Error in MODIMA:", e$message))
      NA  # Return NA if an error occurs
    })
  })
  print(paste("MODIMA completed in:", modima_time[3], "seconds"))
  
  # MEDTEST
  medtest_time <- system.time({
    medtest_p_value <- tryCatch({
      m_list <- list(euc = dist(M_matrix))
      medtest_result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)
      medtest_result$permP
    }, error = function(e) {
      message(paste("Error in MEDTEST:", e$message))
      NA  # Return NA if an error occurs
    })
  })
  print(paste("MEDTEST completed in:", medtest_time[3], "seconds"))
  
  
  # CCMM
  ccmm_time <- system.time({
    
    tryCatch({
      ccmm_result <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
                          sig.level = 0.05, tol = 1e-06, max.iter = 5000)
      ccmm_CI_lower <- as.numeric(ccmm_result$TIDE.CI["2.5%"])
      ccmm_CI_upper <- as.numeric(ccmm_result$TIDE.CI["97.5%"])
    }, error = function(e) {
      
      ccmm_CI_lower <- -1
      ccmm_CI_upper <- 1
      
      print(paste("Error in CCMM:", e$message))
    })
  })
  print(paste("CCMM completed in:", ccmm_time[3], "seconds")) 
  
  # LDM-MED
  ldm_time <- system.time({
    data_ldm <- data.frame(Y = Y_vector, T = T_vector, M_nz_matrix = M_nz_matrix)
    ldm_p_value <- tryCatch({
      ldm_result <- ldm(
        formula = M_nz_matrix ~ T + Y,
        data = data_ldm,
        seed = 1234,
        test.mediation = TRUE
      )
      ldm_result$med.p.global.omni
    }, error = function(e) {
      NA  
    })
  })
  print(paste("LDM-MED completed in:", ldm_time[3], "seconds"))
  
  # PERMANOVA-MED
  permanova_time <- system.time({
    permanova_p_value <- tryCatch({
      data_ldm <- data.frame(Y = Y_vector, T = T_vector)
      formula <- as.formula("M_nz_matrix ~ T + Y")
      res_perm <- permanovaFL(
        formula = formula,
        data = data_ldm,
        dist.method = "bray",  
        seed = 1234,
        n.perm.max = 1000,
        n.cores = 4,
        test.mediation = TRUE,
        verbose = TRUE
      )
      res_perm$med.p.permanova
    }, error = function(e) {
      NA  
    })
  })
  print(paste("PERMANOVA-MED completed in:", permanova_time[3], "seconds"))
  
  loop_end_time <- Sys.time()  # End timing for the entire loop
  print(paste("Loop", i, "completed in:", loop_end_time - loop_start_time, "mins"))
  
  # Store results
  results[[i]] <- list(
    modima_p_value = modima_p_value,
    medtest_p_value = medtest_p_value,
    ccmm_CI_lower = ccmm_CI_lower,
    ccmm_CI_upper = ccmm_CI_upper,
    ldm_p_value = ldm_p_value,
    permanova_p_value = permanova_p_value
  )
}

# Convert results to a data frame
results_df <- do.call(rbind, lapply(results, function(x) {
  c(modima_p_value = x$modima_p_value,
    medtest_p_value = x$medtest_p_value,
    ccmm_CI_lower = x$ccmm_CI_lower,
    ccmm_CI_upper = x$ccmm_CI_upper,
    ldm_p_value = x$ldm_p_value,
    permanova_p_value = x$permanova_p_value)
}))

# Convert the list of results into a data frame
results_df <- as.data.frame(results_df)

# Save the data frame to a CSV file
write.csv(results_df, file = paste0("per_valid_taxa", ".csv"), row.names = FALSE)

print(paste("Results for per", "% zero values have been saved"))


print("All results have been saved")

















# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Function to generate power results from a CSV file for `per_valid_taxa`
generate_taxa_power_results <- function(file_name) {
  results_df <- read.csv(file_name)
  
  # 保留前 800 行数据以对齐更新后的 num_spiked_features
  results_df <- head(results_df, 800)
  
  # Update feature spiked proportions (removed 1)
  num_spiked_features <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  results_df$feature_spiked <- rep(num_spiked_features, each = 80)
  
  # Replace missing values with 1 (acceptance)
  results_df[is.na(results_df)] <- 1
  
  # Rename methods for better readability
  results_long <- melt(results_df, id.vars = "feature_spiked", 
                       measure.vars = c("modima_p_value", "medtest_p_value", "ldm_p_value", "permanova_p_value"),
                       variable.name = "Method", value.name = "p_value")
  
  results_long$Method <- recode(results_long$Method,
                                "modima_p_value" = "MODIMA",
                                "medtest_p_value" = "MedTest",
                                "ldm_p_value" = "LDM",
                                "permanova_p_value" = "PERMANOVA")
  
  # Calculate power (p-value < 0.05)
  power_results <- results_long %>%
    group_by(feature_spiked, Method) %>%
    summarise(power = mean(p_value < 0.05, na.rm = TRUE))
  
  # Calculate rejection rate based on confidence intervals (reject if not containing 0)
  results_df <- results_df %>%
    mutate(ccmm_reject = ifelse(ccmm_CI_lower > 0 | ccmm_CI_upper < 0, 1, 0))
  
  # Add ccmm_reject to power_results
  ccmm_power_results <- results_df %>%
    group_by(feature_spiked) %>%
    summarise(Method = "CCMM", power = mean(ccmm_reject, na.rm = TRUE))
  
  # Combine all methods' power_results
  power_results <- bind_rows(power_results, ccmm_power_results)
  
  # Replace zero power values with a small value
  power_results$power[power_results$power == 0] <- 0.001
  
  return(power_results)
}

# Generate power results for `per_valid_taxa`
power_results_taxa <- generate_taxa_power_results("per_valid_taxa.csv")

# Combine all power results (if multiple datasets are used, append them here)
all_power_results_taxa <- power_results_taxa

# Plot power analysis for `per_valid_taxa` with a line chart
ggplot(all_power_results_taxa, aes(x = feature_spiked, y = power, color = Method, group = Method)) +
  geom_line(size = 0.8) +  # Make lines thinner
  geom_point(size = 1.5) +  # Smaller points
  theme_minimal() +
  labs(title = "Power Analysis by Feature Spiked Proportion",
       x = "Proportion of Features Spiked",
       y = "Power") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",  # Move legend to bottom
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 10),  # Adjust legend text size
        plot.title = element_text(hjust = 0.5)) +  # Center the plot title
  scale_color_brewer(palette = "Set1") +  # Use a color palette for better distinction
  coord_cartesian(ylim = c(0, 1))  # Fixed y-axis range
