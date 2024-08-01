rm(list = ls())

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
##
source("modima.R")
source("MedOmniTest.R")
# install_github("chanw0/SparseMCMM")
library(SparseMCMM)
# install.packages("ccmm")
library(ccmm)
# install_github("yijuanhu/LDM")
library(LDM)
# install_github( "mkoslovsky/MicroBVS")
library(MicroBVS)
# devtools::install_github("quranwu/MedZIM")
library(MedZIM)
# install.packages("rrBLUP")
# install.packages("CompQuadForm")
source("highmed2019.R")
source('fromSKAT.R')
set.seed(1234)

#mediation effect

# read the results
output_dir <- "simulation_data_sample_new"

# List all .rds files in the directory
rds_files <- list.files(output_dir, pattern = "*.rds", full.names = TRUE)
rds_files <- mixedsort(rds_files)
# Function to preprocess a single result
preprocess_single_result <- function(result) {
  Y_vector <- as.numeric(result$Y)
  T_vector <- as.numeric(result$T)
  M_a_matrix <- as.matrix(t(result$absolute_M))
  colnames(M_a_matrix) <- paste0("taxon_", 1:ncol(M_a_matrix))
  M_matrix <- M_a_matrix / rowSums(M_a_matrix)
  temp <- M_a_matrix
  temp[temp == 0] <- 0.5
  M_nz_matrix <- temp / rowSums(temp)
  
  list(Y_vector = Y_vector, T_vector = T_vector, M_a_matrix = M_a_matrix, M_matrix = M_matrix, M_nz_matrix = M_nz_matrix, feature = result$feature)
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

# Split the results into categories
split_results <- split_results(rds_files)

# Preprocess data for each category
preprocessed_results_by_sample <- lapply(split_results$n, preprocess_results)

# Mediation analysis for preprocessed data
results <- list()

preprocessed_data <- preprocessed_results_by_sample[["sample"]]
#i=3
for (i in 1:length(preprocessed_data)) {
  
  loop_start_time <- Sys.time()  # Start timing for the entire loop
  
  data_list <- preprocessed_data[i]
  data = data_list[[1]]
   # data$feature
  # Extract necessary vectors and matrices
  Y_vector <- as.vector(data$Y_vector)
  T_vector <- as.vector(data$T_vector)
  M_matrix <- as.matrix(data$M_matrix)
  M_a_matrix <- as.matrix(data$M_a_matrix)
  M_nz_matrix <- as.matrix(data$M_nz_matrix)
  
  # MODIMA
  modima_time <- system.time({
    T_dist <- dist(T_vector)
    Y_dist <- dist(Y_vector)
    M_dist <- dist(M_matrix)
    modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
    modima_p_value <- modima_result$p.value
  })
  print(paste("MODIMA completed in:", modima_time[3], "seconds"))
  
  # MEDTEST
  medtest_time <- system.time({
    m_list <- list(euc = dist(M_matrix))
    medtest_result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)
    medtest_p_value <- medtest_result$permP
  })
  print(paste("MEDTEST completed in:", medtest_time[3], "seconds"))
  
  # CCMM
  ccmm_time <- system.time({
    ccmm_result <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
                        sig.level = 0.05, tol = 1e-06, max.iter = 5000)
    ccmm_CI_lower <- as.numeric(ccmm_result$TIDE.CI["2.5%"])
    ccmm_CI_upper <- as.numeric(ccmm_result$TIDE.CI["97.5%"])
  })
  print(paste("CCMM completed in:", ccmm_time[3], "seconds"))
  
  # SparseMCMM
  sparsemcmm_time <- system.time({
    sparsemcmm_ome_p_value <- tryCatch({
      res <- SparseMCMM(T_vector, M_nz_matrix, Y_vector, n.split = 1, num.per = 50)
      ress <- res$Test
      ress["OME"]
    }, error = function(e) {
      NULL  
    })
  })
  print(paste("SparseMCMM completed in:", sparsemcmm_time[3], "seconds"))
  
  # LDM-MED
  ldm_time <- system.time({
    data_ldm <- data.frame(Y = Y_vector, T = T_vector,M_nz_matrix = M_nz_matrix)
    ldm_p_value <- tryCatch({
     ldm_result <- ldm(
        formula = M_nz_matrix ~ T + Y,
        data = data_ldm,
        seed = 1234,
        test.mediation = TRUE
      )
      ldm_result$med.p.global.omni
    }, error = function(e) {
      NULL  
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
      NULL  
    })
  })
  print(paste("PERMANOVA-MED completed in:", permanova_time[3], "seconds"))
  
  # microbvs
  microbvs_time <- system.time({
    model_real <- MCMC_Med(iterations = 3000, trt = T_vector, Y = Y_vector, Z = M_matrix, taxa = 2)
    result_real <- Selection_Med2(model = model_real)
    microbvs_result <- result_real$global_quantile
  })
  print(paste("microbvs completed in:", microbvs_time[3], "seconds"))
  
  # # Medzim
  # medzim_time <- system.time({
  #   taxon_name <- "taxon_"
  #   libsize <- rowSums(M_a_matrix)
  #   dat <- data.frame(Y_vector, T_vector, libsize)
  #   dat <- cbind(dat, M_matrix)
  #   MedZIM_results <- MedZIM_func(
  #     dat = dat,
  #     xVar = "T_vector",
  #     yVar = "Y_vector",
  #     taxon_name = taxon_name,
  #     libSize_name = "libsize",
  #     obs_gt_0 = 2,
  #     obs_eq_0 = 2,
  #     inter_x_mg0 = TRUE,
  #     inter_x_m = FALSE,
  #     eval.max = 200,
  #     iter.max = 200,
  #     x_from = 0,
  #     x_to = 1,
  #     type1error = 0.05,
  #     paraJobs = 2
  #   )
  
  #   first_taxon_result <- MedZIM_results$fullList
  #   first_taxon_result <- first_taxon_result[[1]]
  #   first_taxon_result <- first_taxon_result[[1]]
  #   MedZIM_pvalue <- first_taxon_result["pNIE"]
  # })
  # print(paste("Medzim completed in:", medzim_time[3], "seconds"))
  
  loop_end_time <- Sys.time()  # End timing for the entire loop
  print(paste("Loop", i, "completed in:", loop_end_time - loop_start_time, "mins"))
  
  # Store results
  results[[i]] <- list(
    modima_p_value = modima_p_value,
    medtest_p_value = medtest_p_value,
    ccmm_CI_lower = ccmm_CI_lower,
    ccmm_CI_upper = ccmm_CI_upper,
    sparsemcmm_ome_p_value = sparsemcmm_ome_p_value,
    ldm_p_value = ldm_p_value,
    permanova_p_value = permanova_p_value,
    microbvs_result = microbvs_result
   
  )
}

# Convert results to a data frame
results_df <- do.call(rbind, lapply(results, function(x) {
  c(modima_p_value = x$modima_p_value,
    medtest_p_value = x$medtest_p_value,
    ccmm_CI_lower = x$ccmm_CI_lower,
    ccmm_CI_upper = x$ccmm_CI_upper,
    sparsemcmm_ome_p_value = x$sparsemcmm_ome_p_value,
    ldm_p_value = x$ldm_p_value,
    permanova_p_value = x$permanova_p_value,
    microbvs_result = x$microbvs_result
    
    )
}))

# Convert the list of results into a data frame
results_df <- as.data.frame(results_df)

# Save the data frame to a CSV file
write.csv(results_df, file = "mediation_analysis_results_sample_new.csv", row.names = FALSE)

print("Results have been saved")

############## end mediation effect

##################  no mediation effect
# read the results
output_dir <- "simulation_data_sample_new_no_mediation"

# List all .rds files in the directory
rds_files <- list.files(output_dir, pattern = "*.rds", full.names = TRUE)
rds_files <- mixedsort(rds_files)
# Function to preprocess a single result
preprocess_single_result <- function(result) {
  Y_vector <- as.numeric(result$Y)
  T_vector <- as.numeric(result$T)
  M_a_matrix <- as.matrix(t(result$absolute_M))
  colnames(M_a_matrix) <- paste0("taxon_", 1:ncol(M_a_matrix))
  M_matrix <- M_a_matrix / rowSums(M_a_matrix)
  temp <- M_a_matrix
  temp[temp == 0] <- 0.5
  M_nz_matrix <- temp / rowSums(temp)
  
  list(Y_vector = Y_vector, T_vector = T_vector, M_a_matrix = M_a_matrix, M_matrix = M_matrix, M_nz_matrix = M_nz_matrix, feature = result$feature)
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

# Split the results into categories
split_results <- split_results(rds_files)

# Preprocess data for each category
preprocessed_results_by_sample <- lapply(split_results$n, preprocess_results)

# Mediation analysis for preprocessed data
results <- list()

preprocessed_data <- preprocessed_results_by_sample[["sample"]]
#i=1
for (i in 1:length(preprocessed_data)) {
  
  loop_start_time <- Sys.time()  # Start timing for the entire loop
  
  data_list <- preprocessed_data[i]
  data = data_list[[1]]
  
  # Extract necessary vectors and matrices
  Y_vector <- as.vector(data$Y_vector)
  T_vector <- as.vector(data$T_vector)
  M_matrix <- as.matrix(data$M_matrix)
  M_a_matrix <- as.matrix(data$M_a_matrix)
  M_nz_matrix <- as.matrix(data$M_nz_matrix)
  
  # MODIMA
  modima_time <- system.time({
    T_dist <- dist(T_vector)
    Y_dist <- dist(Y_vector)
    M_dist <- dist(M_matrix)
    modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
    modima_p_value <- modima_result$p.value
  })
  print(paste("MODIMA completed in:", modima_time[3], "seconds"))
  
  # MEDTEST
  medtest_time <- system.time({
    m_list <- list(euc = dist(M_matrix))
    medtest_result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)
    medtest_p_value <- medtest_result$permP
  })
  print(paste("MEDTEST completed in:", medtest_time[3], "seconds"))
  
  # CCMM
  ccmm_time <- system.time({
    ccmm_result <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
                        sig.level = 0.05, tol = 1e-06, max.iter = 5000)
    ccmm_CI_lower <- as.numeric(ccmm_result$TIDE.CI["2.5%"])
    ccmm_CI_upper <- as.numeric(ccmm_result$TIDE.CI["97.5%"])
  })
  print(paste("CCMM completed in:", ccmm_time[3], "seconds"))
  
  # SparseMCMM
  sparsemcmm_time <- system.time({
    sparsemcmm_ome_p_value <- tryCatch({
      res <- SparseMCMM(T_vector, M_nz_matrix, Y_vector, n.split = 1, num.per = 10)
      ress <- res$Test
      ress["OME"]
    }, error = function(e) {
      NULL  
    })
  })
  print(paste("SparseMCMM completed in:", sparsemcmm_time[3], "seconds"))
  
  # LDM-MED
  ldm_time <- system.time({
    data_ldm <- data.frame(Y = Y_vector, T = T_vector)
    ldm_p_value <- tryCatch({
      ldm_result <- ldm(
        formula = M_nz_matrix ~ T + Y,
        data = data_ldm,
        seed = 1234,
        test.mediation = TRUE
      )
      ldm_result$med.p.global.omni
    }, error = function(e) {
      NULL  
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
      NULL  
    })
  })
  print(paste("PERMANOVA-MED completed in:", permanova_time[3], "seconds"))
  
  # microbvs
  microbvs_time <- system.time({
    model_real <- MCMC_Med(iterations = 3000, trt = T_vector, Y = Y_vector, Z = M_matrix, taxa = 2)
    result_real <- Selection_Med2(model = model_real)
    microbvs_result <- result_real$global_quantile
  })
  print(paste("microbvs completed in:", microbvs_time[3], "seconds"))
  
  # # Medzim
  # medzim_time <- system.time({
  #   taxon_name <- "taxon_"
  #   libsize <- rowSums(M_a_matrix)
  #   dat <- data.frame(Y_vector, T_vector, libsize)
  #   dat <- cbind(dat, M_matrix)
  #   MedZIM_results <- MedZIM_func(
  #     dat = dat,
  #     xVar = "T_vector",
  #     yVar = "Y_vector",
  #     taxon_name = taxon_name,
  #     libSize_name = "libsize",
  #     obs_gt_0 = 2,
  #     obs_eq_0 = 2,
  #     inter_x_mg0 = TRUE,
  #     inter_x_m = FALSE,
  #     eval.max = 200,
  #     iter.max = 200,
  #     x_from = 0,
  #     x_to = 1,
  #     type1error = 0.05,
  #     paraJobs = 2
  #   )
  #   first_taxon_result <- MedZIM_results$fullList
  #   first_taxon_result <- first_taxon_result[[1]]
  #   first_taxon_result <- first_taxon_result[[1]]
  #   MedZIM_pvalue <- first_taxon_result["pNIE"]
  # })
  # print(paste("Medzim completed in:", medzim_time[3], "seconds"))
  
  loop_end_time <- Sys.time()  # End timing for the entire loop
  print(paste("Loop", i, "completed in:", loop_end_time - loop_start_time, "mins"))
  
  # Store results
  results[[i]] <- list(
    modima_p_value = modima_p_value,
    medtest_p_value = medtest_p_value,
    ccmm_CI_lower = ccmm_CI_lower,
    ccmm_CI_upper = ccmm_CI_upper,
    sparsemcmm_ome_p_value = sparsemcmm_ome_p_value,
    ldm_p_value = ldm_p_value,
    permanova_p_value = permanova_p_value,
    microbvs_result = microbvs_result
    
  )
}

# Convert results to a data frame
results_df <- do.call(rbind, lapply(results, function(x) {
  c(modima_p_value = x$modima_p_value,
    medtest_p_value = x$medtest_p_value,
    ccmm_CI_lower = x$ccmm_CI_lower,
    ccmm_CI_upper = x$ccmm_CI_upper,
    sparsemcmm_ome_p_value = x$sparsemcmm_ome_p_value,
    ldm_p_value = x$ldm_p_value,
    permanova_p_value = x$permanova_p_value,
    microbvs_result = x$microbvs_result
    # MedZIM_pvalue = x$MedZIM_pvalue
    )
}))

# Convert the list of results into a data frame
results_df <- as.data.frame(results_df)

# Save the data frame to a CSV file
write.csv(results_df, file = "mediation_analysis_results_new_no_mediation_effect_sample.csv", row.names = FALSE)

print("Results have been saved")
#####################################################################################3

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Read the saved CSV file
results_df <- read.csv("mediation_analysis_results_sample_new.csv")
#colnames(results_df)
# Add sample size information
results_df$sample_size <- rep(c(50, 100, 300), each = 10)

# Replace missing values with 1 (acceptance)
results_df[is.na(results_df)] <- 1

# Rename methods for better readability
results_long <- melt(results_df, id.vars = "sample_size", 
                     measure.vars = c("modima_p_value", "medtest_p_value", 
                                      "sparsemcmm_ome_p_value.OME", "ldm_p_value", 
                                      "permanova_p_value"
                                      # , "MedZIM_pvalue.pNIE"
                                      ),
                     #MedZIM_pvalue
                     variable.name = "Method", value.name = "p_value")

results_long$Method <- recode(results_long$Method,
                              "modima_p_value" = "MODIMA",
                              "medtest_p_value" = "MedTest",
                              "sparsemcmm_ome_p_value.OME" = "SparseMCMM",
                              "ldm_p_value" = "LDM",
                              "permanova_p_value" = "PERMANOVA"
                              # ,
                              # "MedZIM_pvalue.pNIE" = "MedZIM"
                              )

# Calculate power (p-value < 0.05)
power_results <- results_long %>%
  group_by(sample_size, Method) %>%
  summarise(power = mean(p_value < 0.05, na.rm = TRUE))

# Calculate rejection rate based on confidence intervals (reject if not containing 0)
results_df <- results_df %>%
  mutate(ccmm_reject = ifelse(ccmm_CI_lower > 0 | ccmm_CI_upper < 0, 1, 0),
         microbvs_reject = ifelse(`microbvs_result.2.5.` > 0 | `microbvs_result.97.5.` < 0, 1, 0))

# Add ccmm_reject and microbvs_reject to power_results
ccmm_power_results <- results_df %>%
  group_by(sample_size) %>%
  summarise(Method = "CCMM", power = mean(ccmm_reject, na.rm = TRUE))

microbvs_power_results <- results_df %>%
  group_by(sample_size) %>%
  summarise(Method = "MicroBVS", power = mean(microbvs_reject, na.rm = TRUE))

# Combine all methods' power_results
power_results <- bind_rows(power_results, ccmm_power_results, microbvs_power_results)

# Replace zero power values with a small value
power_results$power[power_results$power == 0] <- 0.001

# Power plot
ggplot(power_results, aes(x = Method, y = power, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Power Analysis",
       x = "Method",
       y = "Power") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample Size")

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Read the saved CSV file
results_df <- read.csv("mediation_analysis_results_sample.csv")
#colnames(results_df)
# Add sample size information
results_df$sample_size <- rep(c(50, 100, 300), each = 20)

# Replace missing values with 1 (acceptance)
results_df[is.na(results_df)] <- 1

# Rename methods for better readability
results_long <- melt(results_df, id.vars = "sample_size", 
                     measure.vars = c("modima_p_value", "medtest_p_value", 
                                      "sparsemcmm_ome_p_value.OME", "ldm_p_value", 
                                      "permanova_p_value"
                                      # , "MedZIM_pvalue.pNIE"
                     ),
                     #MedZIM_pvalue
                     variable.name = "Method", value.name = "p_value")

results_long$Method <- recode(results_long$Method,
                              "modima_p_value" = "MODIMA",
                              "medtest_p_value" = "MedTest",
                              "sparsemcmm_ome_p_value.OME" = "SparseMCMM",
                              "ldm_p_value" = "LDM",
                              "permanova_p_value" = "PERMANOVA"
                              # ,
                              # "MedZIM_pvalue.pNIE" = "MedZIM"
)

# Calculate power (p-value < 0.05)
power_results <- results_long %>%
  group_by(sample_size, Method) %>%
  summarise(power = mean(p_value < 0.05, na.rm = TRUE))

# Calculate rejection rate based on confidence intervals (reject if not containing 0)
results_df <- results_df %>%
  mutate(ccmm_reject = ifelse(ccmm_CI_lower > 0 | ccmm_CI_upper < 0, 1, 0),
         microbvs_reject = ifelse(`microbvs_result.2.5.` > 0 | `microbvs_result.97.5.` < 0, 1, 0))

# Add ccmm_reject and microbvs_reject to power_results
ccmm_power_results <- results_df %>%
  group_by(sample_size) %>%
  summarise(Method = "CCMM", power = mean(ccmm_reject, na.rm = TRUE))

microbvs_power_results <- results_df %>%
  group_by(sample_size) %>%
  summarise(Method = "MicroBVS", power = mean(microbvs_reject, na.rm = TRUE))

# Combine all methods' power_results
power_results <- bind_rows(power_results, ccmm_power_results, microbvs_power_results)

# Replace zero power values with a small value
power_results$power[power_results$power == 0] <- 0.001

# Power plot
ggplot(power_results, aes(x = Method, y = power, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Power Analysis",
       x = "Method",
       y = "Power") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample Size")
############### type-1 error
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Read the saved CSV file for no mediation effect
results_no_mediation_df <- read.csv("mediation_analysis_results_new_no_mediation_effect_sample.csv")

# Add sample size information
results_no_mediation_df$sample_size <- rep(c(50, 100, 300), each = 10)

# Replace missing values with 1 (acceptance)
results_no_mediation_df[is.na(results_no_mediation_df)] <- 1

# Rename methods for better readability
results_no_mediation_long <- melt(results_no_mediation_df, id.vars = "sample_size", 
                                  measure.vars = c("modima_p_value", "medtest_p_value", 
                                                   "sparsemcmm_ome_p_value.OME", "ldm_p_value", 
                                                   "permanova_p_value"
                                                   # , "MedZIM_pvalue.pNIE"
                                                   ),
                                  variable.name = "Method", value.name = "p_value")

results_no_mediation_long$Method <- recode(results_no_mediation_long$Method,
                                           "modima_p_value" = "MODIMA",
                                           "medtest_p_value" = "MedTest",
                                           "sparsemcmm_ome_p_value.OME" = "SparseMCMM",
                                           "ldm_p_value" = "LDM",
                                           "permanova_p_value" = "PERMANOVA"
                                           # ,
                                           # "MedZIM_pvalue.pNIE" = "MedZIM"
                                           )

# Calculate Type-1 error (p-value < 0.05)
type1_error_results <- results_no_mediation_long %>%
  group_by(sample_size, Method) %>%
  summarise(type1_error = mean(p_value < 0.05, na.rm = TRUE))

# Calculate rejection rate based on confidence intervals (reject if not containing 0)
results_no_mediation_df <- results_no_mediation_df %>%
  mutate(ccmm_reject = ifelse(ccmm_CI_lower > 0 | ccmm_CI_upper < 0, 1, 0),
         microbvs_reject = ifelse(`microbvs_result.2.5.` > 0 | `microbvs_result.97.5.` < 0, 1, 0))

# Add ccmm_reject and microbvs_reject to type1_error_results
ccmm_type1_error_results <- results_no_mediation_df %>%
  group_by(sample_size) %>%
  summarise(Method = "CCMM", type1_error = mean(ccmm_reject, na.rm = TRUE))

microbvs_type1_error_results <- results_no_mediation_df %>%
  group_by(sample_size) %>%
  summarise(Method = "MicroBVS", type1_error = mean(microbvs_reject, na.rm = TRUE))

# Combine all methods' type1_error_results
type1_error_results <- bind_rows(type1_error_results, ccmm_type1_error_results, microbvs_type1_error_results)

# Replace zero type1_error values with a small value
type1_error_results$type1_error[type1_error_results$type1_error == 0] <- 0.001

# Type-1 error plot with a horizontal line at 0.05
ggplot(type1_error_results, aes(x = Method, y = type1_error, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Type-1 Error Analysis",
       x = "Method",
       y = "Type-1 Error") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample Size")
##########
library(ggplot2)
library(reshape2)
library(dplyr)

# Read the saved CSV file for no mediation effect
results_no_mediation_df <- read.csv("mediation_analysis_results_no_mediation_effect__sample.csv")

# Add sample size information
results_no_mediation_df$sample_size <- rep(c(50, 100, 300), each = 20)

# Replace missing values with 1 (acceptance)
results_no_mediation_df[is.na(results_no_mediation_df)] <- 1

# Rename methods for better readability
results_no_mediation_long <- melt(results_no_mediation_df, id.vars = "sample_size", 
                                  measure.vars = c("modima_p_value", "medtest_p_value", 
                                                   "sparsemcmm_ome_p_value.OME", "ldm_p_value", 
                                                   "permanova_p_value"
                                                   # , "MedZIM_pvalue.pNIE"
                                  ),
                                  variable.name = "Method", value.name = "p_value")

results_no_mediation_long$Method <- recode(results_no_mediation_long$Method,
                                           "modima_p_value" = "MODIMA",
                                           "medtest_p_value" = "MedTest",
                                           "sparsemcmm_ome_p_value.OME" = "SparseMCMM",
                                           "ldm_p_value" = "LDM",
                                           "permanova_p_value" = "PERMANOVA"
                                           # ,
                                           # "MedZIM_pvalue.pNIE" = "MedZIM"
)

# Calculate Type-1 error (p-value < 0.05)
type1_error_results <- results_no_mediation_long %>%
  group_by(sample_size, Method) %>%
  summarise(type1_error = mean(p_value < 0.05, na.rm = TRUE))

# Calculate rejection rate based on confidence intervals (reject if not containing 0)
results_no_mediation_df <- results_no_mediation_df %>%
  mutate(ccmm_reject = ifelse(ccmm_CI_lower > 0 | ccmm_CI_upper < 0, 1, 0),
         microbvs_reject = ifelse(`microbvs_result.2.5.` > 0 | `microbvs_result.97.5.` < 0, 1, 0))

# Add ccmm_reject and microbvs_reject to type1_error_results
ccmm_type1_error_results <- results_no_mediation_df %>%
  group_by(sample_size) %>%
  summarise(Method = "CCMM", type1_error = mean(ccmm_reject, na.rm = TRUE))

microbvs_type1_error_results <- results_no_mediation_df %>%
  group_by(sample_size) %>%
  summarise(Method = "MicroBVS", type1_error = mean(microbvs_reject, na.rm = TRUE))

# Combine all methods' type1_error_results
type1_error_results <- bind_rows(type1_error_results, ccmm_type1_error_results, microbvs_type1_error_results)

# Replace zero type1_error values with a small value
type1_error_results$type1_error[type1_error_results$type1_error == 0] <- 0.001

# Type-1 error plot with a horizontal line at 0.05
ggplot(type1_error_results, aes(x = Method, y = type1_error, fill = as.factor(sample_size))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Type-1 Error Analysis",
       x = "Method",
       y = "Type-1 Error") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample Size")

#########33

