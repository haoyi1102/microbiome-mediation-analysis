# Function to preprocess a single result
preprocess_single_result <- function(result) {
  Y_vector <- as.numeric(result$Y)
  T_vector <- as.numeric(result$T)
  M_a_matrix <- as.matrix(t(result$absolute_M))
  colnames(M_a_matrix) <- paste0("taxon_", 1:ncol(M_a_matrix))
  M_matrix <- M_a_matrix / rowSums(M_a_matrix)
  
  M_nz_matrix <- apply(M_matrix, 1, function(row) {
    min_value <- min(row[row > 0])
    pseudo_count <- min_value / 100
    row[row == 0] <- pseudo_count
    return(row)
  })
  M_nz_matrix <- t(M_nz_matrix)
  
  list(Y_vector = Y_vector, T_vector = T_vector, M_a_matrix = M_a_matrix, M_matrix = M_matrix, M_nz_matrix = M_nz_matrix, feature = result$feature)
}

# Function to preprocess all results in a category
preprocess_results <- function(results) {
  lapply(results, preprocess_single_result)
}