calculate_bootstrap <- function(epsilon_fitted, Bootstrap_samples,
                                quantiles_vector){
  N <- dim(epsilon_fitted)[1]
  T_grid <- dim(epsilon_fitted)[2]
  
  bootstrap_indices <- sample(1:N, size = Bootstrap_samples, replace = TRUE)
  
  bootstrap_values <- epsilon_fitted[bootstrap_indices, , drop = FALSE]
  
  quantile_matrix <- matrix(NA, nrow = length(quantiles_vector), ncol = T_grid)
  for (i in 1:T_grid){
    quantile_matrix[,i] <- quantile(bootstrap_values[,i], probs = quantiles_vector)
  }
  return(quantile_matrix)
}