calculate_fitted_values <- function(y_matrix_squared, y_matrix, PC,
                                    delta, alpha, beta) {
  
  # --------------------------------------------------
  # Ensure PC is a matrix
  # --------------------------------------------------
  if (is.vector(PC)) {
    PC <- matrix(PC, ncol = 1)
  }
  M <- dim(PC)[2]
  T_grid <- dim(PC)[1]
  dt <- 1 / (T_grid - 1)
  N <- dim(y_matrix_squared)[1]
  delta_vector <- as.vector(PC %*% delta)
  alpha_kernel_matrix <- PC %*% alpha %*% t(PC)
  beta_kernel_matrix <- PC %*% beta %*% t(PC)
  sigma_matrix <- matrix(NA, nrow = N, ncol = T_grid)
  sigma_matrix[1,]<- delta_vector+as.vector(alpha_kernel_matrix%*%colMeans(y_matrix_squared))*dt
  for (i in 2:N) {
    sigma_matrix[i,]<- delta_vector + as.vector(alpha_kernel_matrix%*%y_matrix_squared[i-1,])*dt +
      as.vector(beta_kernel_matrix%*%sigma_matrix[i-1,])*dt
  }
  epsilon_fitted <- y_matrix/sqrt(sigma_matrix)
  
  return(list(
    epsilon_fitted = epsilon_fitted,
    sigma_fitted = sigma_matrix,
    delta = delta_vector,
    alpha = alpha_kernel_matrix,
    beta = beta_kernel_matrix
  ))
}