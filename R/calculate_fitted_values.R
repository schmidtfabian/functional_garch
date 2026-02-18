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
  delta_vector <- rep(0,T_grid)
  alpha_kernel_matrix <- matrix(0, nrow = T_grid, ncol = T_grid)
  beta_kernel_matrix <- matrix(0, nrow = T_grid, ncol = T_grid)
  for (i in 1:M) {
    delta_vector <- delta_vector+delta[i]*PC[,i]
    for (j in 1:M) {
      alpha_kernel_matrix <- alpha_kernel_matrix + alpha[i,j]*PC[,i]%o%PC[,j]
      beta_kernel_matrix <- beta_kernel_matrix + beta[i,j]*PC[,i]%o%PC[,j]
    }
  }
  sigma_matrix <- matrix(NA, nrow = N, ncol = T_grid)
  sigma_matrix[1,]<- delta_vector+as.vector(alpha_kernel_matrix%*%colMeans(y_matrix_squared))*dt
  for (i in 2:N) {
    sigma_matrix[i,]<- delta_vector + alpha_kernel_matrix%*%y_matrix_squared[i-1,]*dt +
      beta_kernel_matrix%*%sigma_matrix[i-1,]*dt
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