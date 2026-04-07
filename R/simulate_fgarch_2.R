# --------------------------------------------------
# simulate_fgarch_2()
#
# Simulates functional GARCH process with OU-type innovations
#
# Args:
#   burn_in      : integer, burn-in sample size
#   train_size   : integer, training sample size
#   eval_size    : integer, evaluation sample size
#   T_grid_points: integer, number of discretization points
#   epsilon_function: function, innovation/error function
#   seed         : optional integer; if provided, RNG seed is set
#
# Returns:
#   list with:
#     y_train         : training sample matrix
#     y_eval          : evaluation sample matrix
#     sigma_squared   : full volatility paths
#     grid            : discretization grid
#
# Notes:
#   - If seed = NULL, randomness is controlled externally (main.R)
#   - If seed is provided, simulation is reproducible in isolation
# --------------------------------------------------
simulate_fgarch_2 <- function(
    burn_in,
    train_size,
    eval_size,
    T_grid_points,
    epsilon_function,
    seed = NULL,
    ...
) {
  
  # --------------------------------------------------
  # Optional reproducibility
  # --------------------------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # --------------------------------------------------
  # Setup
  # --------------------------------------------------
  N_observations <- burn_in + train_size + eval_size
  
  t_grid <- seq(0, 1, length.out = T_grid_points)
  t_difference <- 1 / (T_grid_points - 1)
  
  delta_function <- function(t) exp(-(t)^2)
  delta_vector <- delta_function(t_grid)
  alpha_function <- function(s, t) 0.5*exp(-(s-t)^2)
  beta_function <- function(s, t) 12*t*(1-t)*s*(1-s)
  
  alpha_kernel_matrix <- outer(t_grid, t_grid, alpha_function)
  beta_kernel_matrix  <- outer(t_grid, t_grid, beta_function)
  
  sigma_squared_matrix <- matrix(0, N_observations, T_grid_points)
  y_matrix <- matrix(0, N_observations, T_grid_points)
  
  sigma_squared_matrix[1, ] <- delta_vector
  
  # --------------------------------------------------
  # Simulation loop
  # --------------------------------------------------
  for (i in 1:N_observations) {
    
    if (i > 1) {
      sigma_squared_matrix[i, ] <- delta_vector +
        as.vector(alpha_kernel_matrix %*% y_matrix[i-1, ]^2) * t_difference +
        as.vector(beta_kernel_matrix  %*% sigma_squared_matrix[i-1, ]) * t_difference
    }
    
    y_matrix[i, ] <- sqrt(sigma_squared_matrix[i, ]) *
      epsilon_function(t_grid, ...)
  }
  
  # --------------------------------------------------
  # Split samples, discard burn-in
  # --------------------------------------------------
  y_eval <- y_matrix[
    (burn_in + train_size + 1):(burn_in + train_size + eval_size), 
  ]
  
  y_train <- y_matrix[
    (burn_in + 1):(burn_in + train_size), 
  ]
  
  # --------------------------------------------------
  # Return results
  # --------------------------------------------------
  return(list(
    y_train = y_train,
    y_eval  = y_eval,
    sigma_squared = sigma_squared_matrix,
    grid = t_grid,
    delta = delta_vector,
    alpha = alpha_kernel_matrix,
    beta = beta_kernel_matrix
  ))
}