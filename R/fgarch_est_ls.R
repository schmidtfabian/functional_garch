library(nloptr)

# --------------------------------------------------
# Function: fgarch_est_ls()
#
# Description:
#   Estimates Functional GARCH(1,1) model using LS.
#
# Args:
#   y_matrix_squared : matrix (n x T)
#   PC               : vector or matrix (T x M)
#
# Returns:
#   list with estimated delta, alpha, beta and nloptr output
# --------------------------------------------------

fgarch_est_ls <- function(y_matrix_squared, PC, check_points = 40, maxeval = 1000) {
  
  # --------------------------------------------------
  # Ensure PC is a matrix
  # --------------------------------------------------
  if (is.vector(PC)) {
    PC <- matrix(PC, ncol = 1)
  }
  
  # --------------------------------------------------
  # Projection step
  # --------------------------------------------------
  T_grid <- ncol(y_matrix_squared)
  dt <- 1 / (T_grid - 1)
  
  PC_t <- t(PC) #avoids computing inside of the loop
  PC_dt <- PC*dt
  
  scores <- y_matrix_squared %*% PC_dt
  n <- as.numeric(nrow(scores))
  M <- as.numeric(ncol(scores))
  
  # select indices for pointwise constraints (sparse)
  k <- min(as.integer(check_points), T_grid)
  idx <- unique(as.integer(round(seq(1, T_grid, length.out = k))))
  
  # --------------------------------------------------
  # Parameter vector structure
  # theta = (delta, vec(alpha), vec(beta))
  # --------------------------------------------------
  
  # --------------------------------------------------
  # Conditional variance recursion
  # --------------------------------------------------
  garch_variance <- function(delta, A, B, y, n, M) {
    # we produce s_i for i = 2,...,n  -> total n-1 rows
    sigma <- matrix(NA_real_, nrow = n-1, ncol = M)
    
    # first recursion step: corresponds to i = 2
    sigma[1, ] <- delta + A %*% as.numeric(y[1, ])
    
    # remaining steps
    if (n > 2) {
      for (i in 2:(n-1)) {
        sigma[i, ] <- delta + A %*% as.numeric(y[i-1, ]) + B %*% as.numeric(sigma[i-1, ])
      }
    }
    
    return(sigma)
  }
  
  
  # --------------------------------------------------
  # Objective functions
  # --------------------------------------------------
  
  obj_wrapper <- function(params) {
    
    # Unpack
    delta <- params[1:M]
    A     <- matrix(params[(M + 1):(M + M^2)], nrow = M, ncol = M)
    B     <- matrix(params[(M + M^2 + 1):(M + 2 * M^2)], nrow = M, ncol = M)
    
    # Calculate Variance
    sigma <- garch_variance(delta, A, B, scores, n, M)
    
    # Calculate Objective
    resid <- scores[-1, ] - sigma
    return(sum(resid^2))
  }
  
  stationarity_positivity_constraint <- function(params) {
    delta <- params[1:M]
    A <- matrix(params[(M + 1):(M + M^2)], nrow = M)
    B <- matrix(params[(M + M^2 + 1):(M + 2*M^2)], nrow = M)
    delta_vector <- as.vector(PC%*%delta)
    alpha_kernel_matrix <- PC %*% A %*% PC_t
    beta_kernel_matrix <- PC %*% B %*% PC_t
    norm_operators <- sqrt(sum((alpha_kernel_matrix + beta_kernel_matrix)^2) * dt^2)
    c_norm <- norm_operators - 1
    # enforce non-negativity only at idx
    c_delta <- -delta_vector[idx]
    c_alpha <- as.vector(-alpha_kernel_matrix[idx, idx])
    c_beta  <- as.vector(-beta_kernel_matrix[idx, idx])
    return(c(c_norm, c_delta, c_alpha, c_beta))
  }
  
  # ------------------------------------------------------------------
  # Initialization
  # ------------------------------------------------------------------
  init_delta <- rep(0.01/M, M)
  init_alpha <- matrix(0.5/M, M,M)
  init_beta  <- matrix(0.5/M, M,M)               
  
  # Flatten
  init_theta <- c(init_delta, as.vector(init_alpha), as.vector(init_beta))
  
  # ------------------------------------------------------------------
  # Optimization
  # ------------------------------------------------------------------
  optim_res <- nloptr(
    x0 = init_theta,
    eval_f = obj_wrapper,
    lb = c(rep(1e-4,M),rep(0,length(init_theta)-M)),
    ub = rep(1, length(init_theta)),
    eval_g_ineq = stationarity_positivity_constraint,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 1e-6,
      maxeval = maxeval
    )
  )
  
  # ------------------------------------------------------------------
  # Reconstruction
  # ------------------------------------------------------------------
  theta_hat <- optim_res$solution
  
  delta <- theta_hat[1:M]
  beta <- matrix(theta_hat[(M + M^2 + 1):(M + 2 * M^2)], nrow = M, ncol = M)
  alpha <- matrix(theta_hat[(M + 1):(M + M^2)], nrow = M, ncol = M)
  
  
  return(list(
    delta = delta,
    alpha = alpha,
    beta  = beta,
    optim_output = optim_res
  ))
}