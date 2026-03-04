library(nloptr)

# --------------------------------------------------
# Function: fgarch_est_qml()
#
# Description:
#   Estimates Functional GARCH(1,1) model using QML.
#
# Args:
#   y_matrix_squared : matrix (n x T)
#   PC               : vector or matrix (T x M)
#
# Returns:
#   list with estimated delta, alpha, beta and nloptr output
# --------------------------------------------------

fgarch_est_qml <- function(y_matrix_squared, PC) {
  
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
  scores <- y_matrix_squared %*% PC * dt
  n <- nrow(scores)
  M <- ncol(scores)
  
  # --------------------------------------------------
  # Conditional variance recursion
  # --------------------------------------------------
  garch_variance <- function(delta, A, B, y, n, M, T, PC, PC_t) {
    
    sigma <- matrix(NA_real_, nrow = n, ncol = T)
    
    delta_vector <- as.vector(PC%*%delta)
    alpha_kernel_matrix <- PC %*% A %*% PC_t
    beta_kernel_matrix <- PC %*% B %*% PC_t
    
    sigma[1, ] <- delta_vector
    
    # remaining steps
    for (i in 2:n) {
      sigma[i, ] <- delta_vector + alpha_kernel_matrix %*% y[i-1, ] *dt +
        beta_kernel_matrix %*% sigma[i-1, ] *dt
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
    sigma <- garch_variance(delta, A, B, y_matrix_squared, n, M, T_grid, PC, PC_t)
    if (any(is.na(sigma)) || any(sigma <= 0)) return(1e10)
    
    y_trim <- scores[-1, ]
    sigma_scores <- sigma %*% PC_dt
    sigma_scores_trim <- sigma_scores[-1,]
    lik <- (y_trim / sigma_scores_trim) + log(sigma_scores_trim)
    return(mean(rowSums(lik)))
    
  }
  
  # ------------------------------------------------------------------
  # Stationarity constraint on kernel matrices
  # ------------------------------------------------------------------
  stationarity_constraint <- function(params) {
    A <- matrix(params[(M + 1):(M + M^2)], nrow = M)
    B <- matrix(params[(M + M^2 + 1):(M + 2*M^2)], nrow = M)
    alpha_kernel_matrix <- PC %*% A %*% PC_t
    beta_kernel_matrix <- PC %*% B %*% PC_t
    norm_operators <- sqrt(sum((alpha_kernel_matrix + beta_kernel_matrix)^2) * dt^2)
    return(norm_operators - 1)
  }
  
  
  # ------------------------------------------------------------------
  # Initialization
  # ------------------------------------------------------------------
  
  norms <- sqrt(dt * colSums(PC^2))
  max_norm <- max(norms)
  
  restrict_beta <- 0.99 * (M^2 * max_norm)^(-1)
  init_delta <- rep(0.01/M, M)
  init_alpha <- matrix(0.3/M, M,M)
  init_beta  <- matrix(0.1*restrict_beta, M,M)               
  
  # Flatten
  init_theta <- c(init_delta, as.vector(init_alpha), as.vector(init_beta))
  
  # ------------------------------------------------------------------
  # Optimization
  # ------------------------------------------------------------------
  optim_res <- nloptr(
    x0 = init_theta,
    eval_f = obj_wrapper,
    lb = c(rep(1e-4,M),rep(0,length(init_theta)-M)),
    ub = c(rep(1,M+M^2), rep(restrict_beta,M^2)),
    eval_g_ineq = stationarity_constraint,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 1e-6,
      maxeval = 1000
    )
  )
  
  # ------------------------------------------------------------------
  # Reconstruction
  # ------------------------------------------------------------------
  theta_hat <- optim_res$solution
  
  B <- matrix(theta_hat[(M + M^2 + 1):(M + 2 * M^2)], nrow = M, ncol = M)
  A <- matrix(theta_hat[(M + 1):(M + M^2)], nrow = M, ncol = M)
  
  
  return(list(
    delta = theta_hat[1:M],
    alpha = A,
    beta  = B,
    optim_output = optim_res,
    restrict_beta = restrict_beta
  ))
}


