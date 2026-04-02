# --------------------------------------------------
# Function: pc_garch_estimate()
#
# Description:
#   Estimates multivariate PC-GARCH model using QML.
#
# Args:
#   y_matrix_squared : matrix (n x T)
#   PC               : vector or matrix (T x M)
#
# Returns:
#   list with estimated delta, alpha, beta and optim output
# --------------------------------------------------

pc_garch_estimate <- function(y, PC, maxeval = 1000) {
  
  # --------------------------------------------------
  # Ensure PC is a matrix
  # --------------------------------------------------
  if (is.vector(PC)) {
    PC <- matrix(PC, ncol = 1)
  }
  scores <- y %*% PC
  scores_squared <- scores^2
  
  n <- nrow(scores)
  M <- ncol(scores)
  
  garch_variance <- function(delta, A, B, y, n, M) {
    lambda <- matrix(NA_real_, nrow = n, ncol = M)
    
    lambda[1,] <- solve(diag(nrow = M) - A - B, delta)
    
    # remaining steps
    if (n > 1) {
      for (i in 2:(n)) {
        lambda[i,] <- delta + A %*% as.vector(y[i-1, ]) + B %*% as.vector(lambda[i-1, ])
      }
    }
    
    return(lambda)
  }
  eval_f <- function(params) {
    
    # Unpack
    delta <- params[1:M]
    a_diag <- params[(M + 1):(2 * M)]
    b_diag <- params[(2 * M + 1):(3 * M)]
    
    A <- diag(a_diag, nrow = M)
    B <- diag(b_diag, nrow = M)
    
    
    # Calculate Variance
    lambda <- garch_variance(delta, A, B, scores_squared, n, M)
    if (any(!is.finite(lambda)) || any(lambda <= 1e-8)) {
      return(1e10)
    }
    lik <- (scores_squared / lambda) + log(lambda)
    return(mean(rowSums(lik)))
  }
  
  eval_g_ineq <- function(params) {
    a_diag <- params[(M + 1):(2 * M)]
    b_diag <- params[(2 * M + 1):(3 * M)]
    
    a_diag + b_diag - (1 - 1e-6)
  }
  
  init_delta <- rep(0.01, M)
  init_alpha <- rep(0.1,M)
  init_beta  <- rep(0.1,M)               
  
  # Flatten
  init_theta <- c(init_delta, init_alpha, init_beta)
  
  optim_res <- nloptr::nloptr(
    x0          = init_theta,
    eval_f      = eval_f,
    eval_g_ineq = eval_g_ineq,
    lb          = rep(0, 3 * M),
    ub          = rep(1, 3 * M),
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel  = 0,
      maxeval   = maxeval
    )
  )
  
  theta_hat <- optim_res$solution
  
  delta  <- theta_hat[1:M]
  a_diag <- theta_hat[(M + 1):(2 * M)]
  b_diag <- theta_hat[(2 * M + 1):(3 * M)]
  
  alpha <- diag(a_diag, nrow = M)
  beta  <- diag(b_diag, nrow = M)
  
  lambda <- garch_variance(delta, alpha, beta, scores_squared, n, M)
  
  T <- nrow(PC)
  H_sqrt_fitted <- array(NA_real_, dim = c(n, T, T))
  
  for (i in 1:n) {
    Lambda_matrix <- diag(1 / sqrt(lambda[i, ]), nrow = M)
    H_sqrt_fitted[i, , ] <- PC %*% Lambda_matrix %*% t(PC)
  }
  
  return(list(
    delta = delta,
    alpha = alpha,
    beta = beta,
    optim_res = optim_res,
    H_sqrt_fitted = H_sqrt_fitted,
    lambda_fitted = lambda
  ))
}