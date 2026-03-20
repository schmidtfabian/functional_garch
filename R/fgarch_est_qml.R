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

fgarch_est_qml <- function(y_matrix_squared, PC, maxeval = 1000) {
  
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
  scores <- as.matrix(scores)   # safeguard for M = 1
  
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
    
    sigma[1, ] <- solve(diag(T) -
                          (alpha_kernel_matrix + beta_kernel_matrix)*dt)%*% delta_vector
    
    # remaining steps
    if (n > 1) {
      for (i in 2:n) {
        sigma[i, ] <- delta_vector + alpha_kernel_matrix %*% y[i-1, ] *dt +
          beta_kernel_matrix %*% sigma[i-1, ] *dt
      }
    }
    return(sigma)
  }
  
  
  # --------------------------------------------------
  # Objective functions
  # --------------------------------------------------
  
  obj_wrapper <- function(params) {
    
    delta <- params[1:M]
    A <- matrix(params[(M + 1):(M + M^2)], nrow = M, ncol = M)
    B <- matrix(params[(M + M^2 + 1):(M + 2 * M^2)], nrow = M, ncol = M)
    
    sigma <- garch_variance(delta, A, B, y_matrix_squared, n, M, T_grid, PC, PC_t)
    if (any(is.na(sigma)) || any(!is.finite(sigma))) return(1e10)
    
    sigma_scores <- sigma %*% PC_dt
    sigma_scores <- as.matrix(sigma_scores)   # safeguard for M = 1
    
    y_trim <- scores[-1, , drop = FALSE]
    sigma_scores_trim <- sigma_scores[-1, , drop = FALSE]
    
    if (any(is.na(sigma_scores_trim)) || any(!is.finite(sigma_scores_trim)) ||
        any(sigma_scores_trim <= 0)) {
      return(1e8)
    }
    
    lik <- (y_trim / sigma_scores_trim) + log(sigma_scores_trim)
    mean(rowSums(lik))
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
    return(norm_operators - 0.9999)
  }
  
  
  # ------------------------------------------------------------------
  # Initialization
  # ------------------------------------------------------------------
  
  norms <- sqrt(dt * colSums(PC^2))
  max_norm <- max(norms)
  
  restrict_beta <- 0.99 * (M^2 * max_norm)^(-1)
  init_delta <- c(0.01,rep(0, M-1))
  init_alpha <- matrix(0, M,M)
  init_beta  <- matrix(0, M,M)               
  
  # Flatten
  init_theta <- c(init_delta, as.vector(init_alpha), as.vector(init_beta))
  
  # ------------------------------------------------------------------
  # Optimization
  # ------------------------------------------------------------------
  optim_res <- nloptr::nloptr(
    x0 = init_theta,
    eval_f = obj_wrapper,
    lb = c(1e-4,rep(0,length(init_theta)-1)),
    ub = c(rep(1,M+M^2), rep(restrict_beta,M^2)),
    eval_g_ineq = stationarity_constraint,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      maxeval = maxeval
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


