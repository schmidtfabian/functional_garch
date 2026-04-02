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

fgarch_est_ls <- function(y_matrix_squared, PC, maxeval = 1000,
                          positivity_grid_size = NULL, constrain_positivity = TRUE) {
  
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
  
  # --------------------------------------------------
  # Positivity grid
  # --------------------------------------------------
  if (is.null(positivity_grid_size) || positivity_grid_size >= T_grid) {
    positivity_idx <- seq_len(T_grid)
  } else {
    positivity_idx <- unique(round(seq(1, T_grid, length.out = positivity_grid_size)))
  }
  
  T_pos <- length(positivity_idx)
  PC_pos <- PC[positivity_idx, , drop = FALSE]
  
  # --------------------------------------------------
  # Kronecker products
  # --------------------------------------------------
  K_PC_pos <- kronecker(PC_pos, PC_pos)  # (T_pos^2) x (M^2)
  
  scores <- y_matrix_squared %*% PC * dt
  n <- as.numeric(nrow(scores))
  M <- as.numeric(ncol(scores))
  
  # --------------------------------------------------
  # Conditional variance recursion
  # --------------------------------------------------
  garch_variance <- function(delta, A, B, y, n, M) {
    # we produce s_i for i = 2,...,n  -> total n-1 rows
    sigma <- matrix(NA_real_, nrow = n-1, ncol = M)
    
    sigma_initial <- solve(diag(M) - A - B) %*% delta
    # first recursion step: corresponds to i = 2
    sigma[1, ] <- delta + A %*% as.numeric(y[1, ]) + B %*% sigma_initial
    
    # remaining steps
    if (n > 2) {
      for (i in 2:(n-1)) {
        sigma[i, ] <- delta + A %*% as.numeric(y[i, ]) + B %*% as.numeric(sigma[i-1, ])
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
    
    if (norm(A + B, type = "F") >= 0.9999) {
      return(1e12)
    }
    
    # Calculate Variance
    sigma <- garch_variance(delta, A, B, scores, n, M)
    
    # Calculate Objective
    resid <- scores[-1, , drop = FALSE] - sigma
    return(sum(resid^2))
  }
  
  stationarity_positivity_constraint <- function(params) {
    delta <- params[1:M]
    Avec  <- params[(M + 1):(M + M^2)]
    Bvec  <- params[(M + M^2 + 1):(M + 2 * M^2)]
    
    A <- matrix(Avec, nrow = M, ncol = M)
    B <- matrix(Bvec, nrow = M, ncol = M)
    
    c_norm <- norm(A + B, type = "F") - 0.9999
    if (constrain_positivity){
      delta_vector <- as.vector(PC_pos%*%delta)
      alpha_vec <- as.vector(K_PC_pos %*% Avec)
      beta_vec  <- as.vector(K_PC_pos %*% Bvec)
      
      c_delta <- 1e-4 - delta_vector
      c_alpha <- -alpha_vec
      c_beta  <- -beta_vec
      return(c(c_norm, c_delta, c_alpha, c_beta))
    }else{
      return(c_norm)
    }
  }
  
  # ------------------------------------------------------------------
  # Initialization
  # ------------------------------------------------------------------
  init_delta <- c(0.01, rep(0, M-1))
  init_alpha <- matrix(c(0.4, rep(0,M^2-1)), M,M)
  init_beta  <- matrix(c(0.4, rep(0,M^2-1)), M,M)             
  
  # Flatten
  init_theta <- c(init_delta, as.vector(init_alpha), as.vector(init_beta))
  
  # ------------------------------------------------------------------
  # Optimization
  # ------------------------------------------------------------------
  optim_res <- nloptr::nloptr(
    x0 = init_theta,
    eval_f = obj_wrapper,
    lb = rep(-1, length(init_theta)),
    ub = rep(1, length(init_theta)),
    eval_g_ineq = stationarity_positivity_constraint,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      maxeval = maxeval,
      xtol_rel = 0,
      print_level = 0
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