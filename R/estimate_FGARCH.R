# --------------------------------------------------
# Function: garch_estimate_scores()
#
# Description:
#   Estimates multivariate GARCH-type model on projected
#   squared functional data scores using either LS or QML.
#
# Args:
#   y_matrix_squared : matrix (n x T)
#   PC               : vector or matrix (T x M)
#   method           : "ls" or "qml"
#
# Returns:
#   list with estimated delta, alpha, beta and optim output
# --------------------------------------------------

garch_estimate_scores <- function(y_matrix_squared, PC, method = "ls") {
  
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
  
  scores <- y_matrix_squared %*% PC * dt
  n <- nrow(scores)
  M <- ncol(scores)
  
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
    
    # --- Constraint Handling ---
    
    norm_B <- norm(B, "F")
    if (norm_B >= 0.999) {
      B_final <- B / (norm_B + 0.001)
    } else {
      B_final <- B
    }
    
    L <- diag(nrow = M)
    L[lower.tri(L)] <- A[lower.tri(A)]
    U <- matrix(0,M,M)
    U[upper.tri(U, diag = TRUE)] <- A[upper.tri(A, diag = TRUE)]
    for (i in 1:M) {
      if (U[i, i] == 0) {
        U[i, i] <- U[i, i] + rnorm(1, mean = 0, sd = 1e-8)
        # ensure still nonzero (rare edge case)
        if (U[i, i] == 0) U[i, i] <- eps
      }
    }
    A_final <- L%*%U
    
    # Calculate Variance
    sigma <- garch_variance(delta, A_final, B_final, scores, n, M)
    
    # Calculate Objective
    if (method == "ls") {
      resid <- scores[-1, ] - sigma
      return(sum(resid^2))
      
    } else if (method == "qml") {
      if (any(sigma <= 1e-8) || any(is.na(sigma)) || any(is.infinite(sigma))) {
        P_CONST <- 1e10
        # Calculate violation magnitude for the slope
        # We sum the absolute magnitude of all negative/zero elements
        bad_vals <- sigma[sigma <= 1e-8]
        # Handle NA/Inf in the sum safely
        bad_vals <- bad_vals[!is.na(bad_vals) & is.finite(bad_vals)]
        
        violation <- sum(abs(bad_vals)) + 0.1 # Add 0.1 to ensure we are above P_CONST
        
        return(P_CONST + 1e5 * violation)
      }
      y_trim <- scores[-1, ]
      lik <- (y_trim / sigma) + log(sigma)
      return(mean(rowSums(lik)))
    }
  }
  
  
  # ------------------------------------------------------------------
  # Initialization
  # ------------------------------------------------------------------
  # Initialize diagonal matrices to ensure full rank and stability, scale by M to ensure norm<1
  init_delta <- rep(0.01/M, M)
  init_alpha <- matrix(0.5, M,M)
  init_beta  <- matrix(0.5, M,M)               
  
  # Flatten
  init_theta <- c(init_delta, as.vector(init_alpha), as.vector(init_beta))
  
  # ------------------------------------------------------------------
  # Optimization
  # ------------------------------------------------------------------
  optim_res <- optim(
    par     = init_theta,
    fn      = obj_wrapper,
    method  = "L-BFGS-B",
    lower   = 1e-8,
    upper   = 1000,
  )
  
  # ------------------------------------------------------------------
  # Reconstruction
  # ------------------------------------------------------------------
  theta_hat <- optim_res$par
  
  B <- matrix(theta_hat[(M + M^2 + 1):(M + 2 * M^2)], nrow = M, ncol = M)
  A <- matrix(theta_hat[(M + 1):(M + M^2)], nrow = M, ncol = M)
  
  norm_B <- norm(B, "F")
  if (norm_B >= 0.999) {
    B_final <- B / (norm_B + 0.001)
  } else {
    B_final <- B
  }
  L <- diag(nrow = M)
  L[lower.tri(L)] <- A[lower.tri(A)]
  U <- matrix(0,M,M)
  U[upper.tri(U, diag = TRUE)] <- A[upper.tri(A, diag = TRUE)]
  A_final <- L%*%U
  
  
  return(list(
    delta = theta_hat[1:M],
    alpha = A_final,
    beta  = B_final,
    optim_output = optim_res
  ))
}


