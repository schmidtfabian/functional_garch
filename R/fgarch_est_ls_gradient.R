# --------------------------------------------------
# Function: fgarch_est_ls()
#
# Description:
#   Estimates Functional GARCH(1,1) model using LS
#   with
#     (i) positivity constraints imposed on a coarse grid
#     (ii) smooth stationarity constraint
#          ||A + B||_F^2 <= tau^2
#     (iii) explicit Jacobians for all constraints
#     (iv) numerical gradient for the objective
#
# Args:
#   y_matrix_squared     : matrix (n x T_grid)
#   PC                   : matrix (T_grid x M)
#   maxeval              : maximum number of optimizer evaluations
#   tau                  : stationarity threshold, 0 < tau < 1
#   positivity_grid_size : number of grid points used for positivity
#                          constraints; if NULL or >= T_grid, the full
#                          grid is used
#   lb_theta             : lower bound for parameters
#   ub_theta             : upper bound for parameters
#   print_level          : nloptr print level
#
# Returns:
#   list with estimated delta, alpha, beta,
#   reconstructed delta_vec, alpha_kernel, beta_kernel on full grid,
#   coarse-grid objects, optimizer output, and constraint info
# --------------------------------------------------

fgarch_est_ls_test <- function(y_matrix_squared,
                          PC,
                          maxeval = 1000,
                          tau = 0.9999,
                          positivity_grid_size = 100,
                          lb_theta = -1,
                          ub_theta = 1,
                          print_level = 0) {
  
  # --------------------------------------------------
  # Ensure PC is a matrix
  # --------------------------------------------------
  if (is.vector(PC)) {
    PC <- matrix(PC, ncol = 1)
  }
  
  # --------------------------------------------------
  # Basic dimensions
  # --------------------------------------------------
  n <- nrow(y_matrix_squared)
  T_grid <- ncol(y_matrix_squared)
  M <- ncol(PC)
  
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
  # Projection step on full grid
  # --------------------------------------------------
  dt <- 1 / (T_grid - 1)
  PC_dt <- PC * dt
  scores <- y_matrix_squared %*% PC_dt   # n x M
  
  # --------------------------------------------------
  # Kronecker products
  # --------------------------------------------------
  K_PC_pos <- kronecker(PC_pos, PC_pos)  # (T_pos^2) x (M^2)
  
  # --------------------------------------------------
  # Helper: unpack parameter vector
  # --------------------------------------------------
  unpack_params <- function(params, M) {
    delta <- params[1:M]
    
    A <- matrix(
      params[(M + 1):(M + M^2)],
      nrow = M,
      ncol = M
    )
    
    B <- matrix(
      params[(M + M^2 + 1):(M + 2 * M^2)],
      nrow = M,
      ncol = M
    )
    
    list(delta = delta, A = A, B = B)
  }
  
  # --------------------------------------------------
  # Conditional variance recursion
  # --------------------------------------------------
  garch_variance <- function(delta, A, B, y, n, M) {
    sigma <- matrix(NA_real_, nrow = n - 1, ncol = M)
    
    init_mat <- diag(M) - A - B
    sigma_initial <- solve(init_mat, delta)
    
    sigma[1, ] <- as.numeric(delta + A %*% y[1, ] + B %*% sigma_initial)
    
    if (n > 2) {
      for (i in 2:(n - 1)) {
        sigma[i, ] <- as.numeric(delta + A %*% y[i, ] + B %*% sigma[i - 1, ])
      }
    }
    
    sigma
  }
  
  # --------------------------------------------------
  # Objective function
  # --------------------------------------------------
  obj_wrapper <- function(params) {
    pars <- unpack_params(params, M)
    delta <- pars$delta
    A <- pars$A
    B <- pars$B
    
    init_mat <- diag(M) - A - B
    if (rcond(init_mat) < 1e-12) {
      return(1e12)
    }
    
    sigma <- tryCatch(
      garch_variance(delta, A, B, scores, n, M),
      error = function(e) NULL
    )
    
    if (is.null(sigma) || any(!is.finite(sigma))) {
      return(1e12)
    }
    
    resid <- scores[-1, , drop = FALSE] - sigma
    sum(resid^2)
  }
  
  # --------------------------------------------------
  # Numerical gradient of objective
  # --------------------------------------------------
  obj_grad_wrapper <- function(params) {
    numDeriv::grad(
      func = obj_wrapper,
      x = params,
      method = "simple"
    )
  }
  
  # --------------------------------------------------
  # Positivity constraints on coarse grid
  #
  # nloptr requires g(x) <= 0
  # so positivity is imposed as negatives <= 0
  # --------------------------------------------------
  positivity_constraint <- function(params) {
    pars <- unpack_params(params, M)
    delta <- pars$delta
    A <- pars$A
    B <- pars$B
    
    delta_vec_pos <- as.vector(PC_pos %*% delta)
    alpha_kernel_pos <- PC_pos %*% A %*% t(PC_pos)
    beta_kernel_pos  <- PC_pos %*% B %*% t(PC_pos)
    
    c(
      -delta_vec_pos,
      -as.vector(alpha_kernel_pos),
      -as.vector(beta_kernel_pos)
    )
  }
  
  # --------------------------------------------------
  # Constant Jacobian of positivity constraints
  # (This follows from the Kronecker representation)
  # --------------------------------------------------
  build_positivity_jacobian <- function(PC_pos, K_PC_pos) {
    T_pos <- nrow(PC_pos)
    M <- ncol(PC_pos)
    
    n_constr <- T_pos + 2 * T_pos^2
    n_param <- M + 2 * M^2
    
    J <- matrix(0, nrow = n_constr, ncol = n_param)
    
    # delta block
    J[1:T_pos, 1:M] <- -PC_pos
    
    # alpha block
    row_alpha <- (T_pos + 1):(T_pos + T_pos^2)
    col_alpha <- (M + 1):(M + M^2)
    J[row_alpha, col_alpha] <- -K_PC_pos
    
    # beta block
    row_beta <- (T_pos + T_pos^2 + 1):(T_pos + 2 * T_pos^2)
    col_beta <- (M + M^2 + 1):(M + 2 * M^2)
    J[row_beta, col_beta] <- -K_PC_pos
    
    J
  }
  
  positivity_jacobian_const <- build_positivity_jacobian(PC_pos, K_PC_pos)
  
  positivity_constraint_jac <- function(params) {
    positivity_jacobian_const
  }
  
  # --------------------------------------------------
  # Smooth stationarity constraint:
  #   ||A + B||_F^2 - tau^2 <= 0
  # A smoothed constraint is used so that the gradient is simpler to compute
  # --------------------------------------------------
  stationarity_constraint <- function(params) {
    pars <- unpack_params(params, M)
    S <- pars$A + pars$B
    sum(S^2) - tau^2
  }
  
  # --------------------------------------------------
  # Gradient of stationarity constraint
  # --------------------------------------------------
  stationarity_constraint_jac <- function(params) {
    pars <- unpack_params(params, M)
    S <- pars$A + pars$B
    
    grad <- matrix(0, nrow = 1, ncol = M + 2 * M^2)
    block <- as.vector(2 * S)
    
    grad[1, (M + 1):(M + M^2)] <- block
    grad[1, (M + M^2 + 1):(M + 2 * M^2)] <- block
    
    grad
  }
  
  # --------------------------------------------------
  # Combined constraints and Jacobian
  # --------------------------------------------------
  eval_g_ineq <- function(params) {
    c(
      stationarity_constraint(params),
      positivity_constraint(params)
    )
  }
  
  eval_jac_g_ineq <- function(params) {
    rbind(
      stationarity_constraint_jac(params),
      positivity_constraint_jac(params)
    )
  }
  
  # --------------------------------------------------
  # Initialization
  # --------------------------------------------------
  init_delta <- c(0.01, rep(0, M - 1))
  init_alpha <- matrix(0, nrow = M, ncol = M)
  init_beta  <- matrix(0, nrow = M, ncol = M)
  
  init_alpha[1, 1] <- 0.2
  init_beta[1, 1]  <- 0.2
  
  init_theta <- c(init_delta, as.vector(init_alpha), as.vector(init_beta))
  
  # --------------------------------------------------
  # Bounds
  # --------------------------------------------------
  n_param <- length(init_theta)
  lb <- rep(lb_theta, n_param)
  ub <- rep(ub_theta, n_param)
  
  # --------------------------------------------------
  # Optimization
  # --------------------------------------------------
  optim_res <- nloptr::nloptr(
    x0 = init_theta,
    eval_f = obj_wrapper,
    eval_grad_f = obj_grad_wrapper,
    eval_g_ineq = eval_g_ineq,
    eval_jac_g_ineq = eval_jac_g_ineq,
    lb = lb,
    ub = ub,
    opts = list(
      algorithm = "NLOPT_LD_SLSQP",
      maxeval = maxeval,
      xtol_rel = 1e-8,
      print_level = print_level
    )
  )
  
  # --------------------------------------------------
  # Reconstruction on full grid
  # --------------------------------------------------
  theta_hat <- optim_res$solution
  pars_hat <- unpack_params(theta_hat, M)
  
  delta_hat <- pars_hat$delta
  alpha_hat <- pars_hat$A
  beta_hat  <- pars_hat$B
  
  delta_vec_full <- as.vector(PC %*% delta_hat)
  alpha_kernel_full <- PC %*% alpha_hat %*% t(PC)
  beta_kernel_full  <- PC %*% beta_hat %*% t(PC)
  
  # --------------------------------------------------
  # Reconstruction on positivity grid
  # --------------------------------------------------
  delta_vec_pos <- as.vector(PC_pos %*% delta_hat)
  alpha_kernel_pos <- PC_pos %*% alpha_hat %*% t(PC_pos)
  beta_kernel_pos  <- PC_pos %*% beta_hat %*% t(PC_pos)
  
  # --------------------------------------------------
  # Return
  # --------------------------------------------------
  list(
    delta = delta_hat,
    alpha = alpha_hat,
    beta  = beta_hat,
    
    delta_vec = delta_vec_full,
    alpha_kernel = alpha_kernel_full,
    beta_kernel  = beta_kernel_full,
    
    delta_vec_pos = delta_vec_pos,
    alpha_kernel_pos = alpha_kernel_pos,
    beta_kernel_pos  = beta_kernel_pos,
    
    positivity_idx = positivity_idx,
    positivity_grid_size_used = T_pos,
    positivity_jacobian = positivity_jacobian_const,
    
    min_delta_full = min(delta_vec_full),
    min_alpha_full = min(alpha_kernel_full),
    min_beta_full  = min(beta_kernel_full),
    
    min_delta_pos = min(delta_vec_pos),
    min_alpha_pos = min(alpha_kernel_pos),
    min_beta_pos  = min(beta_kernel_pos),
    
    optim_output = optim_res
  )
}