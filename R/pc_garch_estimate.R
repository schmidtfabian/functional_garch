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

pc_garch_estimate <- function(y, PC) {
  
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
    # we produce s_i for i = 2,...,n  -> total n-1 rows
    Lambda <- diag(nrow = M)
    
    # first recursion step: corresponds to i = 2
    Lambda[1,1] <- delta + A %*% as.numeric(y[1, ])
    
    # remaining steps
    if (n > 2) {
      for (i in 2:(n-1)) {
        sigma[i, ] <- delta + A %*% as.numeric(y[i-1, ]) + B %*% as.numeric(sigma[i-1, ])
      }
    }
    
    return(sigma)
  }
}