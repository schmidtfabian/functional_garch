# --------------------------------------------------
# Function: baseline_compute()
#
# Description:
#   Computes truncated functional principal components from squared curves.
#   Demeans the curves, constructs covariance matrix, 
#   and returns truncated eigenvectors, eigenvalues, and explained variance.
#
# Args:
#   y_squared      : matrix (observations x grid points) of squared curves
#
# Returns:
#   list containing:
#     eigenvectors          : matrix of functional PCs (columns)
#     eigenvalues           : vector of eigenvalues
#     explained_variance    : vector of proportion of variance explained by each PC
# --------------------------------------------------

baseline_compute <- function(y_squared) {
  
  n_obs <- nrow(y_squared)
  T_grid <- ncol(y_squared)
  dt <- 1 / (T_grid - 1)   # compute dt automatically assuming equidistant grid
  
  # --------------------------------------------------
  # Demean y^2
  # --------------------------------------------------
  mean_function <- colMeans(y_squared)
  y_squared_demeaned <- y_squared - matrix(mean_function, nrow = n_obs, ncol = T_grid, byrow = TRUE)
  
  # --------------------------------------------------
  # Covariance matrix
  # --------------------------------------------------
  cov_matrix <- t(y_squared_demeaned) %*% y_squared_demeaned * (1/n_obs) * dt
  
  # --------------------------------------------------
  # Eigen decomposition
  # --------------------------------------------------
  cov_eigen <- eigen(cov_matrix)
  
  eigenvectors <- cov_eigen$vectors
  eigenvectors <- eigenvectors / sqrt(dt)
  eigenvalues  <- cov_eigen$values
  
  # Change sign of eigenvectors if the first PC is all negative
  if (mean(eigenvectors[,1])<0) {
    eigenvectors <- -eigenvectors
  }
  
  # Ensure positivity of each PC
  for (k in 1:ncol(eigenvectors)) {
    min_value <- min(eigenvectors[,k])
    shift <- min(min_value, 0)
    eigenvectors[,k] <- eigenvectors[,k] - shift
  }
  
  # --------------------------------------------------
  # Explained variance for each PC
  # --------------------------------------------------
  explained_variance <- eigenvalues / sum(eigenvalues)
  
  # --------------------------------------------------
  # Return results
  # --------------------------------------------------
  return(list(
    eigenvectors = eigenvectors,
    eigenvalues = eigenvalues,
    explained_variance = explained_variance
  ))
}
