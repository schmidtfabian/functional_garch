# --------------------------------------------------
# Function: fpca_compute()
#
# Description:
#   Computes functional principal components from squared curves.
#   Demeans the curves, constructs covariance matrix, 
#   and returns eigenvectors, eigenvalues, and explained variance.
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

fpca_compute <- function(y_squared) {
  
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
  eigenvalues  <- cov_eigen$values
  
  # --------------------------------------------------
  # Rescale functional PCs to continuous L2 norm
  # --------------------------------------------------
  eigenvectors <- eigenvectors / sqrt(dt)
  
  # Flip sign of each PC if mean < 0
  for (k in 1:ncol(eigenvectors)) {
    if (mean(eigenvectors[, k]) < 0) {
      eigenvectors[, k] <- -eigenvectors[, k]
    }
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
