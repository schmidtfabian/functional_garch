# --------------------------------------------------
# Function: fpca_compute()
#
# Description:
#   Computes functional principal components.
#   Demeans the matrix, constructs covariance matrix, 
#   and returns eigenvectors, eigenvalues, and explained variance.
#
# Args:
#   y     : matrix (observations x grid points)
#
# Returns:
#   list containing:
#     eigenvectors          : matrix of PCs (columns)
#     eigenvalues           : vector of eigenvalues
#     explained_variance    : vector of proportion of variance explained by each PC
# --------------------------------------------------

pca_compute <- function(y) {
  
  n_obs <- nrow(y)
  T_grid <- ncol(y)
  
  # --------------------------------------------------
  # Demean y^2
  # --------------------------------------------------
  mean_function <- colMeans(y)
  y_demeaned <- y - matrix(mean_function, nrow = n_obs, ncol = T_grid, byrow = TRUE)
  
  # --------------------------------------------------
  # Covariance matrix
  # --------------------------------------------------
  cov_matrix <- t(y_demeaned) %*% y_demeaned * (1/n_obs)
  
  # --------------------------------------------------
  # Eigen decomposition
  # --------------------------------------------------
  cov_eigen <- eigen(cov_matrix)
  
  eigenvectors <- cov_eigen$vectors
  eigenvalues  <- cov_eigen$values
  
  # --------------------------------------------------
  # Rescale functional PCs to continuous L2 norm
  # --------------------------------------------------
  eigenvectors <- eigenvectors
  
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
