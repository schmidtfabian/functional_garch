# Forecast value at risk
#
# This function computes recursive forecasts of conditional variances using
# previously fitted model parameters and then derives VaR forecasts by scaling
# the variance paths with supplied quantiles.
#
# Inputs:
# - y_eval_squared: Numeric matrix (N x T) of squared observations for the
#   evaluation period.
# - y_train_squared: Numeric matrix of squared observations from the training
#   period; the last row is used to initialize the recursion.
# - fitted_values: List containing model parameters and fitted training variance:
#     * sigma_fitted: Matrix of fitted variances from the training sample
#     * delta: Numeric intercept term (scalar or vector over grid points)
#     * alpha: Coefficient matrix for lagged squared observations
#     * beta: Coefficient matrix for lagged conditional variances
# - quantile_matrix: Numeric matrix (Q x T) of quantiles for VaR computation,
#   where Q is the number of quantile levels and T the number of grid points.
# Outputs:
# - Returns a list with:
#     * sigma_fitted: Matrix (N x T) of forecasted conditional variances
#     * VAR_forecast: Array (Q x N x T) of VaR forecasts for each quantile level
forecast_VAR <- function(y_eval_squared, y_train_squared, fitted_values, quantile_matrix){
  # --------------------------------------------------
  # Setup
  # --------------------------------------------------
  sigma_train <- fitted_values$sigma_fitted
  delta <- fitted_values$delta
  alpha <- fitted_values$alpha
  beta <- fitted_values$beta
  
  T_grid_points <- dim(y_eval_squared)[2]
  N <- dim(y_eval_squared)[1]
  dt <- 1/(T_grid_points-1)
  
  # --------------------------------------------------
  # Recursion
  # --------------------------------------------------
  
  sigma_squared <- matrix(data = NA, nrow = N, ncol = T_grid_points)
  sigma_squared[1,] <- delta +
    as.vector(alpha %*% y_train_squared[nrow(y_train_squared),])*dt +
    as.vector(beta %*% sigma_train[nrow(sigma_train),])*dt
  for (i in 2:N){
    sigma_squared[i,] <- delta +
      as.vector(alpha %*% y_eval_squared[i-1,])*dt +
      as.vector(beta %*% sigma_squared[i-1, ])*dt
  }

  # --------------------------------------------------
  # VAR forecasting
  # --------------------------------------------------
  
  Q <- dim(quantile_matrix)[1]
  VAR_forecast <- array(NA, dim = c(Q,N,T_grid_points))
  for (q in 1:Q){
    VAR_forecast[q,,] <- quantile_matrix[q,]*sqrt(sigma_squared)
  }
  
  # --------------------------------------------------
  # Return results
  # --------------------------------------------------
  return(list(
    sigma_forecast = sigma_squared,
    VAR_forecast = VAR_forecast
  ))
}