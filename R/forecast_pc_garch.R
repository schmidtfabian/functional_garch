forecast_pc_garch <- function(lambda_fitted, y_eval, y_train,
                              PC, delta, alpha, beta,
                              quantiles_VAR) {
  # --------------------------------------------------
  # Ensure PC is a matrix
  # --------------------------------------------------
  if (is.vector(PC)) {
    PC <- matrix(PC, ncol = 1)
  }
  
  n <- dim(y_eval)[1]
  T <- dim(y_eval)[2]
  M <- dim(PC)[2]
  q <- dim(quantiles_VAR)[1]
  
  scores_train <- y_train %*%PC
  scores_train_squared <- scores_train^2
  scores <- y_eval %*% PC
  scores_squared <- scores^2
  
  lambda_forecast <- matrix(NA, nrow = n, ncol = M)
  lambda_forecast[1,] <- delta + alpha %*% scores_train_squared[nrow(scores_train_squared),] 
  + beta %*% lambda_fitted[nrow(lambda_fitted),]
  if (n > 2){
    for (i in 2:n) {
      lambda_forecast[i,] <- delta + alpha %*% scores_squared[i-1,] + 
        beta %*% lambda_forecast[i-1,]
    }
  }
  
  H_sqrt_fitted <- array(NA, dim = c(n,T,T))
  for (i in 1:n) {
    Lambda_matrix <- diag(sqrt(lambda_forecast[i,]))
    H_sqrt_fitted[i,,] <- PC %*% Lambda_matrix %*% t(PC)
  }
  y_forecast_VAR <- array(NA,dim = c(q,n,T))
  for (j in 1:q) {
    for (i in 1:n) {
      y_forecast_VAR[j,i,] <- H_sqrt_fitted[i,,]%*%quantiles_VAR[j,]
    }
  }
  return(list(
    y_forecast_VAR = y_forecast_VAR
  ))
}