compute_test_statistic <- function(y_forecast, quantiles_VAR, y_eval,
                                   number_of_simulations = 10000,
                                   delta_sim, alpha_sim, beta_sim,
                                   delta_est, alpha_est, beta_est,
                                   compute_msd = TRUE) {
  
  q <- dim(y_forecast)[1]
  n <- dim(y_forecast)[2]
  T <- dim(y_forecast)[3]
  dt <- 1 / (T - 1)
  
  ## Indicator array: q × n × T
  indicator_arr <- array(NA, dim = c(q, n, T))
  for (i in 1:q) {
    indicator_arr[i, , ] <- y_eval < y_forecast[i, , ]
  }
  
  mean_indicator <- matrix(NA, nrow = q, ncol = T)
  mean_indicator_minus_quantile <- matrix(NA, nrow = q, ncol = T)
  test_T_n <- numeric(q)
  p_values <- numeric(q)
  total_deviation <- numeric(q)
  
  for (i in 1:q) {
    
    ## Mean over cross-section
    mean_indicator[i, ] <- colMeans(indicator_arr[i, , ], na.rm = TRUE)
    mean_indicator_minus_quantile[i, ] <-
      mean_indicator[i, ] - quantiles_VAR[i]
    
    ## Trapezoidal approximation of L2 norm
    diff_vec <- mean_indicator_minus_quantile[i, ]
    norm_sq <- sum((diff_vec[-T]^2 + diff_vec[-1]^2) / 2) * dt
    test_T_n[i] <- n * norm_sq
    
    ## Covariance matrix over t = 1,...,T
    cov_matrix_indicator <- cov(t(indicator_arr[i, , ]), use = "pairwise.complete.obs")
    eigenvals <- eigen(cov_matrix_indicator, only.values = TRUE)$values
    
    ## Simulate asymptotic distribution
    sim_stats <- numeric(number_of_simulations)
    for (j in 1:number_of_simulations) {
      sim_stats[j] <- sum(eigenvals * rnorm(T)^2)
    }
    
    ## p-value
    p_values[i] <- mean(sim_stats > test_T_n[i])
    
    ## Total (unconditional) deviation
    total_deviation[i] <-
      quantiles_VAR[i] - mean(indicator_arr[i, , ], na.rm = TRUE)
  }
  
  if (compute_msd) {
    msd_delta <- sqrt(sum((delta_est - delta_sim)^2) * dt)
    msd_alpha <- (dt * RSpectra::svds(alpha_est - alpha_sim, k = 1, nu = 0, nv = 0)$d[1])^2
    msd_beta  <- (dt * RSpectra::svds(beta_est - beta_sim, k = 1, nu = 0, nv = 0)$d[1])^2
  } else {
    msd_delta <- msd_alpha <- msd_beta <- NA_real_
  }
  
  
  return(list(
    test_statistic = test_T_n,
    p_value = p_values,
    total_deviation = total_deviation,
    msd_delta = msd_delta,
    msd_alpha = msd_alpha,
    msd_beta = msd_beta
  ))
}