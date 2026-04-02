compute_test_statistic <- function(y_forecast, quantiles_VAR, y_eval,
                                   number_of_simulations = 10000,
                                   delta_sim = NULL, alpha_sim = NULL, beta_sim = NULL,
                                   delta_est = NULL, alpha_est = NULL, beta_est = NULL,
                                   compute_tests_other = TRUE, sigma_train = NULL, sigma_fit = NULL,
                                   sigma_forecast = NULL, sigma_eval = NULL) {
  
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
    mean_indicator[i, ] <- colMeans(indicator_arr[i, , ])
    mean_indicator_minus_quantile[i, ] <-
      mean_indicator[i, ] - quantiles_VAR[i]
    
    norm_sq <- sum(mean_indicator_minus_quantile[i,]^2) * dt
    test_T_n[i] <- n * norm_sq
    
    ## Covariance matrix over t = 1,...,T
    cov_matrix_indicator <- cov(indicator_arr[i, , ])*dt
    eigenvals <- eigen(cov_matrix_indicator, only.values = TRUE, symmetric = TRUE)$values
    
    ## Simulate asymptotic distribution
    sim_stats <- numeric(number_of_simulations)
    for (j in 1:number_of_simulations) {
      sim_stats[j] <- sum(eigenvals * rnorm(T)^2)
    }
    
    ## p-value
    p_values[i] <- mean(sim_stats > test_T_n[i])
    
    ## Total (unconditional) deviation
    total_deviation[i] <-
      mean(indicator_arr[i, , ]) - quantiles_VAR[i]
  }
  
  
  if (compute_tests_other) {
    msd_delta <- sum((delta_est - delta_sim)^2) * dt
    msd_alpha <- (dt * RSpectra::svds(alpha_est - alpha_sim, k = 1, nu = 0, nv = 0)$d[1])^2
    msd_beta  <- (dt * RSpectra::svds(beta_est - beta_sim, k = 1, nu = 0, nv = 0)$d[1])^2
    
    ISE_fit  <- dt * rowSums((sigma_fit - sigma_train)^2)
    RISE_fit <- sqrt(ISE_fit)
    mean_RISE_fit <- mean(RISE_fit)
    
    ISE_eval  <- dt * rowSums((sigma_forecast - sigma_eval)^2)
    RISE_eval <- sqrt(ISE_eval)
    mean_RISE_eval <- mean(RISE_eval)
    
    norm_train <- dt * rowSums(sigma_train^2)
    norm_eval <- dt * rowSums(sigma_eval^2)
    
    rel_L2_fit <- ISE_fit / norm_train
    mean_rel_L2_fit <- mean(rel_L2_fit)
    
    rel_L2_eval <- ISE_eval / norm_eval
    mean_rel_L2_eval <- mean(rel_L2_eval)
    
    return(list(
      p_value = p_values,
      total_deviation = total_deviation,
      msd_delta = msd_delta,
      msd_alpha = msd_alpha,
      msd_beta = msd_beta,
      mean_RISE_fit = mean_RISE_fit,
      mean_RISE_eval = mean_RISE_eval,
      mean_rel_L2_fit = mean_rel_L2_fit,
      mean_rel_L2_eval = mean_rel_L2_eval
    ))
  } else {
    return(list(
      p_value = p_values,
      total_deviation = total_deviation
    ))
  }
  
}