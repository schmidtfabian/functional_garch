# --------------------------------------------------
# Master Thesis — FGARCH Simulation
# Author: Fabian Christian Schmidt
# --------------------------------------------------

rm(list=ls())

install.packages(setdiff(c("rstudioapi","ggplot2", "dplyr", "tidyr", "nloptr"),
                         rownames(installed.packages())))
library(ggplot2)

file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

# Create subdirectory "bld"
bld_dir <- file.path(file_dir, "bld")
if (!dir.exists(bld_dir)) {
  dir.create(bld_dir, recursive = TRUE)
}

source_dir <- function(path) {
  files <- list.files(
    path = file.path(file_dir, path),
    pattern = "\\.R$",
    full.names = TRUE
  )
  files <- sort(files)
  
  for (f in files) source(f)
}

source_dir("R")

set.seed(9372)
N_sim <- 100
quantiles_VAR <- c(0.025,0.01,0.005)
deviation_VAR <- test1 <- test2 <- matrix(NA, N_sim, 3*length(quantiles_VAR))
deviation_parameters <- matrix(NA, N_sim, 9)

for (i in 1:N_sim) {
  #### Setup Simulation ####
  size_burn_in_sample <- 1000
  size_training_sample <- 1000
  size_eval_sample <- 200
  
  T_grid_points <- 100
  t_difference <- 1/(T_grid_points+1)
  
  #### Simulation ####
  simulation2 <- simulate_fgarch_2(size_burn_in_sample,size_training_sample,
                                   size_eval_sample,T_grid_points,
                                   epsilon_function = epsilon_function_t_ou,
                                   nu = 5
  )
  
  t_grid <- simulation2$grid
  y_matrix <- simulation2$y_train
  y_matrix_squared <- simulation2$y_train^2
  y_eval <- simulation2$y_eval
  y_eval_squared <- y_eval^2
  sigma_squared <- simulation2$sigma_squared
  
  alpha_simulation <- simulation2$alpha
  beta_simulation <- simulation2$beta
  delta_simulation <- simulation2$delta
  
  #### Functional Principal Components ####
  
  fpca_results <- baseline_compute(y_matrix_squared)
  
  PC_matrix_qml <- as.matrix(fpca_results$eigenvectors[,1:2])
  PC_matrix_qml <- cbind(colMeans(y_matrix_squared),PC_matrix_qml)
  PC_matrix_qml <- cbind(rep(1,T_grid_points),PC_matrix_qml)
  
  #### FGARCH estimation ####
  fgarch_estimates <- fgarch_est_qml(y_matrix_squared, PC_matrix_qml)
  
  estimate_delta <- as.vector(fgarch_estimates$delta)
  estimate_alpha <- fgarch_estimates$alpha
  estimate_beta <- fgarch_estimates$beta
  
  #### Fitted values ####
  
  fitted_values <- calculate_fitted_values(y_matrix_squared, y_matrix, PC_matrix_qml,
                                           estimate_delta, estimate_alpha, estimate_beta)
  
  epsilon_fitted <- fitted_values$epsilon_fitted
  sigma_fitted <- fitted_values$sigma_fitted
  alpha_fitted <- fitted_values$alpha
  beta_fitted <- fitted_values$beta
  delta_fitted <- fitted_values$delta
  
  #### Bootstrap ####
  Bootstrap_samples <- 10000
  
  quantile_matrix <- calculate_bootstrap(epsilon_fitted,Bootstrap_samples,quantiles_VAR)
  Forecasts_VAR <- forecast_VAR(y_eval_squared,y_matrix_squared,fitted_values,quantile_matrix)
  
  test_VaR <- compute_test_statistic(Forecasts_VAR$VAR_forecast, quantiles_VAR, y_eval)
  test1[i,1:length(quantiles_VAR)] <- test_VaR$p_value
  deviation_VAR[i,1:length(quantiles_VAR)] <- test_VaR$total_deviation
  
  #### PC-GARCH ####
  
  PCA_pc_garch <- pca_compute(y_matrix)
  PCs_pc_garch <- PCA_pc_garch$eigenvectors[,1]
  PC_garch_results <- pc_garch_estimate(y_matrix, PCs_pc_garch)
  
  estimate_alpha_pc_garch <- PC_garch_results$alpha
  estimate_beta_pc_garch <- PC_garch_results$beta
  estimate_delta_pc_garch <- PC_garch_results$delta
  H_sqrt_fitted_inverse <- PC_garch_results$H_sqrt_fitted
  
  epsilon_fitted_pc_garch <- matrix(NA, nrow = size_training_sample, ncol = T_grid_points)
  for (i in 1:size_training_sample) {
    epsilon_fitted_pc_garch[i,] <- t(H_sqrt_fitted_inverse[i,,]%*%y_matrix[i,])
  }
  
  quantiles_pc_garch <- calculate_bootstrap(epsilon_fitted_pc_garch, Bootstrap_samples, quantiles_VAR)
  
  forecasts_pc_GARCH <- forecast_pc_garch(PC_garch_results$lambda_fitted, y_eval,
                                          y_matrix, PCs_pc_garch, estimate_delta_pc_garch,
                                          estimate_alpha_pc_garch, estimate_beta_pc_garch,
                                          quantiles_pc_garch)
  
  test_VAR_pc_garch <- compute_test_statistic(forecasts_pc_GARCH$y_forecast_VAR, quantiles_VAR, y_eval)
  test1[i,2*length(quantiles_VAR)+1:3*length(quantiles_VAR)] <- test_VaR$p_value
  deviation_VAR[i,2*length(quantiles_VAR)+1:3*length(quantiles_VAR)] <- test_VaR$total_deviation
}

sqrt(sum((alpha_simulation + beta_simulation)^2) * (1/(T_grid_points-1))^2)
sqrt(sum((alpha_fitted+beta_fitted)^2) * (1/(T_grid_points-1))^2)
matplot(
  x = t_grid,
  y = t(sigma_squared[1500:1505,]),
  type = "l",
  lty = 1,
  col = rainbow(5),
  xlab = "t",
  ylab = "f(t)"
)
matplot(
  x = t_grid,
  y = t(sigma_fitted[500:505,]),
  type = "l",
  lty = 1,
  col = rainbow(5),
  xlab = "t",
  ylab = "f(t)"
)

plot(t_grid, quantile_matrix[3,], type = "l", ylim = c(-10,0), col = "#1B4F72")
lines(t_grid, quantile_matrix[2,], col = "#BA4A00")
lines(t_grid, quantile_matrix[1,], col = "#2C2C2C")

plot(t_grid, quantiles_pc_garch[3,], type = "l", ylim = c(-10,0), col = "#1B4F72")
lines(t_grid, quantiles_pc_garch[2,], col = "#BA4A00")
lines(t_grid, quantiles_pc_garch[1,], col = "#2C2C2C")

plot(t_grid, quantiles_pc_garch[3,], type = "l", ylim = c(-10,0), col = "#1B4F72")
lines(t_grid, quantiles_pc_garch[2,], col = "#BA4A00")
lines(t_grid, quantiles_pc_garch[1,], col = "#2C2C2C")

y_forecast_pc_garch_quantile1 <- forecasts_pc_GARCH$y_forecast_VAR[1,,]
plot(t_grid,y_eval[103,], type = "l",
     ylim = c(min(y_forecast_quantile1[103,],
                  y_forecast_pc_garch_quantile1[103,]),
              max(y_eval[103,])), col = "#1B4F72")


lines(t_grid,y_forecast_quantile1[103,], col = "#BA4A00")
lines(t_grid, y_forecast_pc_garch_quantile1[103,], col = "#2C2C2C")