# --------------------------------------------------
# Master Thesis — FGARCH Simulation
# Author: Fabian Christian Schmidt
# --------------------------------------------------

rm(list=ls())

install.packages(setdiff(c("rstudioapi","ggplot2", "dplyr", "tidyr", "nloptr"),
                         rownames(installed.packages())))
library(ggplot2)
library(tidyr)
library(dplyr)

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

N_sim <- 5
quantiles_VAR <- c(0.025,0.01,0.005)
set.seed(9372)
nu_vec = c(5,10)
PC_vec = c(1,3)
deviation_VAR <- test1 <- test2 <- array(NA, dim = c(N_sim, length(nu_vec),
                                                     length(PC_vec),2,length(quantiles_VAR)))
deviation_parameters <- array(NA, dim = c(N_sim, length(nu_vec),
                                          length(PC_vec),2,3))

#### Setup ####
size_burn_in_sample <- 1000
size_training_sample <- 1000
size_eval_sample <- 200

T_grid_points <- 100
t_difference <- 1/(T_grid_points-1)

Bootstrap_samples <- 10000

for (i in 1:N_sim) {
  #### Simulation ####
  simulation2 <- simulate_fgarch_3(size_burn_in_sample,size_training_sample,
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
  
  fpca_results <- fpca_compute(y_matrix_squared)
  
  PC_matrix_ls <- as.matrix(fpca_results$eigenvectors[,1:5])
  
  baseline_results <- baseline_compute(y_matrix_squared)
  PC_matrix_qml <- as.matrix(baseline_results$eigenvectors[,1:3])
  PC_matrix_qml <- cbind(rep(1,T_grid_points),colMeans(y_matrix_squared),PC_matrix_qml)
  
  #### FGARCH estimation ####
  
  fgarch_estimates_ls <- fgarch_est_ls(y_matrix_squared,PC_matrix_ls)
  
  estimate_delta_ls <- as.vector(fgarch_estimates_ls$delta)
  estimate_alpha_ls <- fgarch_estimates_ls$alpha
  estimate_beta_ls <- fgarch_estimates_ls$beta
  
  fgarch_estimates_qml <- fgarch_est_qml(y_matrix_squared, PC_matrix_qml)
  
  estimate_delta_qml <- as.vector(fgarch_estimates_qml$delta)
  estimate_alpha_qml <- fgarch_estimates_qml$alpha
  estimate_beta_qml <- fgarch_estimates_qml$beta
  
  #### Fitted values ####
  
  fitted_values_ls <- calculate_fitted_values(y_matrix_squared, y_matrix, PC_matrix_ls,
                                              estimate_delta_ls, estimate_alpha_ls, estimate_beta_ls)
  
  epsilon_fitted_ls <- fitted_values_ls$epsilon_fitted
  sigma_fitted_ls <- fitted_values_ls$sigma_fitted
  alpha_fitted_ls <- fitted_values_ls$alpha
  beta_fitted_ls <- fitted_values_ls$beta
  delta_fitted_ls <- fitted_values_ls$delta
  
  fitted_values_qml <- calculate_fitted_values(y_matrix_squared, y_matrix, PC_matrix_qml,
                                               estimate_delta_qml, estimate_alpha_qml, estimate_beta_qml)
  
  epsilon_fitted_qml <- fitted_values_qml$epsilon_fitted
  sigma_fitted_qml <- fitted_values_qml$sigma_fitted
  alpha_fitted_qml <- fitted_values_qml$alpha
  beta_fitted_qml <- fitted_values_qml$beta
  delta_fitted_qml <- fitted_values_qml$delta
  
  #### Bootstrap ####
  quantile_matrix_ls <- calculate_bootstrap(epsilon_fitted_ls,Bootstrap_samples,quantiles_VAR)
  quantile_matrix_qml <- calculate_bootstrap(epsilon_fitted_qml,Bootstrap_samples,quantiles_VAR)
  
  #### Forecasting ####
  Forecasts_VAR_ls <- forecast_VAR(y_eval_squared,y_matrix_squared,fitted_values_ls,quantile_matrix_ls)
  Forecasts_VAR_qml <- forecast_VAR(y_eval_squared,y_matrix_squared,fitted_values_qml,quantile_matrix_qml)
  
  #### Testing ####
  test_VaR_ls <- compute_test_statistic(Forecasts_VAR_ls$VAR_forecast, quantiles_VAR, y_eval)
  test1[q,i,j,1,] <- test_VaR_ls$p_value
  deviation_VAR[q,i,j,1,] <- test_VaR_ls$total_deviation
  
  test_VaR_qml <- compute_test_statistic(Forecasts_VAR_qml$VAR_forecast, quantiles_VAR, y_eval)
  test1[q,i,j,2,]<- test_VaR_qml$p_value
  deviation_VAR[q,i,j,2,] <- test_VaR_qml$total_deviation
  
  #### PC-GARCH ####
  
  PCA_pc_garch <- pca_compute(y_matrix)
  PCs_pc_garch <- PCA_pc_garch$eigenvectors[,1:3]
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
sqrt(sum((alpha_fitted_ls + beta_fitted_ls)^2) * (1/(T_grid_points-1))^2)
sqrt(sum((alpha_fitted_qml + beta_fitted_qml)^2) * (1/(T_grid_points-1))^2)

# Build data frame
df_1 <- data.frame(
  t = t_grid,
  Actual      = y_eval[100, ],
  VAR_ls      = Forecasts_VAR_ls$VAR_forecast[1, 100, ],
  VAR_qml     = Forecasts_VAR_qml$VAR_forecast[1, 100, ],
  PC_GARCH    = forecasts_pc_GARCH$y_forecast_VAR[1, 100, ]
)

# Convert to long format
df_1_long <- df_1 %>%
  pivot_longer(
    cols = -t,
    names_to = "Series",
    values_to = "value"
  )

# Plot
ggplot(df_1_long, aes(x = t, y = value,
                    linetype = Series,
                    shape = Series)) +
  geom_line(linewidth = 1, color = "black") +
  scale_linetype_manual(values = c(
    "Actual"   = "solid",
    "VAR_ls"   = "dashed",
    "VAR_qml"  = "dotted",
    "PC_GARCH" = "dotdash"
  )) +
  labs(x = "t", y = "f(t)", linetype = "Series", shape = "Series") +
  theme_bw()


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

persp(t_grid, t_grid, alpha_simulation,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "alpha(s,t)",
      ticktype = "detailed", main = "Real alpha")
persp(t_grid, t_grid, alpha_fitted_ls,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "alpha(s,t)",
      ticktype = "detailed", main = "Estimated alpha ls")
persp(t_grid, t_grid, alpha_fitted_qml,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "alpha(s,t)",
      ticktype = "detailed", main = "Estimated alpha qml")

persp(t_grid, t_grid, beta_simulation,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "beta(s,t)",
      ticktype = "detailed", main = "Real beta")
persp(t_grid, t_grid, beta_fitted_ls,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "beta(s,t)",
      ticktype = "detailed", main = "Estimated beta ls")
persp(t_grid, t_grid, beta_fitted_qml,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "beta(s,t)",
      ticktype = "detailed", main = "Estimated beta qml")