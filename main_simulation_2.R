# --------------------------------------------------
# Master Thesis — FGARCH Simulation
# Author: Fabian Christian Schmidt
# --------------------------------------------------

rm(list=ls())

install.packages(setdiff(c("rstudioapi","ggplot2", "dplyr", "tidyr"),
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
N_sim <- 100
quantiles_VAR <- c(0.025,0.01,0.005)
set.seed(9372)
deviation_VAR <- test1 <- test2 <- matrix(NA, N_sim, 2*length(quantiles_VAR))
deviation_parameters <- matrix(NA, N_sim, 6)
nu = 5

#### Setup ####
size_burn_in_sample <- 1000
size_training_sample <- 1000
size_eval_sample <- 200

T_grid_points <- 100
t_difference <- 1/(T_grid_points+1)

#### Simulation ####
simulation2 <- simulate_fgarch_2(size_burn_in_sample,size_training_sample,
                                 size_eval_sample,T_grid_points,
                                 epsilon_function = epsilon_function_t_ou,
                                 nu = nu
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

PC_matrix_ls <- as.matrix(fpca_results$eigenvectors[,1:2])

baseline_results <- baseline_compute(y_matrix_squared)
PC_matrix_qml <- as.matrix(baseline_results$eigenvectors[,1:1])
PC_matrix_qml <- cbind(colMeans(y_matrix_squared),PC_matrix_qml)
PC_matrix_qml <- cbind(rep(1,T_grid_points),PC_matrix_qml)

#### FGARCH estimation ####

fgarch_estimates_ls <- fgarch_est_ls(y_matrix_squared,PC_matrix_ls)

estimate_delta_ls <- as.vector(fgarch_estimates_ls$delta)
estimate_alpha_ls <- fgarch_estimates_ls$alpha
estimate_beta_ls <- fgarch_estimates_ls$beta

#### Fitted values ####

fitted_values_ls <- calculate_fitted_values(y_matrix_squared, y_matrix, PC_matrix_ls,
                                         estimate_delta_ls, estimate_alpha_ls, estimate_beta_ls)

epsilon_fitted_ls <- fitted_values_ls$epsilon_fitted
sigma_fitted_ls <- fitted_values_ls$sigma_fitted
alpha_fitted_ls <- fitted_values_ls$alpha
beta_fitted_ls <- fitted_values_ls$beta
delta_fitted_ls <- fitted_values_ls$delta

Bootstrap_samples <- 10000

quantile_matrix_ls <- calculate_bootstrap(epsilon_fitted_ls,Bootstrap_samples,quantiles_VAR)

Forecasts_VAR_ls <- forecast_VAR(y_eval_squared,y_matrix_squared,fitted_values_ls,quantile_matrix_ls)

test_VaR <- compute_test_statistic(Forecasts_VAR_ls$VAR_forecast, quantiles_VAR, y_eval)
#test1[i,1:length(quantiles_VAR)] <- test_VaR$p_value
#deviation_VAR[i,1:length(quantiles_VAR)] <- test_VaR$total_deviation


sqrt(sum((alpha_simulation + beta_simulation)^2) * (1/(T_grid_points-1))^2)
sqrt(sum((alpha_fitted_ls+beta_fitted_ls)^2) * (1/(T_grid_points-1))^2)
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
  y = t(sigma_fitted_ls[500:505,]),
  type = "l",
  lty = 1,
  col = rainbow(5),
  xlab = "t",
  ylab = "f(t)"
)

plot(t_grid, quantile_matrix_ls[3,], type = "l", ylim = c(-10,0), col = "#1B4F72")
lines(t_grid, quantile_matrix_ls[2,], col = "#BA4A00")
lines(t_grid, quantile_matrix_ls[1,], col = "#2C2C2C")

persp(t_grid, t_grid, alpha_simulation,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "alpha(s,t)",
      ticktype = "detailed")
persp(t_grid, t_grid, alpha_fitted_ls,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "alpha(s,t)",
      ticktype = "detailed")

persp(t_grid, t_grid, beta_simulation,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "alpha(s,t)",
      ticktype = "detailed")
persp(t_grid, t_grid, beta_fitted_ls,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "alpha(s,t)",
      ticktype = "detailed")

