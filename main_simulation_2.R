# --------------------------------------------------
# Master Thesis — FGARCH Simulation
# Author: Fabian Christian Schmidt
# --------------------------------------------------

rm(list=ls())

install.packages(setdiff(c("rstudioapi","ggplot2", "dplyr", "tidyr", "nloptr",
                           "RSpectra"
                           ),
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
nu_vec = c(5,10,25)
PC_vec = c(1,3,5)
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

for (i in 1:length(nu_vec)) {
  for (j in 1:length(PC_vec)) {
    for (q in 1:N_sim) {
      nu <- nu_vec[i]
      no_PCs <- PC_vec[j]
      
      #### Simulation ####
      simulation2 <- simulate_fgarch_3(size_burn_in_sample,size_training_sample,
                                       size_eval_sample,T_grid_points,
                                       epsilon_function = epsilon_function_t_ou,
                                       nu = nu
      )
      
      y_matrix <- simulation2$y_train
      y_matrix_squared <- simulation2$y_train^2
      y_eval <- simulation2$y_eval
      y_eval_squared <- y_eval^2
      sigma_squared <- simulation2$sigma_squared
      
      #### Functional Principal Components ####
      
      fpca_results <- fpca_compute(y_matrix_squared)
      
      PC_matrix_ls <- as.matrix(fpca_results$eigenvectors[,1:(no_PCs+2)])
      
      baseline_results <- baseline_compute(y_matrix_squared)
      PC_matrix_qml <- as.matrix(baseline_results$eigenvectors[,1:no_PCs])
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
      
      fitted_values_qml <- calculate_fitted_values(y_matrix_squared, y_matrix, PC_matrix_qml,
                                                   estimate_delta_qml, estimate_alpha_qml, estimate_beta_qml)
      
      epsilon_fitted_qml <- fitted_values_qml$epsilon_fitted
      
      #### Bootstrap ####
      quantile_matrix_ls <- calculate_bootstrap(epsilon_fitted_ls,Bootstrap_samples,quantiles_VAR)
      quantile_matrix_qml <- calculate_bootstrap(epsilon_fitted_qml,Bootstrap_samples,quantiles_VAR)
      
      #### Forecasting ####
      Forecasts_VAR_ls <- forecast_VAR(y_eval_squared,y_matrix_squared,fitted_values_ls,quantile_matrix_ls)
      Forecasts_VAR_qml <- forecast_VAR(y_eval_squared,y_matrix_squared,fitted_values_qml,quantile_matrix_qml)
      
      #### Testing ####
      test_VaR_ls <- compute_test_statistic(Forecasts_VAR_ls$VAR_forecast, quantiles_VAR, y_eval,
                                            delta_sim = simulation2$delta, alpha_sim = simulation2$alpha,
                                            beta_sim = simulation2$beta, delta_est = fitted_values_ls$delta,
                                            alpha_est = fitted_values_ls$alpha, beta_est = fitted_values_ls$beta)
      test1[q,i,j,1,] <- test_VaR_ls$p_value
      deviation_VAR[q,i,j,1,] <- test_VaR_ls$total_deviation
      
      test_VaR_qml <- compute_test_statistic(Forecasts_VAR_qml$VAR_forecast, quantiles_VAR, y_eval,
                                             delta_sim = simulation2$delta, alpha_sim = simulation2$alpha,
                                             beta_sim = simulation2$beta, delta_est = fitted_values_qml$delta,
                                             alpha_est = fitted_values_qml$alpha, beta_est = fitted_values_qml$beta)
      test1[q,i,j,2,]<- test_VaR_qml$p_value
      deviation_VAR[q,i,j,2,] <- test_VaR_qml$total_deviation
    }
  }
}

mean_dev <- apply(deviation_VAR, c(2,3,4,5), mean, na.rm = TRUE)
sd_dev   <- apply(deviation_VAR, c(2,3,4,5), sd,   na.rm = TRUE)
mean_p_value <- apply(test1, c(2,3,4,5), mean, na.rm = TRUE)
sd_p_value   <- apply(test1, c(2,3,4,5), sd,   na.rm = TRUE)

mean_dev_ls <- mean_dev[,,1,]
sd_dev_ls <- sd_dev[,,1,]
mean_dev_qml <- mean_dev[,,2,]
sd_dev_qml <- sd_dev[,,2,]

mean_p_value_ls <- mean_p_value[,,1,]
sd_p_value_ls <- sd_p_value[,,1,]
mean_p_value_qml <- mean_p_value[,,2,]
sd_p_value_qml <- sd_p_value[,,2,]

t_grid <- simulation2$grid
alpha_simulation <- simulation2$alpha
beta_simulation <- simulation2$beta
delta_simulation <- simulation2$delta
sigma_fitted_ls <- fitted_values_ls$sigma_fitted
alpha_fitted_ls <- fitted_values_ls$alpha
beta_fitted_ls <- fitted_values_ls$beta
delta_fitted_ls <- fitted_values_ls$delta
sigma_fitted_qml <- fitted_values_qml$sigma_fitted
alpha_fitted_qml <- fitted_values_qml$alpha
beta_fitted_qml <- fitted_values_qml$beta
delta_fitted_qml <- fitted_values_qml$delta
sqrt(sum((alpha_simulation + beta_simulation)^2) * (1/(T_grid_points-1))^2)
sqrt(sum((alpha_fitted_ls+beta_fitted_ls)^2) * (1/(T_grid_points-1))^2)
sqrt(sum((alpha_fitted_qml+beta_fitted_qml)^2) * (1/(T_grid_points-1))^2)
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
matplot(
  x = t_grid,
  y = t(sigma_fitted_qml[500:505,]),
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


plot(t_grid,delta_simulation, type = "l")
plot(t_grid, delta_fitted_ls, type = "l")
plot(t_grid, delta_fitted_qml, type = "l")


unc_sigma_sim <- solve(diag(T_grid_points) -
                     alpha_simulation *t_difference -
                     beta_simulation *t_difference)%*% delta_simulation
unc_sigma_ls <- solve(diag(T_grid_points) -
                        alpha_fitted_ls *t_difference -
                        beta_fitted_ls *t_difference)%*% delta_fitted_ls

unc_sigma_qml <- solve(diag(T_grid_points) -
                        alpha_fitted_qml *t_difference -
                        beta_fitted_qml *t_difference)%*% delta_fitted_qml
plot(t_grid, unc_sigma, type = "l")
plot(t_grid, unc_sigma_ls, type = "l")
plot(t_grid, unc_sigma_qml, type = "l")


unc_sigma_ls_proj <- solve(diag(nrow(estimate_alpha_ls)) -
                             estimate_alpha_ls -
                             estimate_beta_ls) %*% estimate_delta_ls

plot(1:nrow(estimate_alpha_ls), unc_sigma_ls_proj)

