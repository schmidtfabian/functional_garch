# --------------------------------------------------
# Master Thesis — FGARCH Simulation
# Author: Fabian Christian Schmidt
# --------------------------------------------------

rm(list=ls())

install.packages(setdiff(c("rstudioapi","ggplot2", "dplyr", "tidyr", "nloptr",
                           "numDeriv"),
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

N_sim <- 100
quantiles_VAR <- c(0.025,0.01,0.005)
set.seed(9372)
nu_vec = c(5,25)
PC_vec = c(2,4)
test1 <- array(
  NA,
  dim = c(N_sim, length(nu_vec), length(PC_vec), 3, length(quantiles_VAR))
)

deviation_parameters <- array(
  NA,
  dim = c(N_sim, length(nu_vec), length(PC_vec), 2, 3)
)

#### Setup ####
size_burn_in_sample <- 1000
size_training_sample <- 500
size_eval_sample <- 200

T_grid_points <- 100
t_difference <- 1/(T_grid_points-1)

Bootstrap_samples <- 10000
# --------------------------------------------------
# Simple text progress bar
# --------------------------------------------------
total_tasks <- N_sim * length(nu_vec) * length(PC_vec)
pb <- txtProgressBar(min = 0, max = total_tasks, style = 3)
task_counter <- 0

# --------------------------------------------------
# Nested for loop
# --------------------------------------------------
for (q in 1:N_sim) {
  for (i in 1:length(nu_vec)) {
    for (j in 1:length(PC_vec)) {
      
      nu <- nu_vec[i]
      no_PCs <- PC_vec[j]
      
      #### Simulation ####
      simulation2 <- simulate_fgarch_3(
        size_burn_in_sample,
        size_training_sample,
        size_eval_sample,
        T_grid_points,
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
      PC_matrix_ls <- as.matrix(fpca_results$eigenvectors[, 1:(no_PCs)])
      
      baseline_results <- baseline_compute(y_matrix_squared)
      PC_matrix_qml <- as.matrix(baseline_results$eigenvectors[, 1:no_PCs])
      PC_matrix_qml <- cbind(rep(1, T_grid_points), colMeans(y_matrix_squared), PC_matrix_qml)
      
      #### FGARCH estimation ####
      fgarch_estimates_ls <- fgarch_est_ls(y_matrix_squared, PC_matrix_ls,
                                                  positivity_grid_size = 75)
      
      estimate_delta_ls <- as.vector(fgarch_estimates_ls$delta)
      estimate_alpha_ls <- fgarch_estimates_ls$alpha
      estimate_beta_ls <- fgarch_estimates_ls$beta
      
      fgarch_estimates_qml <- fgarch_est_qml(y_matrix_squared, PC_matrix_qml)
      
      estimate_delta_qml <- as.vector(fgarch_estimates_qml$delta)
      estimate_alpha_qml <- fgarch_estimates_qml$alpha
      estimate_beta_qml <- fgarch_estimates_qml$beta
      
      #### Fitted values ####
      fitted_values_ls <- calculate_fitted_values(
        y_matrix_squared, y_matrix, PC_matrix_ls,
        estimate_delta_ls, estimate_alpha_ls, estimate_beta_ls
      )
      
      epsilon_fitted_ls <- fitted_values_ls$epsilon_fitted
      
      fitted_values_qml <- calculate_fitted_values(
        y_matrix_squared, y_matrix, PC_matrix_qml,
        estimate_delta_qml, estimate_alpha_qml, estimate_beta_qml
      )
      
      epsilon_fitted_qml <- fitted_values_qml$epsilon_fitted
      
      #### Bootstrap ####
      quantile_matrix_ls <- calculate_bootstrap(
        epsilon_fitted_ls, Bootstrap_samples, quantiles_VAR
      )
      quantile_matrix_qml <- calculate_bootstrap(
        epsilon_fitted_qml, Bootstrap_samples, quantiles_VAR
      )
      
      #### Forecasting ####
      Forecasts_VAR_ls <- forecast_VAR(
        y_eval_squared, y_matrix_squared, fitted_values_ls, quantile_matrix_ls
      )
      Forecasts_VAR_qml <- forecast_VAR(
        y_eval_squared, y_matrix_squared, fitted_values_qml, quantile_matrix_qml
      )
      
      #### Testing ####
      test_VaR_ls <- compute_test_statistic(
        Forecasts_VAR_ls$VAR_forecast, quantiles_VAR, y_eval,
        delta_sim = simulation2$delta, alpha_sim = simulation2$alpha,
        beta_sim = simulation2$beta, delta_est = fitted_values_ls$delta,
        alpha_est = fitted_values_ls$alpha, beta_est = fitted_values_ls$beta,
        sigma_train = sigma_squared[(size_burn_in_sample+1):(size_burn_in_sample+size_training_sample),],
        sigma_eval = sigma_squared[(size_burn_in_sample+size_training_sample+1):nrow(sigma_squared),],
        sigma_fit = fitted_values_ls$sigma_fitted, sigma_forecast = Forecasts_VAR_ls$sigma_forecast
      )
      
      test_VaR_qml <- compute_test_statistic(
        Forecasts_VAR_qml$VAR_forecast, quantiles_VAR, y_eval,
        delta_sim = simulation2$delta, alpha_sim = simulation2$alpha,
        beta_sim = simulation2$beta, delta_est = fitted_values_qml$delta,
        alpha_est = fitted_values_qml$alpha, beta_est = fitted_values_qml$beta,
        sigma_train = sigma_squared[(size_burn_in_sample+1):(size_burn_in_sample+size_training_sample),],
        sigma_eval = sigma_squared[(size_burn_in_sample+size_training_sample+1):nrow(sigma_squared),],
        sigma_fit = fitted_values_qml$sigma_fitted, sigma_forecast = Forecasts_VAR_qml$sigma_forecast
      )
      
      #### PC-GARCH ####
      PCA_pc_garch <- pca_compute(y_matrix)
      PCs_pc_garch <- PCA_pc_garch$eigenvectors[, 1:(no_PCs)]
      PC_garch_results <- pc_garch_estimate(y_matrix, PCs_pc_garch)
      
      estimate_alpha_pc_garch <- PC_garch_results$alpha
      estimate_beta_pc_garch <- PC_garch_results$beta
      estimate_delta_pc_garch <- PC_garch_results$delta
      H_sqrt_fitted_inverse <- PC_garch_results$H_sqrt_fitted
      
      epsilon_fitted_pc_garch <- matrix(
        NA, nrow = size_training_sample, ncol = T_grid_points
      )
      
      for (k in 1:size_training_sample) {
        epsilon_fitted_pc_garch[k, ] <- t(H_sqrt_fitted_inverse[k, , ] %*% y_matrix[k, ])
      }
      
      quantiles_pc_garch <- calculate_bootstrap(
        epsilon_fitted_pc_garch, Bootstrap_samples, quantiles_VAR
      )
      
      forecasts_pc_GARCH <- forecast_pc_garch(
        PC_garch_results$lambda_fitted, y_eval,
        y_matrix, PCs_pc_garch, estimate_delta_pc_garch,
        estimate_alpha_pc_garch, estimate_beta_pc_garch,
        quantiles_pc_garch
      )
      
      test_VAR_pc_garch <- compute_test_statistic(
        forecasts_pc_GARCH$y_forecast_VAR, quantiles_VAR, y_eval,
        compute_tests_other = FALSE
      )
      
      #### Write results directly into arrays ####
      test1[q, i, j, 1, ] <- test_VaR_ls$p_value
      deviation_parameters[q, i, j, 1, 1] <- test_VaR_ls$msd_delta
      deviation_parameters[q, i, j, 1, 2] <- test_VaR_ls$msd_alpha
      deviation_parameters[q, i, j, 1, 3] <- test_VaR_ls$msd_beta
      
      test1[q, i, j, 2, ] <- test_VaR_qml$p_value
      deviation_parameters[q, i, j, 2, 1] <- test_VaR_qml$msd_delta
      deviation_parameters[q, i, j, 2, 2] <- test_VaR_qml$msd_alpha
      deviation_parameters[q, i, j, 2, 3] <- test_VaR_qml$msd_beta
      
      test1[q, i, j, 3, ] <- test_VAR_pc_garch$p_value
      
      #### Update progress bar ####
      task_counter <- task_counter + 1
      setTxtProgressBar(pb, task_counter)
    }
  }
}

close(pb)

t_grid <- simulation2$grid
alpha_simulation <- simulation2$alpha
beta_simulation <- simulation2$beta
delta_simulation <- simulation2$delta

delta_norm_sim <- sqrt(sum((delta_simulation)^2) * t_difference)
alpha_norm_sim <- (t_difference * RSpectra::svds(alpha_simulation, k = 1, nu = 0, nv = 0)$d[1])^2
beta_norm_sim  <- (t_difference * RSpectra::svds(beta_simulation, k = 1, nu = 0, nv = 0)$d[1])^2

#### Create mean and standard deviation arrays #####
mean_p_value <- apply(test1, c(2, 3, 4, 5), mean, na.rm = TRUE)
sd_p_value   <- apply(test1, c(2, 3, 4, 5), sd,   na.rm = TRUE)

mean_dev_par <- apply(deviation_parameters, c(2, 3, 4, 5), mean, na.rm = TRUE)

true_norms <- c(delta_norm_sim, alpha_norm_sim, beta_norm_sim)

rel_mean_dev_par <- array(NA, dim = dim(mean_dev_par))

for (i in 1:length(true_norms)) {
  rel_mean_dev_par[,,,i] <- sqrt(mean_dev_par[,,,i])/true_norms[i]
}


mean_p_value_ls <- mean_p_value[, , 1, ]
sd_p_value_ls <- sd_p_value[, , 1, ]
mean_p_value_qml <- mean_p_value[, , 2, ]
sd_p_value_qml <- sd_p_value[, , 2, ]
mean_p_value_pc_garch <- mean_p_value[, , 3, ]
sd_p_value_pc_garch <- sd_p_value[, , 3, ]

rel_mean_dev_par_ls <- rel_mean_dev_par[,,1,]
rel_mean_dev_par_qml <- rel_mean_dev_par[,,2,]

#### Plots (to test functionality) ####
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
sigma_forecast_ls <- Forecasts_VAR_ls$sigma_forecast
sigma_forecast_qml <- Forecasts_VAR_qml$sigma_forecast
sqrt(sum((alpha_simulation + beta_simulation)^2) * (1/(T_grid_points-1))^2)
sqrt(sum((alpha_fitted_ls + beta_fitted_ls)^2) * (1/(T_grid_points-1))^2)
sqrt(sum((alpha_fitted_qml + beta_fitted_qml)^2) * (1/(T_grid_points-1))^2)

matplot(
  x = t_grid,
  y = t(sigma_squared[1200:1205,]),
  type = "l",
  lty = 1,
  col = rainbow(5),
  xlab = "t",
  ylab = "f(t)"
)
matplot(
  x = t_grid,
  y = t(sigma_fitted_ls[200:205,]),
  type = "l",
  lty = 1,
  col = rainbow(5),
  xlab = "t",
  ylab = "f(t)"
)
matplot(
  x = t_grid,
  y = t(sigma_fitted_qml[200:205,]),
  type = "l",
  lty = 1,
  col = rainbow(5),
  xlab = "t",
  ylab = "f(t)"
)
matplot(
  x = t_grid,
  y = t(sigma_squared[1600:1605,]),
  type = "l",
  lty = 1,
  col = rainbow(5),
  xlab = "t",
  ylab = "f(t)"
)
matplot(
  x = t_grid,
  y = t(sigma_forecast_ls[100:105,]),
  type = "l",
  lty = 1,
  col = rainbow(5),
  xlab = "t",
  ylab = "f(t)"
)
matplot(
  x = t_grid,
  y = t(sigma_forecast_qml[100:105,]),
  type = "l",
  lty = 1,
  col = rainbow(5),
  xlab = "t",
  ylab = "f(t)"
)

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

epsilon_fitted_simulation <- y_matrix / sqrt(sigma_squared[1001:1500,])
epsilon_fitted_simulation[is.na(epsilon_fitted_simulation)] <- 0
quantile_matrix_simulation <- calculate_bootstrap(
  epsilon_fitted_simulation, Bootstrap_samples, quantiles_VAR
)

plot(t_grid, quantile_matrix_ls[1,], type = "l", ylim = c(-4,-1), col = "#1B4F72")
lines(t_grid, quantile_matrix_qml[1,], col = "#BA4A00")
lines(t_grid, quantile_matrix_simulation[1,], col = "#2C2C2C")

plot(t_grid, quantile_matrix_ls[2,], type = "l", ylim = c(-4,-1), col = "#1B4F72")
lines(t_grid, quantile_matrix_qml[2,], col = "#BA4A00")
lines(t_grid, quantile_matrix_simulation[2,], col = "#2C2C2C")

plot(t_grid, quantile_matrix_ls[3,], type = "l", ylim = c(-5,-2), col = "#1B4F72")
lines(t_grid, quantile_matrix_qml[3,], col = "#BA4A00")
lines(t_grid, quantile_matrix_simulation[3,], col = "#2C2C2C")

#### Plotting ####
# Build data frame
df_1 <- data.frame(
  t = t_grid,
  Actual      = y_eval[5, ],
  VAR_ls      = Forecasts_VAR_ls$VAR_forecast[1, 5, ],
  VAR_qml     = Forecasts_VAR_qml$VAR_forecast[1, 5, ],
  PC_GARCH    = forecasts_pc_GARCH$y_forecast_VAR[1, 5, ]
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



kernel_df_simulation <- expand.grid(
  t = t_grid,
  s = t_grid
)

kernel_df_simulation$alpha <- as.vector(alpha_simulation)
p_simulation <- ggplot(kernel_df_simulation, aes(x = t, y = s)) +
  geom_raster(aes(fill = alpha)) +
  geom_contour(aes(z = alpha), color = "black", alpha = 0.35, linewidth = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  coord_fixed() +
  scale_fill_gradient(
    low = "white",
    high = "grey30",
    name = expression(alpha(t, s))
  ) +
  labs(
    x = "t",
    y = "s",
    title = "Simulation kernel alpha"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

p_simulation

#ggsave("kernel_bw_plot.pdf", plot = p_simulation, width = 6, height = 5)

kernel_df_ls <- expand.grid(
  t = t_grid,
  s = t_grid
)

kernel_df_ls$alpha <- as.vector(alpha_fitted_ls)
p_ls <- ggplot(kernel_df_ls, aes(x = t, y = s)) +
  geom_raster(aes(fill = alpha)) +
  geom_contour(aes(z = alpha), color = "black", alpha = 0.35, linewidth = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  coord_fixed() +
  scale_fill_gradient(
    low = "white",
    high = "grey30",
    name = expression(alpha(t, s))
  ) +
  labs(
    x = "t",
    y = "s",
    title = "Estimated kernel alpha LS"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

p_ls

kernel_df_qml <- expand.grid(
  t = t_grid,
  s = t_grid
)

kernel_df_qml$alpha <- as.vector(alpha_fitted_qml)
p_qml <- ggplot(kernel_df_qml, aes(x = t, y = s)) +
  geom_raster(aes(fill = alpha)) +
  geom_contour(aes(z = alpha), color = "black", alpha = 0.35, linewidth = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  coord_fixed() +
  scale_fill_gradient(
    low = "white",
    high = "grey30",
    name = expression(alpha(t, s))
  ) +
  labs(
    x = "t",
    y = "s",
    title = "Estimated kernel alpha QML"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

p_qml

