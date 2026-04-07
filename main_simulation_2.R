# --------------------------------------------------
# Master Thesis — FGARCH Simulation
# Author: Fabian Christian Schmidt
# --------------------------------------------------

rm(list=ls())

install.packages(setdiff(c("rstudioapi","ggplot2", "dplyr", "tidyr", "nloptr",
                           "numDeriv", "RSpectra"),
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
nu_vec = c(3,5,25)
PC_vec = c(2,4)
test1 <- array(
  NA,
  dim = c(N_sim, length(nu_vec), length(PC_vec), 3, length(quantiles_VAR))
)

deviation_parameters <- array(
  NA,
  dim = c(N_sim, length(nu_vec), length(PC_vec), 2, 3)
)
RISE_array <- array(
  NA,
  dim = c(N_sim, length(nu_vec), length(PC_vec), 2, 2)
)
rel_L2_array <- array(
  NA,
  dim = c(N_sim, length(nu_vec), length(PC_vec), 2, 2)
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
      simulation2 <- simulate_fgarch_2(
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
        sigma_train = y_matrix_squared, # set either to sigma_train or y_train
        sigma_eval = y_eval_squared, # set either to sigma_eval or y_eval
        sigma_fit = fitted_values_ls$sigma_fitted, sigma_forecast = Forecasts_VAR_ls$sigma_forecast
      )
      
      test_VaR_qml <- compute_test_statistic(
        Forecasts_VAR_qml$VAR_forecast, quantiles_VAR, y_eval,
        delta_sim = simulation2$delta, alpha_sim = simulation2$alpha,
        beta_sim = simulation2$beta, delta_est = fitted_values_qml$delta,
        alpha_est = fitted_values_qml$alpha, beta_est = fitted_values_qml$beta,
        sigma_train = y_matrix_squared,
        sigma_eval = y_eval_squared,
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
      RISE_array[q, i, j, 1, 1] <- test_VaR_ls$mean_RISE_fit
      RISE_array[q, i, j, 1, 2] <- test_VaR_ls$mean_RISE_eval
      rel_L2_array[q, i, j, 1, 1] <- test_VaR_ls$mean_rel_L2_fit
      rel_L2_array[q, i, j, 1, 2] <- test_VaR_ls$mean_rel_L2_eval
      
      test1[q, i, j, 2, ] <- test_VaR_qml$p_value
      deviation_parameters[q, i, j, 2, 1] <- test_VaR_qml$msd_delta
      deviation_parameters[q, i, j, 2, 2] <- test_VaR_qml$msd_alpha
      deviation_parameters[q, i, j, 2, 3] <- test_VaR_qml$msd_beta
      RISE_array[q, i, j, 2, 1] <- test_VaR_qml$mean_RISE_fit
      RISE_array[q, i, j, 2, 2] <- test_VaR_qml$mean_RISE_eval
      rel_L2_array[q, i, j, 2, 1] <- test_VaR_qml$mean_rel_L2_fit
      rel_L2_array[q, i, j, 2, 2] <- test_VaR_qml$mean_rel_L2_eval
      
      
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

mean_rel_L2 <- apply(rel_L2_array, c(2, 3, 4, 5), mean, na.rm = TRUE)
sd_rel_L2 <- apply(rel_L2_array, c(2, 3, 4, 5), sd, na.rm = TRUE)

mean_RISE <- apply(RISE_array, c(2, 3, 4, 5), mean, na.rm = TRUE)
sd_RISE <- apply(RISE_array, c(2, 3, 4, 5), sd, na.rm = TRUE)

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

mean_RISE_ls <- mean_RISE[,,1,]
sd_RISE_ls <- sd_RISE[,,1,]
mean_rel_L2_ls <- mean_rel_L2[,,1,]
sd_rel_L2_ls <- sd_rel_L2[,,1,]

mean_RISE_qml <- mean_RISE[,,2,]
sd_RISE_qml <- sd_RISE[,,2,]
mean_rel_L2_qml <- mean_rel_L2[,,2,]
sd_rel_L2_qml <- sd_rel_L2[,,2,]

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

Gram_matrix_ls <- t(PC_matrix_ls)%*%(PC_matrix_ls * (1/(T_grid_points-1)))
det(Gram_matrix_ls)
Gram_matrix_qml <- t(PC_matrix_qml)%*%(PC_matrix_qml * (1/(T_grid_points-1)))
det(Gram_matrix_qml)


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
pdf(file.path(bld_dir, "real_alpha.pdf"), width = 7, height = 7)
persp(t_grid, t_grid, alpha_simulation,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "alpha(s,t)",
      ticktype = "detailed", main = "Real Alpha")
dev.off()

pdf(file.path(bld_dir, "estimated_alpha_ls.pdf"), width = 7, height = 7)
persp(t_grid, t_grid, alpha_fitted_ls,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "alpha(s,t)",
      ticktype = "detailed", main = "Estimated Alpha LS")
dev.off()

pdf(file.path(bld_dir, "estimated_alpha_qml.pdf"), width = 7, height = 7)
persp(t_grid, t_grid, alpha_fitted_qml,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "alpha(s,t)",
      ticktype = "detailed", main = "Estimated Alpha QML")
dev.off()

pdf(file.path(bld_dir, "real_beta.pdf"), width = 7, height = 7)
persp(t_grid, t_grid, beta_simulation,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "beta(s,t)",
      ticktype = "detailed", main = "Real Beta")
dev.off()

pdf(file.path(bld_dir, "estimated_beta_ls.pdf"), width = 7, height = 7)
persp(t_grid, t_grid, beta_fitted_ls,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "beta(s,t)",
      ticktype = "detailed", main = "Estimated Beta LS")
dev.off()

pdf(file.path(bld_dir, "estimated_beta_qml.pdf"), width = 7, height = 7)
persp(t_grid, t_grid, beta_fitted_qml,
      theta = 30, phi = 30, expand = 0.5,
      xlab = "s", ylab = "t", zlab = "beta(s,t)",
      ticktype = "detailed", main = "Estimated Beta QML")
dev.off()

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
  Actual      = y_eval[50, ],
  FGARCH_LS   = Forecasts_VAR_ls$VAR_forecast[1, 50, ],
  FGARCH_QML  = Forecasts_VAR_qml$VAR_forecast[1, 50, ],
  PC_GARCH    = forecasts_pc_GARCH$y_forecast_VAR[1, 50, ]
)

# Convert to long format
df_1_long <- df_1 %>%
  pivot_longer(
    cols = -t,
    names_to = "Series",
    values_to = "value"
  )

# Plot
comparison_forecasts <- ggplot(df_1_long, aes(x = t, y = value,
                    linetype = Series,
                    shape = Series)) +
  geom_line(linewidth = 1, color = "black") +
  scale_linetype_manual(values = c(
    "Actual"   = "solid",
    "FGARCH_LS"   = "dashed",
    "FGARCH_QML"  = "dotted",
    "PC_GARCH" = "dotdash"
  )) +
  labs(x = "t", y = "f(t)", linetype = "Series", shape = "Series") +
  theme_bw()

ggsave(file.path(bld_dir,"comparison_forecasts.pdf"), plot = comparison_forecasts, width = 6, height = 5)

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


#### Create p-value tables ####
dimnames(mean_p_value_ls) <- list(
  nu = nu_vec,
  no_PCs = PC_vec,
  parameter = c("q_1", "q_2", "q_3")
)
dimnames(sd_p_value_ls) <- list(
  nu = nu_vec,
  no_PCs = PC_vec,
  parameter = c("q_1", "q_2", "q_3")
)
noquote(format(round(mean_p_value_ls, 4), nsmall = 4, scientific = FALSE))
noquote(format(round(sd_p_value_ls, 4), nsmall = 4, scientific = FALSE))

dimnames(mean_p_value_qml) <- list(
  nu = nu_vec,
  no_PCs = PC_vec,
  parameter = c("q_1", "q_2", "q_3")
)
dimnames(sd_p_value_qml) <- list(
  nu = nu_vec,
  no_PCs = PC_vec,
  parameter = c("q_1", "q_2", "q_3")
)
noquote(format(round(mean_p_value_qml, 4), nsmall = 4, scientific = FALSE))
noquote(format(round(sd_p_value_qml, 4), nsmall = 4, scientific = FALSE))

dimnames(mean_p_value_pc_garch) <- list(
  nu = nu_vec,
  no_PCs = PC_vec,
  parameter = c("q_1", "q_2", "q_3")
)
dimnames(sd_p_value_pc_garch) <- list(
  nu = nu_vec,
  no_PCs = PC_vec,
  parameter = c("q_1", "q_2", "q_3")
)
noquote(format(round(mean_p_value_pc_garch, 4), nsmall = 4, scientific = FALSE))
noquote(format(round(sd_p_value_pc_garch, 4), nsmall = 4, scientific = FALSE))

par_names <- c("q_1", "q_2", "q_3")

# formatter: 4 decimals, no scientific notation
fmt <- function(x) {
  out <- format(round(x, 4), nsmall = 4, scientific = FALSE, trim = TRUE)
  out <- sub("^(-?)0\\.", "\\1.", out)   # .1234 instead of 0.1234
  out
}

cell_fun <- function(mu, sd) {
  paste0("\\makecell[c]{", fmt(mu), " \\\\[-2pt] (", fmt(sd), ")}")
}

# build rows
rows <- character(length(PC_vec))

for (j in seq_along(PC_vec)) {   # sample size
  vals <- c()
  
  for (i in seq_along(nu_vec)) {                # nu
    for (k in seq_along(par_names)) {           # parameter
      vals <- c(vals, cell_fun(mean_p_value_ls[i, j, k], sd_p_value_ls[i, j, k]))
    }
  }
  
  rows[j] <- paste0(
    PC_vec[j], " & ",
    paste(vals, collapse = " & "),
    " \\\\"
  )
}

latex_tab <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{p-value for LS estimation by number of basis functions used, parameter $\nu$, and quantile}",
  "\\label{TablepvalueLS2}",
  "\\begin{tabular}{c ccc ccc ccc}",
  "\\toprule",
  "\\multirow{2}{*}{\\makecell[c]{M}} & \\multicolumn{3}{c}{$\\nu = 3$} & \\multicolumn{3}{c}{$\\nu = 5$} & \\multicolumn{3}{c}{$\\nu = 25$} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}",
  "& $q_{0.025}$ & $q_{0.01}$ & $q_{0.005}$ & $q_{0.025}$ & $q_{0.01}$ & $q_{0.005}$ & $q_{0.025}$ & $q_{0.01}$ & $q_{0.005}$ \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

#cat(paste(latex_tab, collapse = "\n"))
writeLines(latex_tab, file.path(bld_dir, "mean_p_value_ls_table_2.tex"))


# build rows
rows <- character(length(PC_vec))

for (j in seq_along(PC_vec)) {   # sample size
  vals <- c()
  
  for (i in seq_along(nu_vec)) {                # nu
    for (k in seq_along(par_names)) {           # parameter
      vals <- c(vals, cell_fun(mean_p_value_qml[i, j, k], sd_p_value_qml[i, j, k]))
    }
  }
  
  rows[j] <- paste0(
    PC_vec[j], " & ",
    paste(vals, collapse = " & "),
    " \\\\"
  )
}

latex_tab <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Mean p-value and standard deviation for QML}",
  "\\label{TablepvalueQML2}",
  "\\begin{tabular}{c ccc ccc ccc}",
  "\\toprule",
  "\\multirow{2}{*}{\\makecell[c]{M}} & \\multicolumn{3}{c}{$\\nu = 3$} & \\multicolumn{3}{c}{$\\nu = 5$} & \\multicolumn{3}{c}{$\\nu = 25$} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}",
  "& $q_{0.025}$ & $q_{0.01}$ & $q_{0.005}$ & $q_{0.025}$ & $q_{0.01}$ & $q_{0.005}$ & $q_{0.025}$ & $q_{0.01}$ & $q_{0.005}$ \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

#cat(paste(latex_tab, collapse = "\n"))
writeLines(latex_tab, file.path(bld_dir, "mean_p_value_qml_table_2.tex"))

# build rows
rows <- character(length(PC_vec))

for (j in seq_along(PC_vec)) {   # sample size
  vals <- c()
  
  for (i in seq_along(nu_vec)) {                # nu
    for (k in seq_along(par_names)) {           # parameter
      vals <- c(vals, cell_fun(mean_p_value_pc_garch[i, j, k], sd_p_value_pc_garch[i, j, k]))
    }
  }
  
  rows[j] <- paste0(
    PC_vec[j], " & ",
    paste(vals, collapse = " & "),
    " \\\\"
  )
}

latex_tab <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Mean p-value and standard deviation for PC-GARCH}",
  "\\label{TablepvaluePCGARCH}",
  "\\begin{tabular}{c ccc ccc ccc}",
  "\\toprule",
  "\\multirow{2}{*}{\\makecell[c]{M}} & \\multicolumn{3}{c}{$\\nu = 3$} & \\multicolumn{3}{c}{$\\nu = 5$} & \\multicolumn{3}{c}{$\\nu = 25$} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}",
  "& $q_{0.025}$ & $q_{0.01}$ & $q_{0.005}$ & $q_{0.025}$ & $q_{0.01}$ & $q_{0.005}$ & $q_{0.025}$ & $q_{0.01}$ & $q_{0.005}$ \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

#cat(paste(latex_tab, collapse = "\n"))
writeLines(latex_tab, file.path(bld_dir, "mean_p_value_pc_garch_table.tex"))

#### Tables of relative squared mean deviation of the parameters ####
dimnames(rel_mean_dev_par_ls) <- list(
  nu = nu_vec,
  no_PCs = PC_vec,
  parameter = c("$\\hat{\\delta}$", "$\\hat{\\alpha}$", "$\\hat{\\beta}$")
)

noquote(format(round(rel_mean_dev_par_ls, 4), nsmall = 4, scientific = FALSE))

# build rows
rows <- character(length(PC_vec))

for (j in seq_along(PC_vec)) {   # no_PCs
  vals <- c()
  
  for (i in seq_along(nu_vec)) {          # nu
    for (k in seq_along(par_names)) {     # parameter
      vals <- c(vals, fmt(rel_mean_dev_par_ls[i, j, k]))
    }
  }
  
  rows[j] <- paste0(
    PC_vec[j], " & ",
    paste(vals, collapse = " & "),
    " \\\\"
  )
}

latex_tab <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Relative mean squared deviation for LS estimation by number of basis functions used and parameter $\\nu$}",
  "\\label{TableRMSDLS}",
  "\\begin{tabular}{c ccc ccc ccc}",
  "\\toprule",
  "\\multirow{2}{*}{\\makecell[c]{M}} & \\multicolumn{3}{c}{$\\nu = 3$} & \\multicolumn{3}{c}{$\\nu = 5$} & \\multicolumn{3}{c}{$\\nu = 25$} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}",
  "& $\\hat{\\delta}$ & $\\hat{\\alpha}$ & $\\hat{\\beta}$ & $\\hat{\\delta}$ & $\\hat{\\alpha}$ & $\\hat{\\beta}$ & $\\hat{\\delta}$ & $\\hat{\\alpha}$ & $\\hat{\\beta}$ \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(latex_tab, file.path(bld_dir, "rel_mean_dev_par_ls.tex"))

dimnames(rel_mean_dev_par_qml) <- list(
  nu = nu_vec,
  no_PCs = PC_vec,
  parameter = c("$\\hat{\\delta}$", "$\\hat{\\alpha}$", "$\\hat{\\beta}$")
)

noquote(format(round(rel_mean_dev_par_qml, 4), nsmall = 4, scientific = FALSE))

# build rows
rows <- character(length(PC_vec))

for (j in seq_along(PC_vec)) {   # no_PCs
  vals <- c()
  
  for (i in seq_along(nu_vec)) {          # nu
    for (k in seq_along(par_names)) {     # parameter
      vals <- c(vals, fmt(rel_mean_dev_par_qml[i, j, k]))
    }
  }
  
  rows[j] <- paste0(
    PC_vec[j], " & ",
    paste(vals, collapse = " & "),
    " \\\\"
  )
}

latex_tab <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Relative mean squared deviation for QML estimation by number of basis functions used and parameter $\\nu$}",
  "\\label{TableRMSDQML}",
  "\\begin{tabular}{c ccc ccc ccc}",
  "\\toprule",
  "\\multirow{2}{*}{\\makecell[c]{M}} & \\multicolumn{3}{c}{$\\nu = 3$} & \\multicolumn{3}{c}{$\\nu = 5$} & \\multicolumn{3}{c}{$\\nu = 25$} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}",
  "& $\\hat{\\delta}$ & $\\hat{\\alpha}$ & $\\hat{\\beta}$ & $\\hat{\\delta}$ & $\\hat{\\alpha}$ & $\\hat{\\beta}$ & $\\hat{\\delta}$ & $\\hat{\\alpha}$ & $\\hat{\\beta}$ \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(latex_tab, file.path(bld_dir, "rel_mean_dev_par_qml.tex"))