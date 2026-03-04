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

set.seed(9372)


#### Setup ####
size_burn_in_sample <- 1000
size_training_sample <- 1000
size_eval_sample <- 200

T_grid_points <- 100
t_difference <- 1/(T_grid_points+1)

#### Simulation ####
simulation1 <- simulate_fgarch_1(size_burn_in_sample,size_training_sample,
                                 size_eval_sample,T_grid_points,
                                 epsilon_function = epsilon_function_ou
                                 )

t_grid <- simulation1$grid
y_matrix <- simulation1$y_train
y_matrix_squared <- simulation1$y_train^2
y_eval <- simulation1$y_eval
y_eval_squared <- y_eval^2

#### Functional Principal Components ####

fpca_results <- fpca_compute(y_matrix_squared)

PC_matrix <- as.matrix(fpca_results$eigenvectors[,1:3])
PC_matrix <- cbind(colMeans(y_matrix_squared),PC_matrix)
PC_matrix <- cbind(rep(1,T_grid_points),PC_matrix)

Gram_matrix <- t(PC_matrix[,3:ncol(PC_matrix)])%*%(PC_matrix[,3:ncol(PC_matrix)] * (1/(T_grid_points-1)))
det(Gram_matrix)

optimal_pc_function <- function(t) sqrt(30)*t*(1-t)
optimal_pc <- optimal_pc_function(t_grid)

caption1 <- paste0(
  "Explained variance (estimated PC): ",
  round(fpca_results$explained_variance[1] * 100, 1), "%"
)

df_fpca_comparison <- data.frame(
  t = t_grid,
  True = optimal_pc,
  Estimated = (-PC_matrix[, 3] + max(PC_matrix[,3]))/sqrt(t_difference)
)

ggplot(df_fpca_comparison, aes(x = t)) +
  geom_line(aes(y = Estimated, linetype = "Estimated FPC"), linewidth = 1) +
  geom_line(aes(y = True, linetype = "True PC"), linewidth = 1) +
  labs(
    x = "t",
    y = "Principal component value",
    title = "Estimated vs. true principal component",
    linetype = NULL,
    caption = caption1
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, size = 10),
    plot.caption.position = "plot",
    legend.position = "top"
  )

output_file <- file.path(bld_dir, "pc_comparison.pdf")
ggsave(
  filename = output_file,
  plot = last_plot(),
  width = 6,
  height = 4
)

#norms <- sqrt(colSums(PC_matrix^2) * t_difference)
#PC_matrix_normalized <- sweep(PC_matrix, 2, norms, "/")

#### FGARCH estimation ####

fgarch_estimates <- garch_estimate_scores(y_matrix_squared,PC_matrix, method = "qml")

estimate_delta <- as.vector(fgarch_estimates$delta)
estimate_alpha <- fgarch_estimates$alpha
estimate_beta <- fgarch_estimates$beta

#### Fitted values ####

fitted_values <- calculate_fitted_values(y_matrix_squared, y_matrix, PC_matrix,
                                         estimate_delta, estimate_alpha, estimate_beta)

epsilon_fitted <- fitted_values$epsilon_fitted
sigma_fitted <- fitted_values$sigma_fitted
alpha_fitted <- fitted_values$alpha
beta_fitted <- fitted_values$beta
delta_fitted <- fitted_values$delta

sqrt(sum((alpha_fitted+beta_fitted)^2) * (1/(T_grid_points-1))^2)

quantiles_VAR <- c(0.025,0.01,0.005)
Bootstrap_samples <- 10000

quantile_matrix <- calculate_bootstrap(epsilon_fitted,Bootstrap_samples,quantiles_VAR)

plot(t_grid, quantile_matrix[3,], type = "l", ylim = c(-10,0), col = "#1B4F72")
lines(t_grid, quantile_matrix[2,], col = "#BA4A00")
lines(t_grid, quantile_matrix[1,], col = "#2C2C2C")

Forecasts_VAR <- forecast_VAR(y_eval_squared,y_matrix_squared,fitted_values,quantile_matrix)

y_forecast_quantile1 <- Forecasts_VAR$VAR_forecast[1,,]
plot(t_grid,y_eval[103,], type = "l", ylim = c(-0.6,0.6), col = "#1B4F72")
lines(t_grid,y_forecast_quantile1[103,], col = "#BA4A00")

indicator_matrix <- y_eval <= y_forecast_quantile1
mean(indicator_matrix)

#### PC-GARCH ####

pca_results <- pca_compute(y_matrix)

#### Testing ####

delta_test <- rep(0,T_grid_points)
alpha_test <- matrix(0, nrow = T_grid_points, ncol = T_grid_points)
beta_test <- matrix(0, nrow = T_grid_points, ncol = T_grid_points)
for (i in 1:1) {
  delta_test <- delta_test+fgarch_estimates$delta[i]*optimal_pc
  for (j in 1:1) {
    alpha_test <- alpha_test + estimate_alpha[i,j]*optimal_pc%o%optimal_pc
    beta_test <- beta_test + estimate_beta[i,j]*optimal_pc%o%optimal_pc
  }
}

