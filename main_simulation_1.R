# --------------------------------------------------
# Master Thesis — FGARCH Simulation
# Author: Fabian Christian Schmidt
# --------------------------------------------------

rm(list = ls())

install.packages(
  setdiff(
    c("rstudioapi", "ggplot2", "dplyr", "tidyr", "nloptr",
      "future.apply", "progressr"),
    rownames(installed.packages())
  )
)

library(ggplot2)
library(future.apply)
library(progressr)

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
quantiles_VAR <- c(0.025, 0.01, 0.005)
sizes_training_sample <- c(250,500,1000)
set.seed(9372)

nu_vec <- c(5,10,25)

deviation_VAR <- test1 <- array(
  NA_real_,
  dim = c(N_sim, length(nu_vec), length(sizes_training_sample), 2, length(quantiles_VAR))
)

deviation_parameters <- array(
  NA_real_,
  dim = c(N_sim, length(nu_vec),length(sizes_training_sample), 2, 3)
)

#### Setup ####
size_burn_in_sample <- 1000
size_eval_sample <- 200

T_grid_points <- 100
t_difference <- 1 / (T_grid_points - 1)

Bootstrap_samples <- 10000

# --------------------------------------------------
# Parallel setup
# --------------------------------------------------

# Use all but one core
future::plan(future::multisession, workers = max(1, future::availableCores() - 1))

# Simple text progress bar
progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

# Build task grid: one row = one former loop iteration
task_grid <- expand.grid(
  q = 1:N_sim,
  i = 1:length(nu_vec),
  j = 1:length(sizes_training_sample)
)

# --------------------------------------------------
# Function for one simulation task
# --------------------------------------------------
run_one_task <- function(q, i, j) {
  nu <- nu_vec[i]
  size_training_sample <- sizes_training_sample[j]
  
  #### Simulation ####
  simulation1 <- simulate_fgarch_1(
    size_burn_in_sample,
    size_training_sample,
    size_eval_sample,
    T_grid_points,
    epsilon_function = epsilon_function_t_ou,
    nu = nu
  )
  
  y_matrix <- simulation1$y_train
  y_matrix_squared <- simulation1$y_train^2
  y_eval <- simulation1$y_eval
  y_eval_squared <- y_eval^2
  sigma_squared <- simulation1$sigma_squared
  t_grid <- simulation1$grid
  
  
  #### Functional Principal Components ####
  # Assume knowledge of perfect principal component
  perfect_PC <- sqrt(30) * t_grid * (1 - t_grid)
  
  #### FGARCH estimation ####
  fgarch_estimates_ls <- fgarch_est_ls(y_matrix_squared, perfect_PC)
  
  estimate_delta_ls <- as.vector(fgarch_estimates_ls$delta)
  estimate_alpha_ls <- fgarch_estimates_ls$alpha
  estimate_beta_ls <- fgarch_estimates_ls$beta
  
  estimates_ls <- c(estimate_delta_ls[1], estimate_alpha_ls[1,1], estimate_beta_ls[1,1])
  
  fgarch_estimates_qml <- fgarch_est_qml(y_matrix_squared, perfect_PC)
  
  estimate_delta_qml <- as.vector(fgarch_estimates_qml$delta)
  estimate_alpha_qml <- fgarch_estimates_qml$alpha
  estimate_beta_qml <- fgarch_estimates_qml$beta
  
  estimates_qml <- c(estimate_delta_qml[1], estimate_alpha_qml[1,1], estimate_beta_qml[1,1])
  
  #### Fitted values ####
  fitted_values_ls <- calculate_fitted_values(
    y_matrix_squared, y_matrix, perfect_PC,
    estimate_delta_ls, estimate_alpha_ls, estimate_beta_ls
  )
  
  epsilon_fitted_ls <- fitted_values_ls$epsilon_fitted
  
  fitted_values_qml <- calculate_fitted_values(
    y_matrix_squared, y_matrix, perfect_PC,
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
    Forecasts_VAR_ls$VAR_forecast, quantiles_VAR, y_eval, delta_sim = NULL,
    delta_est = NULL, alpha_sim = NULL, alpha_est = NULL, beta_sim = NULL,
    beta_est = NULL, compute_msd = FALSE
  )
  
  test_VaR_qml <- compute_test_statistic(
    Forecasts_VAR_qml$VAR_forecast, quantiles_VAR, y_eval, delta_sim = NULL,
    delta_est = NULL, alpha_sim = NULL, alpha_est = NULL, beta_sim = NULL,
    beta_est = NULL, compute_msd = FALSE
  )
  
  list(
    q = q,
    i = i,
    j = j,
    p_value_ls = test_VaR_ls$p_value,
    deviation_ls = test_VaR_ls$total_deviation,
    p_value_qml = test_VaR_qml$p_value,
    deviation_qml = test_VaR_qml$total_deviation,
    estimates_ls = estimates_ls,
    estimates_qml = estimates_qml
  )
}

future::plan(future::sequential)
# --------------------------------------------------
# Run in parallel with progress bar
# --------------------------------------------------
with_progress({
  p <- progressr::progressor(along = 1:nrow(task_grid))
  
  results <- future.apply::future_lapply(
    X = 1:nrow(task_grid),
    FUN = function(idx) {
      p(sprintf("Task %d / %d", idx, nrow(task_grid)))
      
      row <- task_grid[idx, ]
      run_one_task(q = row$q, i = row$i, j = row$j)
    },
    future.seed = TRUE
  )
})

# --------------------------------------------------
# Write results back into arrays
# --------------------------------------------------
for (res in results) {
  q <- res$q
  i <- res$i
  j <- res$j
  
  test1[q, i, j, 1, ] <- res$p_value_ls
  deviation_VAR[q, i, j, 1, ] <- res$deviation_ls
  deviation_parameters[q, i, j, 1, ] <- res$estimates_ls
  
  test1[q, i, j, 2, ] <- res$p_value_qml
  deviation_VAR[q, i, j, 2, ] <- res$deviation_qml
  deviation_parameters[q, i, j, 2, ] <- res$estimates_qml
}

# Optional: return to sequential plan afterwards
future::plan(future::sequential)

mean_dev <- apply(deviation_VAR, c(2, 3, 4, 5), mean, na.rm = TRUE)
sd_dev   <- apply(deviation_VAR, c(2, 3, 4, 5), sd,   na.rm = TRUE)
mean_p_value <- apply(test1, c(2, 3, 4, 5), mean, na.rm = TRUE)
sd_p_value   <- apply(test1, c(2, 3, 4, 5), sd,   na.rm = TRUE)
mean_dev_par <- apply(deviation_parameters, c(2, 3, 4, 5), mean, na.rm = TRUE)
sd_dev_par   <- apply(deviation_parameters, c(2, 3, 4, 5), sd,   na.rm = TRUE)

mean_dev_ls <- mean_dev[, , 1, ]
sd_dev_ls <- sd_dev[, , 1, ]
mean_dev_qml <- mean_dev[, , 2, ]
sd_dev_qml <- sd_dev[, , 2, ]
mean_dev_par_ls <- mean_dev_par[, , 1, ]
sd_dev_par_ls <- sd_dev_par[, , 1, ]

mean_p_value_ls <- mean_p_value[, , 1, ]
sd_p_value_ls <- sd_p_value[, , 1, ]
mean_p_value_qml <- mean_p_value[, , 2, ]
sd_p_value_qml <- sd_p_value[, , 2, ]
mean_dev_par_qml <- mean_dev_par[, , 2, ]
sd_dev_par_qml <- sd_dev_par[, , 2, ]

true_par <- c(0.009, 0.4, 0.4)

mean_bias_ls <- sweep(mean_dev_par_ls,3,true_par, FUN = "-")
dimnames(mean_bias_ls) <- list(
  nu = nu_vec,
  sample_size = sizes_training_sample,
  parameter = c("d", "a", "b")
)
dimnames(sd_dev_par_ls) <- list(
  nu = nu_vec,
  sample_size = sizes_training_sample,
  parameter = c("d", "a", "b")
)
noquote(format(round(mean_bias_ls, 4), nsmall = 4, scientific = FALSE))
noquote(format(round(sd_dev_par_ls, 4), nsmall = 4, scientific = FALSE))

mean_bias_qml <- sweep(mean_dev_par_qml,3,true_par, FUN = "-")
dimnames(mean_bias_qml) <- list(
  nu = nu_vec,
  sample_size = sizes_training_sample,
  parameter = c("d", "a", "b")
)
dimnames(sd_dev_par_qml) <- list(
  nu = nu_vec,
  sample_size = sizes_training_sample,
  parameter = c("d", "a", "b")
)
noquote(format(round(mean_bias_qml, 4), nsmall = 4, scientific = FALSE))
noquote(format(round(sd_dev_par_qml, 4), nsmall = 4, scientific = FALSE))

par_names <- c("d", "a", "b")

# formatter: 4 decimals, no scientific notation
fmt <- function(x) {
  out <- format(round(x, 4), nsmall = 4, scientific = FALSE, trim = TRUE)
  out <- sub("^(-?)0\\.", "\\1.", out)   # .1234 instead of 0.1234
  out
}

#cell_fun <- function(mu, sd) {
#  paste0("\\makecell[c]{", fmt(mu), " \\\\ (", fmt(sd), ")}")
#}
cell_fun <- function(mu, sd) {
  paste0("\\makecell[c]{", fmt(mu), " \\\\[-2pt] (", fmt(sd), ")}")
}

# build rows
rows <- character(length(sizes_training_sample))

for (j in seq_along(sizes_training_sample)) {   # sample size
  vals <- c()
  
  for (i in seq_along(nu_vec)) {                # nu
    for (k in seq_along(par_names)) {           # parameter
      vals <- c(vals, cell_fun(mean_bias_ls[i, j, k], sd_dev_par_ls[i, j, k]))
    }
  }
  
  rows[j] <- paste0(
    sizes_training_sample[j], " & ",
    paste(vals, collapse = " & "),
    " \\\\"
  )
}

latex_tab <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Mean bias and standard deviation of estimators for LS}",
  "\\begin{tabular}{c ccc ccc ccc}",
  "\\toprule",
  "\\multirow{2}{*}{Sample size} & \\multicolumn{3}{c}{$\\nu = 5$} & \\multicolumn{3}{c}{$\\nu = 10$} & \\multicolumn{3}{c}{$\\nu = 25$} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}",
  "& $d$ & $a$ & $b$ & $d$ & $a$ & $b$ & $d$ & $a$ & $b$ \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

#cat(paste(latex_tab, collapse = "\n"))
writeLines(latex_tab, file.path(bld_dir, "mean_bias_sd_ls_table.tex"))