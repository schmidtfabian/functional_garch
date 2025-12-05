#### Setup ####
rm(list=ls())

#### Simulation ####

set.seed(9372)
size_burn_in_sample <- 1000
size_training_sample <- 1000
size_eval_sample <- 200
N_observations <- size_burn_in_sample + size_training_sample + size_eval_sample
T_grid_points <- 285
t_grid <- seq(0, 1, length.out = T_grid_points)
t_difference <- 1 / (T_grid_points-1)

delta_vector <- rep(0.01, T_grid_points)
alpha_function <- beta_function <- function(t, s) 12*t*(1-t)*s*(1-s)

alpha_kernel_matrix <- outer(t_grid, t_grid, alpha_function)
beta_kernel_matrix <- outer(t_grid, t_grid, beta_function)

sigma_squared_matrix <- matrix(0, N_observations, T_grid_points)
y_matrix <- matrix(0, N_observations, T_grid_points)
sigma_squared_matrix[1,] <- delta_vector

epsilon_function <- function(t) {
  # The epsilon function Aue et al. use is time shifted. Therefore, we first have to shift t
  # to the right B_times. Then the increments of each B(t_i) to B(t_j) is normally distributed
  # with the standard deviation being t_j-t_i
  # Step 1: Compute corresponding Brownian times
  B_times <- (2^(400 * t)) / log(2)
  
  # Step 2: Determine time increments for Brownian motion
  dt <- diff(c(0, B_times))
  
  # Step 3: Generate Brownian motions for each i
  increments <- rnorm(length(dt), mean = 0, sd = sqrt(dt))
  B <- cumsum(increments)
  
  # Step 4: Compute epsilon_i(t)
  scale_factor <- (sqrt(log(2))) / (2^(200 * t))
  epsilon <- scale_factor * B
  
  return(epsilon)
}

epsilon_function_test <- function(t_grid,
                                    lambda = 200 * log(2),
                                    sigma  = 20 * sqrt(log(2)),
                                    X0 = 0) {
  # t_grid: numeric vector of time points (must be sorted and unique)
  # lambda: mean-reversion speed
  # sigma:  diffusion parameter
  # X0:     initial value at t_grid[1]
  
  # Ensure the grid is sorted
  t_grid <- sort(t_grid)
  n <- length(t_grid)
  
  # Pre-allocate output
  X <- numeric(n)
  X[1] <- X0
  
  # Simulate increments
  for (i in 2:n) {
    dt <- t_grid[i] - t_grid[i-1]        # variable time step
    dW <- sqrt(dt) * rnorm(1)            # Brownian increment
    X[i] <- X[i-1] + (-lambda * X[i-1]) * dt + sigma * dW
  }
  
  return(X)
}



for (i in 1:N_observations) {
  if (i > 1) {
    sigma_squared_matrix[i,] <- delta_vector +
      as.vector(alpha_kernel_matrix %*% y_matrix[i-1,]^2) * t_difference +
      as.vector(beta_kernel_matrix %*% sigma_squared_matrix[i-1,]) * t_difference
  }
  
  y_matrix[i,] <- sqrt(sigma_squared_matrix[i,]) * epsilon_function_test(t_grid)
}

#y_matrix_burn_in <- y_matrix[1:size_burn_in_sample,] # remove burn in sample
y_matrix_evaluation <- y_matrix[(size_burn_in_sample + size_training_sample+1):(size_burn_in_sample + size_training_sample + size_eval_sample),]
y_matrix <- y_matrix[(size_burn_in_sample+1):(size_burn_in_sample + size_training_sample),]


#### Principal Components Estimation ####


y_matrix_squared <- y_matrix^2
# Demean y^2 and create covariance matrix
y_matrix_squared_mean_function <- colMeans(y_matrix_squared)
y_matrix_squared_demeaned <- matrix(NA, nrow = size_training_sample, ncol = T_grid_points)
for (i in 1:size_training_sample) {
  y_matrix_squared_demeaned[i,] <- y_matrix_squared[i,] - y_matrix_squared_mean_function 
}
y_squared_covariance_matrix <- t(y_matrix_squared_demeaned)%*%y_matrix_squared_demeaned*(1/size_training_sample)*t_difference
y_squared_covariance_matrix.eigen <- eigen(y_squared_covariance_matrix)
estimated_first_PC <- y_squared_covariance_matrix.eigen$vectors[,1]
variance_explained_first_PC <- y_squared_covariance_matrix.eigen$values[1]/(sum(y_squared_covariance_matrix.eigen$values))
variance_explained_first_PC

# eigen() returns PCs with Euclidean norm = 1 (sum(v^2) = 1). 
# But functional PCs must satisfy the continuous L2 norm: sum(v^2 * dt) = 1.
# Since dt = 1/(T-1), the eigenvector is too small by factor 1/sqrt(dt).
# Therefore we rescale: estimated_first_PC <- estimated_first_PC / sqrt(t_difference) 
# (and optionally flip its sign to match the theoretical direction).
if (mean(estimated_first_PC)< 0) {
  estimated_first_PC <- -1 * estimated_first_PC * 1/sqrt(t_difference)
} else {
  estimated_first_PC <- estimated_first_PC * 1/sqrt(t_difference)
}



optimal_pc_function <- function(t) sqrt(30)*t*(1-t)
optimal_pc <- optimal_pc_function(t_grid)

#### Plotting ####
#Plot the first 5 curves
first_indices_estimation <- 1:5

layout(mat = matrix(data = c(1,2,3,4,5,0), byrow = TRUE, nrow = 3))

for (i in first_indices_estimation) {
  plot(t_grid, y_matrix[i, ], type = "l", col = "black",
       main = paste("y_", i, "(t)", sep=""), xlab = "t", ylab = "y_i(t)",
       xlim = c(0,1), ylim = c(-0.75,0.75))
  abline(h = 0, col = "gray")
}


for (i in first_indices_estimation) {
  plot(t_grid, y_matrix_squared[i, ], type = "l", col = "black",
       main = paste("y_", i, "(t)^2", sep=""), xlab = "t", ylab = "y_i(t)^2",
       xlim = c(0,1))
  abline(h = 0, col = "gray")
}
layout(1)

# Comparing optimal PC with first FPC
plot(t_grid,optimal_pc, type = "l", ylim = c(0,1.75))
lines(t_grid, estimated_first_PC)


#### GARCH model estimation ####

# Project y_matrix_squared into finite dimensional basis
scores_y_squared <- as.vector(y_matrix_squared%*% estimated_first_PC*t_difference)

# Function to compute the conditional variances given parameters
garch_variance <- function(params, y) {
  delta <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  
  n <- length(y)
  sigma2 <- rep(0,n)
  sigma2[1] <- y[1]
  
  for (t in 2:n) {
    sigma2[t] <- delta + alpha * y[t-1] + beta * sigma2[t-1]
  }
  return(sigma2)
}

# Objective function: least squares
garch_ls_obj <- function(params, y) {
  sigma2 <- garch_variance(params, y)
  crossprod(y - sigma2, y - sigma2)
}

# Objective function: Quasi-Maximum-Likelihood
garch_QML_obj <- function(params, y) {
  sigma2 <- garch_variance(params, y)
  n <- length(y)
  l_t <- rep(0,n)
  for (i in 1:n) {
    l_t[i] <- (y[i]/sigma2[i])+log(sigma2[i])
  }
  1/n*sum(l_t)
}

# Initial values (must ensure positivity constraints)
init_params <- c(delta = 0.1, alpha = 0.5, beta = 0.5)
lower_bound_params <- c(0,0,0)
upper_bound_params <- c(1,1,1)

# Estimate using nonlinear minimization
res1ls <- optim(
  par = init_params,
  fn = garch_ls_obj,
  y = scores_y_squared,
  method = "L-BFGS-B",
  lower = lower_bound_params,
  upper = upper_bound_params
)

res1ls$par  # Estimated parameters (omega, alpha, beta)

res1QML <- optim(
  par = init_params,
  fn = garch_QML_obj,
  y = scores_y_squared,
  method = "L-BFGS-B",
  lower = lower_bound_params,
  upper = upper_bound_params
)

res1QML$par


### Same estimation now with perfect PC

scores_y_optimal_pc <- as.vector(y_matrix_squared %*% optimal_pc*t_difference)

# Estimate using nonlinear minimization
res2ls <- optim(
  par = init_params,
  fn = garch_ls_obj,
  y = scores_y_optimal_pc,
  method = "L-BFGS-B",
  lower = lower_bound_params,
  upper = upper_bound_params
)

res2ls$par  # Estimated parameters (omega, alpha, beta)

res2QML <- optim(
  par = init_params,
  fn = garch_QML_obj,
  y = scores_y_optimal_pc,
  method = "L-BFGS-B",
  lower = lower_bound_params,
  upper = upper_bound_params
)

res2QML$par


#### Estimating fitted values ####

alpha_kernel_matrix_estimated_PC <- res1ls$par[2]*estimated_first_PC %o% estimated_first_PC
beta_kernel_matrix_estimated_PC <- res1ls$par[3]*estimated_first_PC %o% estimated_first_PC

delta_vector_estimated_PC <- res1ls$par[1]*estimated_first_PC

sigma_squared_estimated_PC <- matrix(NA, nrow = dim(y_matrix_squared)[1], ncol = T_grid_points)
sigma_squared_estimated_PC[1,] <- delta_vector_estimated_PC + as.vector(alpha_kernel_matrix_estimated_PC %*% y_matrix_squared_mean_function)*t_difference
for (i in 2:dim(y_matrix_squared)[1]){
  sigma_squared_estimated_PC[i,] <- delta_vector_estimated_PC +
    as.vector(alpha_kernel_matrix_estimated_PC %*% y_matrix_squared[i-1,])*t_difference +
    as.vector(beta_kernel_matrix_estimated_PC %*% sigma_squared_estimated_PC[i-1, ])*t_difference
}

epsilon_fitted_estimated_PC <- y_matrix/sqrt(sigma_squared_estimated_PC)


alpha_kernel_matrix_optimal_PC <- res2ls$par[2]*optimal_pc %o% optimal_pc
beta_kernel_matrix_optimal_PC <- res2ls$par[3]*optimal_pc %o% optimal_pc

delta_vector_optimal_PC <- res2ls$par[1]*optimal_pc

sigma_squared_optimal_PC <- matrix(NA, nrow = dim(y_matrix_squared)[1], ncol = T_grid_points)
sigma_squared_optimal_PC[1,] <- delta_vector_optimal_PC + as.vector(alpha_kernel_matrix_optimal_PC %*% y_matrix_squared_mean_function)*t_difference
for (i in 2:dim(y_matrix_squared)[1]){
  sigma_squared_optimal_PC[i,] <- delta_vector_optimal_PC +
    as.vector(alpha_kernel_matrix_optimal_PC %*% y_matrix_squared[i-1,])*t_difference +
    as.vector(beta_kernel_matrix_optimal_PC %*% sigma_squared_optimal_PC[i-1, ])*t_difference
}

epsilon_fitted_optimal_pc <- y_matrix/sqrt(sigma_squared_optimal_PC)

#### Plotting Fitted Values ####

par(mfrow = c(5, 1), mar = c(3, 4, 2, 1))  # layout 5x1

last_indices_estimation <- (dim(y_matrix_squared)[1]-4):dim(y_matrix_squared)[1]

for (i in first_indices_estimation) {
  plot(t_grid, sigma_squared_matrix[(size_burn_in_sample +i), ], type = "l", col = "black",
       main = paste("real sigma_", i, "(t)^2", sep=""), xlab = "t", ylab = "sigma_i(t)^2",
       xlim = c(0,1))
  abline(h = 0, col = "gray")
}
# Plot the estimated curves
for (i in first_indices_estimation) {
  plot(t_grid, sigma_squared_estimated_PC[i,], type = "l", col = "black",
       main = paste("estimated sigma_", i, "(t)^2", sep=""), xlab = "t", ylab = "sigma_i(t)^2",
       xlim = c(0,1))
  abline(h = 0, col = "gray")
}
# Plot the estimated curves
for (i in first_indices_estimation) {
  plot(t_grid, sigma_squared_optimal_PC[i,], type = "l", col = "black",
       main = paste("estimated sigma_", i, "(t)^2", sep=""), xlab = "t", ylab = "sigma_i(t)^2",
       xlim = c(0,1))
  abline(h = 0, col = "gray")
}


par(mfrow=c(1,1))


test_error <- sigma_squared_matrix[1001:2000,]-sigma_squared_estimated_PC
mean(test_error)

plot(t_grid, test_error[500,])
plot(t_grid, 0.01-delta_vector_estimated_PC)


#### Bootstrap ####


Bootstrap_samples <- 10000
bootstrap_matrix_estimated_PC <- matrix(NA, nrow = Bootstrap_samples, ncol = T_grid_points)
for (i in 1:Bootstrap_samples){
  bootstrap_matrix_estimated_PC[i,] <- sample(1:N_observations, size = 1, replace = FALSE)
}

quantiles_VAR <- c(0.025,0.01,0.005)

quantile_matrix_estimated_PC <- matrix(NA, nrow = length(quantiles_VAR), ncol = T_grid_points)
for (i in 1:T_grid_points){
  quantile_matrix_estimated_PC[,i] <- quantile(bootstrap_matrix_estimated_PC[,i], probs = quantiles_VAR)
}

plot(t_grid, quantile_matrix_estimated_PC[3,], type = "l", ylim = c(-10,0), col = "blue")
lines(t_grid, quantile_matrix_estimated_PC[2,], col = "red")
lines(t_grid, quantile_matrix_estimated_PC[1,], col = "yellow")


bootstrap_matrix_optimal_PC <- matrix(NA, nrow = Bootstrap_samples, ncol = T_grid_points)
for (i in 1:Bootstrap_samples){
  bootstrap_matrix_optimal_PC[i,] <- sample(epsilon_fitted_optimal_pc, size = 1, replace = FALSE)
}

quantile_matrix_optimal_PC <- matrix(NA, nrow = length(quantiles_VAR), ncol = T_grid_points)
for (i in 1:T_grid_points){
  quantile_matrix_optimal_PC[,i] <- quantile(bootstrap_matrix_optimal_PC[,i], probs = quantiles_VAR)
}

plot(t_grid, quantile_matrix_optimal_PC[3,], type = "l", ylim = c(-10,0), col = "blue")
lines(t_grid, quantile_matrix_optimal_PC[2,], col = "red")
lines(t_grid, quantile_matrix_optimal_PC[1,], col = "yellow")


#### Forecasting Value-at-risk ####

y_matrix_evaluation_squared <- y_matrix_evaluation^2
sigma_squared_optimal_PC_eval <- matrix(data = NA, nrow = size_eval_sample, ncol = T_grid_points)
sigma_squared_optimal_PC_eval[1,] <- delta_vector_optimal_PC +
  as.vector(alpha_kernel_matrix_optimal_PC %*% y_matrix_squared[size_training_sample,]) +
  as.vector(beta_kernel_matrix_optimal_PC %*% sigma_squared_optimal_PC[size_training_sample,])
for (i in 2:size_eval_sample){
  sigma_squared_optimal_PC_eval[i,] <- delta_vector_optimal_PC +
    as.vector(alpha_kernel_matrix_optimal_PC %*% y_matrix_evaluation_squared[i-1,])*t_difference +
    as.vector(beta_kernel_matrix_optimal_PC %*% sigma_squared_optimal_PC_eval[i-1, ])*t_difference
}

y_matrix_VAR_forecast_quantile_1 <- matrix(data = NA, nrow = size_eval_sample, ncol = T_grid_points)
y_matrix_VAR_forecast_quantile_1 <- quantile_matrix[1,]*sqrt(sigma_squared_optimal_PC_eval)

indicator_matrix_real_VAR <- y_matrix_evaluation > y_matrix_VAR_forecast_quantile_1
total_deviation <- quantiles_VAR[1] - (1 - mean(as.numeric(indicator_matrix_real_VAR), na.rm = TRUE))
1 - mean(as.numeric(indicator_matrix_real_VAR), na.rm = TRUE)
total_deviation

#### Plotting ####

par(mfrow = c(5, 1), mar = c(3, 4, 2, 1))  # layout 5x1

for (i in first_indices_estimation) {
  plot(t_grid, y_matrix_VAR_forecast_quantile_1[5+i,], type = "l", col = "blue",
       main = paste("VAR forecast vs. real y_", i, "(t)", sep=""), xlab = "t", ylab = "y_i(t)",
       xlim = c(0,1))
  abline(h = 0, col = "gray")
  lines(t_grid, y_matrix_evaluation[i,])
}

par(mfrow=c(1,1))


#### Remaining issues ####

# Estimated PC is negative

# Estimated PC only accounts for little of the variance compared to 70% in the paper.

# I am not sure if I understand the Bootstrap passage correctly.
# draw many realizations of epsilon(t) from both functions and summarise
n_samp <- 200  # how many independent epsilon paths to draw
eps_mat <- matrix(NA, nrow = n_samp, ncol = T_grid_points)

for (j in 1:n_samp) {
  eps_mat[j, ] <- epsilon_function_test(t_grid)
}
matplot(t_grid, t(eps_mat[1:6,]), type = "l", lty = 1, main = "6 sample paths: epsilon_function()")