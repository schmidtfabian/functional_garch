#### Setup ####
rm(list=ls())
library(fda)

#### Simulation ####

set.seed(9372)
size_burn_in_sample <- 1000
size_training_sample <- 1000
size_eval_sample <- 200
N_observations <- size_burn_in_sample + size_training_sample + size_eval_sample
T_grid_points <- 285
t_grid <- seq(0, 1, length.out = T_grid_points)
t_difference <- 1 / T_grid_points

delta_vector <- rep(0.01, T_grid_points)
alpha_function <- beta_function <- function(s, t) 12*t*(1-t)*s*(1-s)

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
  scale_factor <- sqrt(log(2)) / (2^(200 * t))
  epsilon <- scale_factor * B
  
  return(epsilon)
}




for (i in 1:N_observations) {
  if (i > 1) {
    sigma_squared_matrix[i,] <- delta_vector +
      sapply(t_grid, function(t)
        sum(alpha_function(t, t_grid) * y_matrix[i-1,]^2) * t_difference) +
      sapply(t_grid, function(t)
        sum(beta_function(t, t_grid) * sigma_squared_matrix[i-1,]) * t_difference)
  }
  
  y_matrix[i,] <- sqrt(sigma_squared_matrix[i,]) * epsilon_function(t_grid)
}

#y_matrix_burn_in <- y_matrix[1:size_burn_in_sample,] # remove burn in sample
y_matrix_evaluation <- y_matrix[(size_burn_in_sample + size_training_sample+1):(size_burn_in_sample + size_training_sample + size_eval_sample),]
y_matrix <- y_matrix[(size_burn_in_sample+1):(size_burn_in_sample + size_training_sample),]


#### Principal Components Estimation ####

# Demean y and create squared matrix
#y_matrix <- sweep(y_matrix,2,colMeans(y_matrix), FUN = "-")
y_matrix_squared <- y_matrix^2

# Smooth Curves
K <- 50
BK.basis <- create.bspline.basis(rangeval=c(0,1), nbasis=K)
BK.basis_Parobj <- fdPar(BK.basis, Lfdobj = 2, # states that the second derivative is penalized 
                         lambda = 1e-3)

y_matrix_squared_fd <- smooth.basis(y=t(y_matrix_squared), fdParobj = BK.basis_Parobj, argvals = t_grid)
y_matrix_fd <- smooth.basis(y=t(y_matrix), fdParobj = BK.basis_Parobj, argvals = t_grid)

# Create smoothed y and optimal Principal Component for plotting
y_eval <- t(eval.fd(t_grid, y_matrix_fd$fd))
optimal_pc <- function(t) sqrt(30)*t*(1-t)
optimal_pc_fd <- smooth.basis(argvals = t_grid, y = optimal_pc(t_grid), fdParobj = BK.basis_Parobj)
optimal_pc_fd <- optimal_pc_fd$fd

# Estimate FPC
y_matrix_squared_fd.pca <- pca.fd(y_matrix_squared_fd$fd, nharm = 1, centerfns = TRUE)
y_matrix_squared_fd.pca$varprop

#### Plotting ####
#Plot the first 5 curves
first_indices_estimation <- 1:5

par(mfrow = c(5, 1), mar = c(3, 4, 2, 1))  # layout 5x1

for (i in first_indices_estimation) {
  plot(t_grid, y_matrix[i, ], type = "l", col = "black",
       main = paste("y_", i, "(t)", sep=""), xlab = "t", ylab = "y_i(t)",
       xlim = c(0,1))
  abline(h = 0, col = "gray")
}
# Plot the smooth curves
for (i in first_indices_estimation) {
  plot(t_grid, y_eval[i,], type = "l", col = "black",
       main = paste("y_", i, "(t)", sep=""), xlab = "t", ylab = "y_i(t)",
       xlim = c(0,1))
  abline(h = 0, col = "gray")
}


par(mfrow=c(1,1))

# Comparing optimal PC with first FPC
plot(y_matrix_squared_fd.pca$harmonics, lwd = 3)
lines(optimal_pc_fd)


#### GARCH model estimation ####

# create basis based on the first principal component
PC_estimated <- y_matrix_squared_fd.pca$harmonics[1]

# Project y_matrix_squared into finite dimensional basis
scores_y_squared <- as.vector(inprod(y_matrix_squared_fd$fd, PC_estimated))

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

scores_y_optimal_pc <- as.vector(inprod(y_matrix_squared_fd$fd,optimal_pc_fd))

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
PC_estimated_basis <- PC_estimated$basis
PC_estimated_coefs <- as.vector(PC_estimated$coefs)

alpha_kernel_estimated_PC_fd <- bifd(coef = res1ls$par[2]*(PC_estimated_coefs%o%PC_estimated_coefs), sbasisobj = PC_estimated_basis, tbasisobj = PC_estimated_basis) 
beta_kernel_estimated_PC_fd <- bifd(coef = res1ls$par[3]*(PC_estimated_coefs%o%PC_estimated_coefs), sbasisobj = PC_estimated_basis, tbasisobj = PC_estimated_basis) 

delta_function_estimated_PC_fd <- res2ls$par[1]*PC_estimated

alpha_kernel_matrix_estimated_PC <- eval.bifd(t_grid, t_grid, alpha_kernel_estimated_PC_fd)
beta_kernel_matrix_estimated_PC <- eval.bifd(t_grid, t_grid, beta_kernel_estimated_PC_fd)

delta_vector_estimated_PC <- as.vector(eval.fd(t_grid, delta_function_estimated_PC_fd))

sigma_squared_estimated_PC <- matrix(NA, nrow = dim(y_matrix_squared)[1], ncol = T_grid_points)
mean_y_matrix_squared <- colMeans(y_matrix_squared)
sigma_squared_estimated_PC[1,] <- delta_vector_estimated_PC + as.vector(alpha_kernel_matrix_estimated_PC %*% y_matrix_squared[1,])
for (i in 2:dim(y_matrix_squared)[1]){
  sigma_squared_estimated_PC[i,] <- delta_vector_estimated_PC +
    as.vector(alpha_kernel_matrix_estimated_PC %*% y_matrix_squared[i-1,])*t_difference +
    as.vector(beta_kernel_matrix_estimated_PC %*% sigma_squared_estimated_PC[i-1, ])*t_difference
}



optimal_PC_basis <- optimal_pc_fd$basis
optimal_PC_coefs <- as.vector(optimal_pc_fd$coefs)

alpha_kernel_optimal_PC_fd <- bifd(coef = res1ls$par[2]*(optimal_PC_coefs%o%optimal_PC_coefs), sbasisobj = optimal_PC_basis, tbasisobj = optimal_PC_basis) 
beta_kernel_optimal_PC_fd <- bifd(coef = res1ls$par[3]*(optimal_PC_coefs%o%optimal_PC_coefs), sbasisobj = optimal_PC_basis, tbasisobj = optimal_PC_basis) 

delta_function_optimal_PC_fd <- res2ls$par[1]*optimal_pc_fd

alpha_kernel_matrix_optimal_PC <- eval.bifd(t_grid, t_grid, alpha_kernel_optimal_PC_fd)
beta_kernel_matrix_optimal_PC <- eval.bifd(t_grid, t_grid, beta_kernel_optimal_PC_fd)

delta_vector_optimal_PC <- as.vector(eval.fd(t_grid, delta_function_optimal_PC_fd))

sigma_squared_optimal_PC <- matrix(NA, nrow = dim(y_matrix_squared)[1], ncol = T_grid_points)
sigma_squared_optimal_PC[1,] <- delta_vector_optimal_PC + as.vector(alpha_kernel_matrix_optimal_PC %*% y_matrix_squared[1,])
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
bootstrap_matrix <- matrix(NA, nrow = Bootstrap_samples, ncol = T_grid_points)
for (i in 1:Bootstrap_samples){
  for (j in 1:T_grid_points){
    bootstrap_matrix[i,j] <- sample(epsilon_fitted_optimal_pc[,j], size = 1, replace = FALSE) # I am not sure if this is correct or the entire curve should be resampled.
  }
}

quantiles_VAR <- c(0.025,0.01,0.005)

quantile_matrix <- matrix(NA, nrow = length(quantiles_VAR), ncol = T_grid_points)
for (i in 1:T_grid_points){
  quantile_matrix[,i] <- quantile(bootstrap_matrix[,i], probs = quantiles_VAR)
}

plot(t_grid, quantile_matrix[3,], type = "l", ylim = c(-10,0), col = "blue")
lines(t_grid, quantile_matrix[2,], col = "red")
lines(t_grid, quantile_matrix[1,], col = "yellow")


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

# Estimated PC is negative at t=1
test_eval <- eval.fd(t_grid, y_matrix_squared_fd.pca$harmonics)

# Estimated PC only accounts for 30% of the variance compared to 70% in the paper.

# I am not sure if I understand the Bootstrap passage correctly.