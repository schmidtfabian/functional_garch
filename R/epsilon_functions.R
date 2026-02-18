# --------------------------------------------------
# Innovation function (Cerovecki et al. 2019)
# --------------------------------------------------
epsilon_function_ou <- function(t) {
  
  B_times <- exp(t)
  dt <- diff(c(0, B_times))
  
  increments <- rnorm(length(dt), mean = 0, sd = sqrt(dt))
  B <- cumsum(increments)
  
  scale_factor <- exp(-t / 2)
  epsilon <- scale_factor * B
  
  return(epsilon)
}

epsilon_function_t_ou <- function(t, nu = 5){
  # Non-Gaussian OU-type innovation using Student-t increments
  # nu > 2 required for finite variance
  
  # Step 1: Compute corresponding Brownian times (same transformation)
  B_times <- exp(t)
  
  # Step 2: Determine time increments
  dt <- diff(c(0, B_times))
  
  # Step 3: Generate non-Gaussian increments with mean 0 and var = dt
  # Standardize Student-t to variance 1, then scale by sqrt(dt)
  increments <- rt(length(dt), df = nu) *
    sqrt((nu - 2) / nu) *
    sqrt(dt)
  
  B <- cumsum(increments)
  
  # Step 4: Construct epsilon_i(t) as in OU case
  scale_factor <- exp(-t/2)
  epsilon <- scale_factor * B
  
  return(epsilon)
}
epsilon_function_laplace_ou <- function(t){
  # OU-type innovation with Laplace increments
  
  B_times <- exp(t)
  dt <- diff(c(0, B_times))
  
  # Laplace(0, b) has variance 2b^2 → choose b = 1/sqrt(2) for variance 1
  u <- runif(length(dt), -0.5, 0.5)
  laplace_std <- -sign(u) * log(1 - 2 * abs(u)) / sqrt(2)
  
  increments <- laplace_std * sqrt(dt)
  B <- cumsum(increments)
  
  scale_factor <- exp(-t/2)
  epsilon <- scale_factor * B
  
  return(epsilon)
}