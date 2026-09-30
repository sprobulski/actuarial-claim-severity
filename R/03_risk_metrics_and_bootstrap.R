#### VaR and ES for the selected model - Frechet-Pareto ####
p_level <- 0.995

# Value-at-Risk (VaR) Calculation
calc_VaR <- function(p, p_func) {
  tryCatch({
    result <- uniroot(function(x) p_func(x) - p, lower = 0, upper = 1e10, tol = 1e-8, extendInt = "yes")
    return(result$root)
  }, error = function(e) return(NA))
}

VaR_995 <- calc_VaR(p_level, p_spl2)
cat("VaR(0.995) =", round(VaR_995, 2), "\n")

# Expected Shortfall (ES) Calculation
calc_ES <- function(p, p_func) {
  tryCatch({
    integrand <- function(u) sapply(u, function(ui) calc_VaR(ui, p_func))
    integral_val <- integrate(integrand, lower = p, upper = 1 - 1e-8, subdivisions = 1000)$value
    ES <- (1 / (1 - p)) * integral_val
    return(ES)
  }, error = function(e) return(NA))
}

ES_995 <- calc_ES(p_level, p_spl2)
cat("ES(0.995) =", round(ES_995, 2), "\n")

#### Parametric and Non-parametric Bootstrap for VaR and ES ####
set.seed(123) 
B <- 100 # Number of iterations 
n_obs <- length(x)
threshold <- 40000

# Result vectors
VaR_np <- numeric(B)
ES_np <- numeric(B)
VaR_p <- numeric(B)
ES_p <- numeric(B)

# Function generating CDF for the model
make_p_spl2_boot <- function(model_boot) {
  function(t) {
    sapply(t, function(v) {
      if(v <= model_boot$u) {
        F_u = exp(-(model_boot$theta_b / model_boot$u)^model_boot$alpha_b)
        model_boot$c_weight * exp(-(model_boot$theta_b / v)^model_boot$alpha_b) / F_u
      } else {
        model_boot$c_weight + (1 - model_boot$c_weight) * (1 - (model_boot$lambda_t / (v - model_boot$u + model_boot$lambda_t))^model_boot$alpha_t)
      }
    })
  }
}

for(i in 1:B) {
  # NON-PARAMETRIC BOOTSTRAP
  # Sampling from historical data (with replacement)
  x_boot_np <- sample(x, n_obs, replace = TRUE)
  
  # Model re-estimation
  model_boot_np <- fit_spliced_frechet_pareto(x_boot_np, u = threshold)
  p_spl2_np <- make_p_spl2_boot(model_boot_np)
  
  # Metrics calculation
  VaR_np[i] <- calc_VaR(p_level, p_spl2_np)
  ES_np[i]  <- calc_ES(p_level, p_spl2_np)
  
  # PARAMETRIC BOOTSTRAP
  u_sim <- runif(n_obs)
  x_boot_p <- sapply(u_sim, function(ui) calc_VaR(ui, p_spl2))
  
  # Model re-estimation
  model_boot_p <- fit_spliced_frechet_pareto(x_boot_p, u = threshold)
  p_spl2_p <- make_p_spl2_boot(model_boot_p)
  
  # Metrics calculation
  VaR_p[i] <- calc_VaR(p_level, p_spl2_p)
  ES_p[i]  <- calc_ES(p_level, p_spl2_p)
  
  if(i %% 10 == 0) cat("Completed:", i, "/", B, "iterations\n")
}

# Confidence intervals calculation
calc_t_bootstrap_CI <- function(theta_est, boot_estimates, alpha_ci = 0.05) {
  sigma_n <- sd(boot_estimates, na.rm = TRUE)
  T_n <- (boot_estimates - theta_est) / sigma_n
  
  H_inv_lower <- quantile(T_n, alpha_ci / 2, na.rm = TRUE)
  H_inv_upper <- quantile(T_n, 1 - alpha_ci / 2, na.rm = TRUE)
  
  CI_lower <- theta_est - sigma_n * H_inv_upper
  CI_upper <- theta_est - sigma_n * H_inv_lower
  names(CI_lower) <- NULL
  names(CI_upper) <- NULL
  
  return(data.frame(Estimate = round(theta_est, 2), Lower_Bound = round(CI_lower, 2), Upper_Bound = round(CI_upper, 2)))
}

cat("\n--- CONFIDENCE INTERVAL RESULTS (VaR) ---\n")
cat("Non-parametric:\n")
print(calc_t_bootstrap_CI(VaR_995, VaR_np))
cat("\nParametric:\n")
print(calc_t_bootstrap_CI(VaR_995, VaR_p))

cat("\n--- CONFIDENCE INTERVAL RESULTS (ES) ---\n")
cat("Non-parametric:\n")
print(calc_t_bootstrap_CI(ES_995, ES_np))
cat("\nParametric:\n")
print(calc_t_bootstrap_CI(ES_995, ES_p))

### VaR and ES based on GPD approximation
p_level <- 0.995
u_gpd <- 40000  
xi <- 0.65      

# Scale parameter calculation
mean_excess_emp <- mean(x[x > u_gpd] - u_gpd)
beta <- mean_excess_emp * (1 - xi)

# Empirical CDF value at point u
F_u <- ecdf(x)(u_gpd)

# VaR calculation
VaR_gpd <- u_gpd + (beta / xi) * (((1 - p_level) / (1 - F_u))^(-xi) - 1)
cat("VaR GPD (0.995) =", round(VaR_gpd, 2), "\n")

# ES calculation
ES_gpd <- (VaR_gpd / (1 - xi)) + ((beta - xi * u_gpd) / (1 - xi))
cat("ES GPD (0.995) =", round(ES_gpd, 2), "\n")