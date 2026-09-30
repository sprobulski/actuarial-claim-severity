#### Distribution Fitting ####

# Otwarcie zapisu do zbiorczego pliku PDF
pdf("plots/02_model_diagnostics.pdf", width = 12, height = 8)

# Pareto Distribution
model_pareto = fit_pareto(x, "ML", "Complete", 0)
print(summary_fit(model_pareto, x))
distribution_fit(model_pareto, x)
actuarial_fit(model_pareto, x)
statistical_fit(model_pareto, x)

# Gamma Distribution
model_gamma = fit_gamma(x, "ML", "Complete", 0)
exp_est = model_gamma$alpha/model_gamma$beta
gradient_h = c(1/model_gamma$beta, -model_gamma$alpha/model_gamma$beta^2)
var_exp_est = t(gradient_h) %*% model_gamma$cov %*% gradient_h

print(summary_fit(model_gamma, x))
distribution_fit(model_gamma, x)
actuarial_fit(model_gamma, x)
statistical_fit(model_gamma, x)

# Frechet Distribution
model_frechet = fit_frechet(x)
print(summary_fit(model_frechet, x))
distribution_fit(model_frechet, x)
actuarial_fit(model_frechet, x)
statistical_fit(model_frechet, x)

# Loglogistic Distribution
model_loglogistic = fit_loglogistic(x)
print(summary_fit(model_loglogistic, x))
distribution_fit(model_loglogistic, x)
actuarial_fit(model_loglogistic, x)
statistical_fit(model_loglogistic, x)

## Spliced Distributions

# Model 1 - Frechet (body), Loglogistic (tail)
model_frechet_loglogistic = fit_spliced_frechet_loglogis(x, u = 40000)
print(summary_fit(model_frechet_loglogistic, x))
distribution_fit(model_frechet_loglogistic, x)
actuarial_fit(model_frechet_loglogistic, x)
statistical_fit(model_frechet_loglogistic, x)

# Model 2 - Frechet (body), Pareto (tail)
model_frechet_pareto = fit_spliced_frechet_pareto(x, u = 40000)
print(summary_fit(model_frechet_pareto, x))
distribution_fit(model_frechet_pareto, x)
actuarial_fit(model_frechet_pareto, x)
statistical_fit(model_frechet_pareto, x)

## Expectation-Maximization Method
model_em_exp = fit_em_exp(x, u = 40000)
print(summary_fit(model_em_exp, x))
distribution_fit(model_em_exp, x)
actuarial_fit(model_em_exp, x)
statistical_fit(model_em_exp, x)

#### GOODNESS-OF-FIT ASSESSMENT ####

# 1. PARETO 
cat("\n--- Diagnostics: PARETO ---\n")
p_par <- function(t) 1 - (model_pareto$lambda / (t + model_pareto$lambda))^model_pareto$alpha
d_par <- function(t) (model_pareto$alpha * model_pareto$lambda^model_pareto$alpha) / ((t + model_pareto$lambda)^(model_pareto$alpha + 1))
plot_pit(x, p_par, "Pareto")
print(paste("A-D test:", calc_ad(x, p_par)$decision))
print(calc_aic(x, d_par, k = 2))
plot_exposure(x, p_par, lev_func = NULL, model_name = "Pareto")

# 2. GAMMA 
cat("\n--- Diagnostics: GAMMA ---\n")
p_gam <- function(t) pgamma(t, shape = model_gamma$alpha, rate = model_gamma$beta)
d_gam <- function(t) dgamma(t, shape = model_gamma$alpha, rate = model_gamma$beta)
plot_pit(x, p_gam, "Gamma")
print(paste("A-D test:", calc_ad(x, p_gam)$decision))
print(calc_aic(x, d_gam, k = 2))
plot_exposure(x, p_gam, lev_func = NULL, model_name = "Gamma")

# 3. FRECHET 
cat("\n--- Diagnostics: FRECHET ---\n")
p_fre <- function(t) exp(-(model_frechet$theta / t)^model_frechet$alpha)
d_fre <- function(t) exp(log(model_frechet$alpha) + model_frechet$alpha * log(model_frechet$theta) - (model_frechet$alpha + 1) * log(t) - (model_frechet$theta / t)^model_frechet$alpha)
plot_pit(x, p_fre, "Frechet")
print(paste("A-D test:", calc_ad(x, p_fre)$decision))
print(calc_aic(x, d_fre, k = 2))
plot_exposure(x, p_fre, lev_func = NULL, model_name = "Frechet")

# 4. LOGLOGISTIC
cat("\n--- Diagnostics: LOGLOGISTIC ---\n")
p_llog <- function(t) {
  u_val = (t / model_loglogistic$theta)^model_loglogistic$gamma
  return(u_val / (1 + u_val))
}
d_llog <- function(t) exp(log(model_loglogistic$gamma) + model_loglogistic$gamma * (log(t) - log(model_loglogistic$theta)) - log(t) - 2 * log(1 + (t / model_loglogistic$theta)^model_loglogistic$gamma))
plot_pit(x, p_llog, "Loglogistic")
print(paste("A-D test:", calc_ad(x, p_llog)$decision))
print(calc_aic(x, d_llog, k = 2))
plot_exposure(x, p_llog, lev_func = NULL, model_name = "Loglogistic")

# 5. Spliced Frechet-Loglogistic
cat("\n--- Diagnostics: SPLICED (Frechet & Loglogistic) ---\n")
p_spl1 <- function(t) {
  sapply(t, function(v) {
    if(v <= model_frechet_loglogistic$u) {
      F_u = exp(-(model_frechet_loglogistic$theta_b / model_frechet_loglogistic$u)^model_frechet_loglogistic$alpha_b)
      model_frechet_loglogistic$c_weight * exp(-(model_frechet_loglogistic$theta_b / v)^model_frechet_loglogistic$alpha_b) / F_u
    } else {
      u_val = ((v - model_frechet_loglogistic$u) / model_frechet_loglogistic$theta_t)^model_frechet_loglogistic$gamma_t
      model_frechet_loglogistic$c_weight + (1 - model_frechet_loglogistic$c_weight) * (u_val / (1 + u_val))
    }
  })
}
d_spl1 <- function(t) {
  sapply(t, function(v) {
    if(v <= model_frechet_loglogistic$u) {
      F_u = exp(-(model_frechet_loglogistic$theta_b / model_frechet_loglogistic$u)^model_frechet_loglogistic$alpha_b)
      pdf_b = exp(log(model_frechet_loglogistic$alpha_b) + model_frechet_loglogistic$alpha_b * log(model_frechet_loglogistic$theta_b) - (model_frechet_loglogistic$alpha_b + 1) * log(v) - (model_frechet_loglogistic$theta_b / v)^model_frechet_loglogistic$alpha_b)
      model_frechet_loglogistic$c_weight * pdf_b / F_u
    } else {
      v_shift = v - model_frechet_loglogistic$u
      pdf_t = exp(log(model_frechet_loglogistic$gamma_t) + model_frechet_loglogistic$gamma_t * (log(v_shift) - log(model_frechet_loglogistic$theta_t)) - log(v_shift) - 2 * log(1 + (v_shift / model_frechet_loglogistic$theta_t)^model_frechet_loglogistic$gamma_t))
      (1 - model_frechet_loglogistic$c_weight) * pdf_t
    }
  })
}
plot_pit(x, p_spl1, "Spliced: Frechet-Loglogistic")
print(paste("A-D test:", calc_ad(x, p_spl1)$decision))
print(calc_aic(x, d_spl1, k = 4))
plot_exposure(x, p_spl1, lev_func = NULL, model_name = "Spliced: Frechet-Loglogistic")

# 6. Spliced Frechet-Pareto
cat("\n--- Diagnostics: SPLICED (Frechet & Pareto) ---\n")
p_spl2 <- function(t) {
  sapply(t, function(v) {
    if(v <= model_frechet_pareto$u) {
      F_u = exp(-(model_frechet_pareto$theta_b / model_frechet_pareto$u)^model_frechet_pareto$alpha_b)
      model_frechet_pareto$c_weight * exp(-(model_frechet_pareto$theta_b / v)^model_frechet_pareto$alpha_b) / F_u
    } else {
      model_frechet_pareto$c_weight + (1 - model_frechet_pareto$c_weight) * (1 - (model_frechet_pareto$lambda_t / (v - model_frechet_pareto$u + model_frechet_pareto$lambda_t))^model_frechet_pareto$alpha_t)
    }
  })
}
d_spl2 <- function(t) {
  sapply(t, function(v) {
    if(v <= model_frechet_pareto$u) {
      F_u = exp(-(model_frechet_pareto$theta_b / model_frechet_pareto$u)^model_frechet_pareto$alpha_b)
      pdf_b = exp(log(model_frechet_pareto$alpha_b) + model_frechet_pareto$alpha_b * log(model_frechet_pareto$theta_b) - (model_frechet_pareto$alpha_b + 1) * log(v) - (model_frechet_pareto$theta_b / v)^model_frechet_pareto$alpha_b)
      model_frechet_pareto$c_weight * pdf_b / F_u
    } else {
      v_shift = v - model_frechet_pareto$u
      pdf_t = (model_frechet_pareto$alpha_t * model_frechet_pareto$lambda_t^model_frechet_pareto$alpha_t) / ((v_shift + model_frechet_pareto$lambda_t)^(model_frechet_pareto$alpha_t + 1))
      (1 - model_frechet_pareto$c_weight) * pdf_t
    }
  })
}
plot_pit(x, p_spl2, "Spliced: Frechet-Pareto")
print(paste("A-D test:", calc_ad(x, p_spl2)$decision))
print(calc_aic(x, d_spl2, k = 4))
plot_exposure(x, p_spl2, lev_func = NULL, model_name = "Spliced: Frechet-Pareto")

# 7. EM EXPONENTIAL
cat("\n--- Diagnostics: EM EXPONENTIAL ---\n")
p_em <- function(t) pexp(t, rate = 1 / model_em_exp$theta)
d_em <- function(t) dexp(t, rate = 1 / model_em_exp$theta)
plot_pit(x, p_em, "EM Exponential")
print(paste("A-D test:", calc_ad(x, p_em)$decision))
print(calc_aic(x, d_em, k = 1)) 
plot_exposure(x, p_em, lev_func = NULL, model_name = "EM Exponential")

# Zamknięcie zapisu do pliku PDF
dev.off()