options(scipen=99)
library(CASdatasets)

# Load frecomfire dataset
data(frecomfire)

x <- frecomfire$ClaimCost2007
log_x <- log(x)

if(!dir.exists("plots")) {
  dir.create("plots")
}

# Otwarcie zapisu do zbiorczego pliku PDF
pdf("plots/01_eda_and_tail_analysis.pdf", width = 12, height = 8)

### DATA ANALYSIS - empirical histogram and CDF ###
par(mfrow = c(2, 2))

# Histogram and empirical CDF for original observations
hist(x, breaks = 50, probability = TRUE, col = "lightblue",
     main = "Histogram (Original)", xlab = "Claim severity", ylab = "Density")
plot(ecdf(x), main = "Empirical CDF (Original)", 
     xlab = "Claim severity", ylab = "P(X <= x)")

# Histogram and empirical CDF for log observations
hist(log_x, breaks = 50, probability = TRUE, col = "lightgreen",
     main = "Histogram (Logarithm)", xlab = "Log(Claim severity)", ylab = "Density")
plot(ecdf(log_x), main = "Empirical CDF (Logarithm)", 
     xlab = "Log(Claim severity)", ylab = "P(X <= x)")

## Kernel density estimation
# Original observations, density estimation and CDF
emp_density = density(x, kernel = "gaussian", bw = 100, from = min(x))
plot(emp_density, type = "l", col = "blue", lwd = 2, 
     xlab = "Claim severity", ylab = "Density", main = "Kernel density estimator")

dx <- diff(emp_density$x[1:2])
cdf_x <- cumsum(emp_density$y * dx)
plot(emp_density$x, cdf_x, type = "l", col = "red", lwd = 2, 
     xlab = "Claim severity", ylab = "P(X <= x)", main = "Kernel CDF")

# Log observations, density estimation and CDF
emp_density_log = density(log(x), kernel = "gaussian", bw = 0.2, from = min(log(x)))
plot(emp_density_log, type = "l", col = "blue", lwd = 2, 
     xlab = "Log(Claim severity)", ylab = "Density", main = "Kernel density estimator (Logarithm)")

dx_log <- diff(emp_density_log$x[1:2])
cdf_log <- cumsum(emp_density_log$y * dx_log)
plot(emp_density_log$x, cdf_log, type = "l", col = "red", lwd = 2, 
     xlab = "Log(Claim severity)", ylab = "P(X <= x)", main = "Kernel CDF (Logarithm)")

par(mfrow = c(1, 1))

## Descriptive statistics

# PARAMETERS AND FUNCTIONS
bw_orig = 100
bw_log = 0.2
log_severity = log(x)

# CALCULATIONS FOR ORIGINAL VALUES
emp_mean_orig <- mean(x)
emp_sd_orig <- sd(x)
emp_q_orig <- quantile(x, probs = c(0.75, 0.95, 0.99))

kernel_mean_orig <- emp_mean_orig
kernel_sd_orig <- sqrt(var(x) + bw_orig^2)
kernel_distr_orig <- function(q, p) { sum(pnorm(q, mean = x, sd = bw_orig)) / length(x) - p }
lower_orig <- min(x) - 10 * bw_orig
upper_orig <- max(x) + 10 * bw_orig

qq_k_orig_075 <- uniroot(kernel_distr_orig, c(lower_orig, upper_orig), p = 0.75)$root
qq_k_orig_095 <- uniroot(kernel_distr_orig, c(lower_orig, upper_orig), p = 0.95)$root
qq_k_orig_099 <- uniroot(kernel_distr_orig, c(lower_orig, upper_orig), p = 0.99)$root
kernel_q_orig <- c(qq_k_orig_075, qq_k_orig_095, qq_k_orig_099)
names(kernel_q_orig) <- c("75%", "95%", "99%")

# CALCULATIONS FOR LOGARITHMS
emp_mean_log <- mean(log_x)
emp_sd_log <- sd(log_x)
emp_q_log <- quantile(log_x, probs = c(0.75, 0.95, 0.99))

kernel_mean_log <- emp_mean_log
kernel_sd_log <- sqrt(var(log_x) + bw_log^2)
kernel_distr_log <- function(q, p) { sum(pnorm(q, mean = log_x, sd = bw_log)) / length(log_x) - p }
lower_log <- min(log_x) - 10 * bw_log
upper_log <- max(log_x) + 10 * bw_log

qq_k_log_075 <- uniroot(kernel_distr_log, c(lower_log, upper_log), p = 0.75)$root
qq_k_log_095 <- uniroot(kernel_distr_log, c(lower_log, upper_log), p = 0.95)$root
qq_k_log_099 <- uniroot(kernel_distr_log, c(lower_log, upper_log), p = 0.99)$root
kernel_q_log <- c(qq_k_log_075, qq_k_log_095, qq_k_log_099)
names(kernel_q_log) <- c("75%", "95%", "99%")

# RESULTS SUMMARY
summary_table <- data.frame(
  Statistic = c("Expected value", "Standard deviation", "Quantile 0.75", "Quantile 0.95", "Quantile 0.99"),
  Empirical_Original = c(emp_mean_orig, emp_sd_orig, emp_q_orig[1], emp_q_orig[2], emp_q_orig[3]),
  Kernel_Original = c(kernel_mean_orig, kernel_sd_orig, kernel_q_orig[1], kernel_q_orig[2], kernel_q_orig[3]),
  Empirical_Log = c(emp_mean_log, emp_sd_log, emp_q_log[1], emp_q_log[2], emp_q_log[3]),
  Kernel_Log = c(kernel_mean_log, kernel_sd_log, kernel_q_log[1], kernel_q_log[2], kernel_q_log[3])
)

colnames(summary_table) <- c("Statistic", "Empirical (Original)", "Kernel (Original)", "Empirical (Log)", "Kernel (Log)")
print(summary_table, row.names = FALSE)

quantile_40k <- ecdf(x)(40000)
print(quantile_40k)

##### Tail Assessment and Hill Estimator, GPD ####

# 1. Tail assessment (quantile and derivative plots)
empirical_tail(x)

# 2. Assessment of Pareto-type tail presence
pareto_plots(y_var = x, u = 0.9)

# 3. Generalized Pareto Distribution (GPD) estimator
fit_gpd(x)

# Zamknięcie pliku PDF
dev.off()