# Libraries
library(moments)
library(bayesmeta)
library(matlib)
library(DescTools)
library(evd)
library(CASdatasets)

# Helper functions

my_kurtosis <- function(u){
  z = (u[4]-4*u[3]*u[1]+6*u[2]*u[1]^2-3*u[1]^4)/(u[2]-u[1]^2)^2
}  

my_skewness <- function(u){
  z = (u[3]-3*u[2]*u[1]+2*u[1]^3)/(u[2]-u[1]^2)^(3/2)
} 

# Log-likelihood functions

gamma_loglik_function <- function(p,x){
  z = sum(log(dgamma(x, shape = p[1], rate = p[2])))
  return(z)
}

weibull_loglik_function <- function(p,x){
  z = sum(log(dweibull(x, shape = p[1], scale = p[2])))
  return(z)
}

pareto_loglik_function <- function(p,x){
  z = sum(log(dlomax(x, shape = p[1], scale = p[2])))
  return(z)
}

gpareto_loglik_function <- function(p,x){
  z = sum(log(dgpd(x, loc = 0, shape = p[1], scale = p[2])))
  return(z)
}

# Estimation methods

fit_gamma <- function(x, m, t, u){
  if(t == "Lower"){x = x + u}
  y = x/mean(x)
  
  loglik_function <- function(p){
    if(t == "Complete"){q = 0}
    if(t == "Upper"){q = length(x)*log(pgamma(u/mean(x), shape = p[1], rate = p[2]))}
    if(t == "Lower"){
      q = length(x)*log(1-pgamma(u/mean(x), shape = p[1], rate = p[2]))}
    z = -gamma_loglik_function(p,y)+q
    return(z)
  }
  
  # Estimation
  alpha_mm = mean(y)^2/var(y)
  beta_mm = mean(y)/var(y)
  
  if(t == "Lower"){
    for (it in c(1:100)){
      dist = pgamma(u/mean(x), shape = alpha_mm, rate = beta_mm)
      z = sum(rgeom(length(x), 1-dist))
      uu = runif(z, 0, 1)
      yy = qgamma(uu*dist, shape = alpha_mm, rate = beta_mm)
      y_new = c(yy,y)
      model_gamma_0 = fit_gamma(y_new, "ML", "Complete",0)
      alpha_mm = model_gamma_0$alpha
      beta_mm = model_gamma_0$beta
    }
  }
  
  model_ml = optim(p = c(alpha_mm,beta_mm), loglik_function)
  alpha_ml = model_ml$par[1]
  beta_ml = model_ml$par[2]
  
  if(m == "MM"){
    alpha = alpha_mm
    beta = beta_mm/mean(x)
  }
  
  if(m == "ML"){
    alpha = alpha_ml
    beta = beta_ml/mean(x)
  }
  
  # Standard errors for "Complete"
  beta_s = beta*mean(x)
  a_alpha = trigamma(alpha)
  a_beta = alpha/beta_s^2
  a_alpha_beta = -1/beta_s
  
  i_fisher = rbind(c(a_alpha,a_alpha_beta), c(a_alpha_beta,a_beta))
  i_fisher = inv(i_fisher)/length(x)
  
  i_fisher[1,2] = i_fisher[1,2]/mean(x)
  i_fisher[2,1] = i_fisher[2,1]/mean(x)
  i_fisher[2,2] = i_fisher[2,2]/(mean(x))^2
  
  # Results
  moments = c(alpha/beta,
              alpha*(alpha+1)/beta^2,
              alpha*(alpha+1)*(alpha+2)/beta^3,
              alpha*(alpha+1)*(alpha+2)*(alpha+3)/beta^4)
  
  moments_truncate <-function(a,k){
    z = pgamma(a*beta, shape = alpha+k, rate = 1)
    return(z)
  }
  
  if(t == "Complete"){
    moments = c(moments[1],
                sqrt(moments[2]-moments[1]^2),
                my_skewness(moments), my_kurtosis(moments))
    
    quantiles = qgamma(c(0.01,0.05,c(1:9)/10,0.95,0.99), shape = alpha, rate = beta)
    distribution_values = pgamma(sort(x), shape = alpha, rate = beta)
    truncated_values = moments_truncate(sort(x),1)
    excess_values = ((moments[1]*(1-moments_truncate(sort(x),1))+10^(-50))/
                       (1-pgamma(sort(x), shape = alpha, rate = beta)+10^(-50))-sort(x))
  }
  
  if(t == "Upper"){
    exp_value = moments[1]
    moments = c(moments[1]*moments_truncate(u,1),
                moments[2]*moments_truncate(u,2),
                moments[3]*moments_truncate(u,3),
                moments[4]*moments_truncate(u,4))
    moments = moments/pgamma(u,shape = alpha, rate = beta)
    moments = c(moments[1],
                sqrt(moments[2]-moments[1]^2),
                my_skewness(moments), my_kurtosis(moments))
    
    quantiles = qgamma(c(0.01,0.05,c(1:9)/10,0.95,0.99)*pgamma(u,shape = alpha, rate = beta),
                       shape = alpha, rate = beta)
    distribution_values = pgamma(sort(x), shape = alpha, rate = beta)/pgamma(u, shape = alpha, rate = beta)
    truncated_values = moments_truncate(sort(x),1)/moments_truncate(u,1)
    excess_values = ((exp_value*(moments_truncate(u,1)-moments_truncate(sort(x),1))+10^(-50))/
                       (pgamma(u, shape = alpha, rate = beta)-pgamma(sort(x), shape = alpha, rate = beta)+10^(-50))-sort(x))
  }
  
  if(t == "Lower"){
    x = x - u
    exp_value = moments[1]  
    moments = c(moments[1]*(1-moments_truncate(u,1)),
                moments[2]*(1-moments_truncate(u,2)),
                moments[3]*(1-moments_truncate(u,3)),
                moments[4]*(1-moments_truncate(u,4)))
    moments = moments/(1-pgamma(u,shape = alpha, rate = beta))
    moments = c(moments[1]-u,
                sqrt(moments[2]-moments[1]^2),
                my_skewness(moments), my_kurtosis(moments))
    
    quantiles = (qgamma(c(0.01,0.05,c(1:9)/10,0.95,0.99)*(1-pgamma(u,shape = alpha, rate = beta))+pgamma(u,shape = alpha, rate = beta),
                        shape = alpha, rate = beta) - u)
    distribution_values = (pgamma(sort(x)+u, shape = alpha, rate = beta)-pgamma(u,shape = alpha, rate = beta))/(1-pgamma(u, shape = alpha, rate = beta))
    truncated_values = ((exp_value*(moments_truncate(u+sort(x),1)-moments_truncate(u,1))
                         -u*(pgamma(u+sort(x),shape = alpha, rate = beta)-pgamma(u,shape = alpha, rate = beta)))/
                          (exp_value*(1-moments_truncate(u,1))
                           -u*(1-pgamma(u,shape = alpha, rate = beta))))
    excess_values = ((exp_value*(1-moments_truncate(sort(x)+u,1))+10^(-50))/
                       (1-pgamma(sort(x)+u, shape = alpha, rate = beta)+10^(-50))-sort(x)-u)
  }
  
  z = list("alpha" = alpha,
           "beta" = beta,
           "moments" = moments,
           "quantiles" = quantiles,
           "distribution_values" = distribution_values,
           "excess_values" = excess_values,
           "truncated_values" = truncated_values,
           "cov" = i_fisher)
  return(z)
}

fit_weibull <- function(x, m, t, u){
  if(t == "Lower"){x = x + u}
  y = x/mean(x)
  
  loglik_function <- function(p){
    if(t == "Complete"){q = 0}
    if(t == "Upper"){q = length(x)*log(pweibull(u/mean(x), shape = p[1], scale = p[2]))}
    if(t == "Lower"){
      q = length(x)*log(1-pweibull(u/mean(x), shape = p[1], scale = p[2]))}
    z = -weibull_loglik_function(p,y)+q
    return(z)
  }
  
  # Estimation
  tau_function <- function(t){
    z = 2*lgamma(1+1/t)-lgamma(1+2/t)+log(1+var(y)/mean(y)^2)
    return(z)
  }
  
  tau_mm = uniroot(tau_function, c(0.001,1000))$root
  c_mm = gamma(1+1/tau_mm)/mean(y)
  
  if(t == "Lower"){
    for (it in c(1:100)){
      dist = pweibull(u/mean(x), shape = tau_mm, scale = 1/c_mm)
      z = sum(rgeom(length(x), 1-dist))
      uu = runif(z, 0, 1)
      yy = qweibull(uu*dist, shape = tau_mm, scale = 1/c_mm)
      y_new = c(yy,y)
      model_weibull_0 = fit_weibull(y_new, "ML", "Complete",0)
      tau_mm = model_weibull_0$tau
      c_mm = model_weibull_0$c
    }
  }
  
  model_ml = optim(p = c(tau_mm,1/c_mm), loglik_function)
  tau_ml = model_ml$par[1]
  c_ml = 1/model_ml$par[2]
  
  if(m == "MM"){
    tau = tau_mm
    c = c_mm/mean(x)
  }
  
  if(m == "ML"){
    tau = tau_ml
    c = c_ml/mean(x)
  }
  
  # Standard errors for "Complete"
  c_s = c*mean(x)
  a_tau = 1/tau^2+mean((c_s*x/mean(x))^tau*(log(c_s*x/mean(x)))^2)
  a_c = tau/c_s^2+tau*(tau-1)*c_s^(tau-2)*mean((x/mean(x))^tau)
  a_tau_c =-1/c_s+c_s^(tau-1)*mean(tau*log(c_s*x/mean(x))*(x/mean(x))^tau+(x/mean(x))^tau)
  
  i_fisher = rbind(c(a_tau,a_tau_c), c(a_tau_c,a_c))
  i_fisher = inv(i_fisher)/length(x)
  i_fisher[1,2] = i_fisher[1,2]/mean(x)
  i_fisher[2,1] = i_fisher[2,1]/mean(x)
  i_fisher[2,2] = i_fisher[2,2]/(mean(x))^2
  
  # Results
  moments = c(gamma(1+1/tau)/c,
              gamma(1+2/tau)/c^2,
              gamma(1+3/tau)/c^3,
              gamma(1+4/tau)/c^4)
  
  moments_truncate <-function(a,k){
    z = pgamma((a*c)^tau, shape = 1+k/tau, rate = 1)
    return(z)
  }
  
  if(t == "Complete"){
    moments = c(moments[1],
                sqrt(moments[2]-moments[1]^2),
                my_skewness(moments), my_kurtosis(moments))
    quantiles = qweibull(c(0.01,0.05,c(1:9)/10,0.95,0.99), shape = tau, scale = 1/c)
    distribution_values = pweibull(sort(x), shape = tau, scale = 1/c)
    truncated_values = moments_truncate(sort(x),1)
    excess_values = ((moments[1]*(1-moments_truncate(sort(x),1))+10^(-50))/
                       (1-pweibull(sort(x), shape = tau, scale = 1/c)+10^(-50))-sort(x))
  }
  
  if(t == "Upper"){
    exp_value = moments[1]
    moments = c(moments[1]*moments_truncate(u,1),
                moments[2]*moments_truncate(u,2),
                moments[3]*moments_truncate(u,3),
                moments[4]*moments_truncate(u,4))
    moments = moments/pweibull(u,shape = tau, scale = 1/c)
    moments = c(moments[1],
                sqrt(moments[2]-moments[1]^2),
                my_skewness(moments), my_kurtosis(moments))
    quantiles = qweibull(c(0.01,0.05,c(1:9)/10,0.95,0.99)*pweibull(u, shape = tau, scale = 1/c),
                         shape = tau, scale = 1/c)
    distribution_values = pweibull(sort(x), shape = tau, scale = 1/c)/pweibull(u, shape = tau, scale = 1/c)
    truncated_values = moments_truncate(sort(x),1)/moments_truncate(u,1)
    excess_values = ((exp_value*(moments_truncate(u,1)-moments_truncate(sort(x),1))+10^(-50))/
                       (pweibull(u, shape = tau, scale = 1/c)-pweibull(sort(x), shape = tau, scale = 1/c)+10^(-50))-sort(x))
  }
  
  if(t == "Lower"){
    x = x - u
    exp_value = moments[1]  
    moments = c(moments[1]*(1-moments_truncate(u,1)),
                moments[2]*(1-moments_truncate(u,2)),
                moments[3]*(1-moments_truncate(u,3)),
                moments[4]*(1-moments_truncate(u,4)))
    moments = moments/(1-pweibull(u,shape = tau, scale = 1/c))
    moments = c(moments[1]-u,
                sqrt(moments[2]-moments[1]^2),
                my_skewness(moments), my_kurtosis(moments))
    quantiles = (qweibull(c(0.01,0.05,c(1:9)/10,0.95,0.99)*(1-pweibull(u, shape = tau, scale = 1/c))+pweibull(u, shape = tau, scale = 1/c),
                          shape = tau, scale = 1/c) - u)
    distribution_values = (pweibull(sort(x)+u, shape = tau, scale = 1/c)-pweibull(u, shape = tau, scale = 1/c))/(1-pweibull(u, shape = tau, scale = 1/c))
    truncated_values = ((exp_value*(moments_truncate(u+sort(x),1)-moments_truncate(u,1))
                         -u*(pweibull(u+sort(x), shape = tau, scale = 1/c)-pweibull(u, shape = tau, scale = 1/c)))/
                          (exp_value*(1-moments_truncate(u,1))
                           -u*(1-pweibull(u,shape = tau, scale = 1/c))))
    excess_values = ((exp_value*(1-moments_truncate(sort(x)+u,1))+10^(-50))/
                       (1-pweibull(sort(x)+u, shape = tau, scale = 1/c)+10^(-50))-sort(x)-u)
  }
  
  z = list("tau" = tau,
           "c" = c,
           "moments" = moments,
           "quantiles" = quantiles,
           "distribution_values" = distribution_values,
           "excess_values" = excess_values,
           "truncated_values" = truncated_values,
           "cov" = i_fisher)
  return(z)
}

fit_pareto <- function(x, m, t, u){
  if(t == "Lower"){x = x + u}
  y = x/mean(x)
  
  loglik_function <- function(p){
    if(t == "Complete"){q = 0}
    if(t == "Upper"){q = length(x)*log(plomax(u/mean(x), shape = p[1], scale = p[2]))}
    if(t == "Lower"){
      q = length(x)*log(1-plomax(u/mean(x), shape = p[1], scale = p[2]))}
    z = -pareto_loglik_function(p,y)+q
    return(z)
  }
  
  # Estimation
  alpha_mm = (2*mean(y)^2-2*(var(y)+mean(y)^2))/(2*mean(y)^2-(var(y)+mean(y)^2))
  lambda_mm = mean(y)*(alpha_mm-1)
  
  if(t == "Lower"){
    mp = fit_pareto(y - u/mean(x), "ML", "Complete",0)
    alpha_mm = mp$alpha
    lambda_mm = mp$lambda}
  
  model_ml = optim(p = c(max(alpha_mm,1),max(lambda_mm,10^(-10))), loglik_function)
  alpha_ml = model_ml$par[1]
  lambda_ml = model_ml$par[2]
  
  if(m == "MM"){
    alpha = ifelse(alpha_mm>2,alpha_mm,"NA")
    lambda = ifelse(alpha_mm>2,lambda_mm*mean(x),"NA")
  }
  
  if(m == "ML"){
    alpha = alpha_ml
    lambda = lambda_ml*mean(x)
  }
  
  # Standard errors for "Complete"
  lambda_s = lambda/mean(x)
  a_alpha = 1/alpha^2
  a_lambda = alpha/lambda_s^2-(alpha+1)*mean(1/(x/mean(x)+lambda_s)^2)
  a_alpha_lambda = -1/lambda_s+mean(1/(x/mean(x)+lambda_s))
  
  i_fisher = rbind(c(a_alpha,a_alpha_lambda), c(a_alpha_lambda,a_lambda))
  i_fisher = inv(i_fisher)/length(x)
  i_fisher[1,2] = i_fisher[1,2]*mean(x)
  i_fisher[2,1] = i_fisher[2,1]*mean(x)
  i_fisher[2,2] = i_fisher[2,2]*(mean(x))^2
  
  # Results
  moments = c(lambda^1*factorial(1)/(alpha-1),
              lambda^2*factorial(2)/(alpha-1)/(alpha-2),
              lambda^3*factorial(3)/(alpha-1)/(alpha-2)/(alpha-3),
              lambda^4*factorial(4)/(alpha-1)/(alpha-2)/(alpha-3)/(alpha-4))
  
  moments_truncate <-function(a,k){
    z = pbeta(a/(a+lambda), shape1 = k+1, shape2 = alpha-k)
    return(z)
  }
  
  if(t == "Complete"){
    moments = c(ifelse(alpha>1,moments[1],0),
                ifelse(alpha>2,sqrt(moments[2]-moments[1]^2),0),
                ifelse(alpha>3,my_skewness(moments),0), 
                ifelse(alpha>4,my_kurtosis(moments),0))
    
    quantiles = qlomax(c(0.01,0.05,c(1:9)/10,0.95,0.99), shape = alpha, scale = lambda)
    distribution_values = plomax(sort(x), shape = alpha, scale = lambda)
    truncated_values = moments_truncate(sort(x),1)
    excess_values = ((moments[1]*(1-moments_truncate(sort(x),1))+10^(-50))/
                       (1-plomax(sort(x), shape = alpha, scale = lambda)+10^(-50))-sort(x))
  }
  
  if(t == "Upper"){
    exp_value = moments[1]
    moments = c(moments[1]*ifelse(alpha>1,moments_truncate(u,1),0),
                moments[2]*ifelse(alpha>2,moments_truncate(u,2),0),
                moments[3]*ifelse(alpha>3,moments_truncate(u,3),0),
                moments[4]*ifelse(alpha>4,moments_truncate(u,4),0))
    moments = moments/plomax(u, shape = alpha, scale = lambda)
    moments = c(ifelse(alpha>1,moments[1],0),
                ifelse(alpha>2,sqrt(moments[2]-moments[1]^2),0),
                ifelse(alpha>3,my_skewness(moments),0), 
                ifelse(alpha>4,my_kurtosis(moments),0))
    
    quantiles = qlomax(c(0.01,0.05,c(1:9)/10,0.95,0.99)*plomax(u, shape = alpha, scale = lambda),
                       shape = alpha, scale = lambda)
    distribution_values = plomax(sort(x), shape = alpha, scale = lambda)/plomax(u, shape = alpha, scale = lambda)
    truncated_values = moments_truncate(sort(x),1)/moments_truncate(u,1)
    excess_values = ((exp_value*(moments_truncate(u,1)-moments_truncate(sort(x),1))+10^(-50))/
                       (plomax(u, shape = alpha, scale = lambda)-plomax(sort(x), shape = alpha, scale = lambda)+10^(-50))-sort(x))
  }
  
  if(t == "Lower"){
    x = x - u
    exp_value = moments[1]  
    moments = c(moments[1]*ifelse(alpha>1,1-moments_truncate(u,1),0),
                moments[2]*ifelse(alpha>2,1-moments_truncate(u,2),0),
                moments[3]*ifelse(alpha>3,1-moments_truncate(u,3),0),
                moments[4]*ifelse(alpha>4,1-moments_truncate(u,4),0))
    moments = moments/(1-plomax(u, shape = alpha, scale = lambda))
    moments = c(ifelse(alpha>1,moments[1]-u,0),
                ifelse(alpha>2,sqrt(moments[2]-moments[1]^2),0),
                ifelse(alpha>3,my_skewness(moments),0), 
                ifelse(alpha>4,my_kurtosis(moments),0))
    
    quantiles = (qlomax(c(0.01,0.05,c(1:9)/10,0.95,0.99)*(1-plomax(u,shape = alpha, scale = lambda))+plomax(u,shape = alpha, scale = lambda),
                        shape = alpha, scale = lambda) - u)
    distribution_values = (plomax(sort(x)+u, shape = alpha, scale = lambda)-plomax(u, shape = alpha, scale = lambda))/(1-plomax(u, shape = alpha, scale = lambda))
    truncated_values = ((exp_value*(moments_truncate(u+sort(x),1)-moments_truncate(u,1))
                         -u*(plomax(u+sort(x), shape = alpha, scale = lambda)-plomax(u, shape = alpha, scale = lambda)))/
                          (exp_value*(1-moments_truncate(u,1))
                           -u*(1-plomax(u, shape = alpha, scale = lambda))))
    excess_values = ((exp_value*(1-moments_truncate(sort(x)+u,1))+10^(-50))/
                       (1-plomax(sort(x)+u, shape = alpha, scale = lambda)+10^(-50))-sort(x)-u)
  }
  
  z = list("alpha" = alpha,
           "lambda" = lambda,
           "moments" = moments,
           "quantiles" = quantiles,
           "distribution_values" = distribution_values,
           "excess_values" = excess_values,
           "truncated_values" = truncated_values,
           "cov" = i_fisher)
  return(z)
}

fit_negbin <- function(x, v, m){
  # Estimation
  freq_mm = sum(x)/sum(v)
  sigma_sq = sum((x/v-freq_mm)^2*v)/(length(x)-1)
  gamma_mm = freq_mm^2/(sigma_sq-freq_mm)
  
  gamma_function <- function(gamma){
    z = sum(digamma(gamma*v+x)*v-digamma(gamma*v)*v+v*log(1/(1+sum(x)/sum(v)/gamma)))
    return(z)
  }
  
  gamma_ml = uniroot(gamma_function, c(gamma_mm/100,gamma_mm*100))$root
  p_ml = 1/(1+sum(x)/sum(v)/gamma_ml)
  freq_ml = (1/p_ml-1)*gamma_ml
  
  if(m == "MM"){
    freq = freq_mm
    gamma = gamma_mm
  }
  
  if(m == "ML"){
    freq = freq_ml
    gamma = gamma_ml
  }
  
  if(m == "Poiss"){
    freq = freq_ml
    gamma = 10^10
  }
  
  # Goodness-of-fit test
  p_residuals = (x/v-freq)/sqrt(freq*(1+freq/gamma)/v)
  plot(p_residuals, main = "Pearson residuals", ylab = "", xlab = "Observation", pch = 19, cex = 1, col = "blue")
  abline(h = 0, col = "red", lwd = 2)
  t = 1-pchisq(sum(p_residuals^2), df = length(x)-1)
  
  z = list("freq" = freq,
           "gamma" = gamma,
           "chi_sq_test_pvalue" = t)
  return(z)
}

# Validation methods

summary_fit <- function(model,x){
  summary_estimation = rbind(
    c(mean(x),model$moments[1]),
    c(sd(x),model$moments[2]),
    c(skewness(x), model$moments[3]),
    c(kurtosis(x),model$moments[4]),
    cbind(quantile(x,c(0.01,0.05,c(1:9)/10,0.95,0.99), type = 6), model$quantiles))
  
  summary_estimation = data.frame(summary_estimation)
  colnames(summary_estimation) = c("Data","Model")
  rownames(summary_estimation) = c("Mean", "Std dev", "Skewness", "Kurtosis", 
                                   "Q1%", "Q5%", "Q10%", "Q20%", "Q30%", "Q40%", "Q50%", 
                                   "Q60%", "Q70%", "Q80%", "Q90%",
                                   "Q95%", "Q99%")
  z = summary_estimation
  return(z)
}

distribution_fit <- function(model,x){
  emp_distr = ecdf(x)
  sort_x = sort(x)
  
  par(mfrow = c(2,2))
  plot(x = sort_x, y = emp_distr(sort_x), col = "blue", type = "s", lwd = 2,
       xlab = "Response", ylab = "", main = "Distribution function")
  lines(x = sort_x, y = model$distribution_values, col = "red", lwd = 2)
  
  quantile_residuals = qnorm(model$distribution_values,0,1)
  normal_quantiles = qnorm(c(1:length(x))/(length(x)+1),0,1)
  
  plot(quantile_residuals ~ normal_quantiles, col = "blue", pch = 19, cex = 0.5,
       xlim = c(-5,5), ylim =c(-5,5),
       ylab = "Theoretical normal quantile", xlab = "Sample normal quantile", 
       main = "QQ plot")
  abline(0,1,col = "red", lwd = 2)
  
  plot(x = log(sort_x), y = log(emp_distr(sort_x)), col = "blue", cex = 0.5, pch = 19,
       xlab = "Logged response", ylab = "", main = "Logged distribution function")
  lines(x = log(sort_x), y = log(model$distribution_values), col = "red", lwd = 2)
  
  plot(x = log(sort_x), y = log(1-emp_distr(sort_x)), col = "blue", cex = 0.5, pch = 19,
       xlab = "Logged response", ylab = "", main = "Logged survival function")
  lines(x = log(sort_x), y = log(1-model$distribution_values), col = "red", lwd = 2)
  
  par(mfrow = c(1,1))
}

actuarial_fit <- function(model,x){
  sort_x = sort(x)
  excess_claims = sapply(2:length(x), function(i){ mean(sort_x[i:length(x)])-sort_x[i-1] })
  truncated_claims = sapply(1:length(x), function(i){ sum(sort_x[1:i]) })/sum(x)
  
  par(mfrow = c(1,2))
  plot(excess_claims ~ sort_x[-length(x)], col = "blue", cex = 0.5, pch = 19,
       xlab = "Threshold", ylab = "", main = "Mean excess loss")
  lines(x = sort_x[-length(x)], y = model$excess_values[-length(x)], col = "red", lwd = 2)
  abline(v = quantile(x,0.95), col = "green", lwd = 2)
  abline(v = quantile(x,0.99), col = "green", lwd = 2)
  
  u = c(1:length(x))/length(x)
  plot(truncated_claims ~ u, col = "blue", cex = 0.5, pch = 19,
       xlab = "Threshold", ylab = "", main = "Loss size index")
  lines(x = u, y = model$truncated_values, col = "red", lwd = 2)
  abline(v = 0.95, col = "green", lwd = 2)
  abline(v = 0.99, col = "green", lwd = 2)
  abline(h = sum(sort_x[1:round(0.95*length(x))])/sum(x), col = "black", lwd = 2)
  abline(h = sum(sort_x[1:round(0.99*length(x))])/sum(x), col = "black", lwd = 2)
  par(mfrow = c(1,1))
}

statistical_fit <-function(model,x){
  emp_distr = ecdf(x)
  sort_x = sort(x)
  
  ks_values = pmax(abs(emp_distr(sort_x)-model$distribution_values), 
                   abs(emp_distr(sort_x)-model$distribution_values-1/length(x)))
  ad_values = ks_values/sqrt(model$distribution_values*(1-model$distribution_values))
  
  par(mfrow = c(1,2))
  plot(ks_values ~ log(sort_x), col = "blue", cex = 0.5, pch = 19,
       xlab = "Logged response", ylab = "", main = "Kolmogorov-Smirnoff test")
  
  plot(ad_values ~ log(sort_x), col = "blue", cex = 0.5, pch = 19,
       xlab = "Logged response", ylab = "", main = "Weighted Kolmogorov-Smirnoff test")
  par(mfrow = c(1,1))
  
  z = ks.test(sort_x,model$distribution_values)$p.value
  return(z)
}

empirical_tail <- function(x){
  emp_distr = ecdf(x)
  sort_x = sort(x)
  
  par(mfrow = c(2,3))
  plot(x = -log(1-emp_distr(sort_x)), y = sort_x, col = "blue", cex = 0.5, pch = 19,
       ylab = "Response", xlab = "Log tail", main = "Exponential QQ-plot")
  plot(x = log(-log(1-emp_distr(sort_x))), y = log(sort_x), col = "blue", cex = 0.5, pch = 19,
       ylab = "Logged response", xlab = "Log-log tail", main = "Weibull QQ-plot")
  plot(x = -log(1-emp_distr(sort_x)), y = log(sort_x), col = "blue", cex = 0.5, pch = 19,
       ylab = "Logged response", xlab = "Log tail", main = "Pareto QQ-plot")
  
  derivative = sapply(c(1:length(x)), function(i){
    -sum(sort_x[i:length(x)]-sort_x[i])/
      sum(log(1-c(i:length(x))/(length(x)+1))-log(1-i/(length(x)+1))) })
  
  plot(x = sort_x, y = derivative, col = "blue", cex = 0.5, pch = 19,
       ylab = "", xlab = "Response", main = "Exponential derivative plot")
  
  derivative = sapply(c(1:length(x)), function(i){
    sum(log(sort_x[i:length(x)])-log(sort_x[i]))/
      sum(log(-log(1-c(i:length(x))/(length(x)+1)))-log(-log(1-i/(length(x)+1)))) })
  
  plot(x = sort_x, y = derivative, col = "blue", cex = 0.5, pch = 19,
       ylab = "", xlab = "Logged response", main = "Weibull derivative plot")
  
  derivative = sapply(c(1:length(x)), function(i){
    -sum(log(sort_x[i:length(x)])-log(sort_x[i]))/
      sum(log(1-c(i:length(x))/(length(x)+1))-log(1-i/(length(x)+1))) })
  
  plot(x = sort_x, y = derivative, col = "blue", cex = 0.5, pch = 19,
       ylab = "", xlab = "Logged response", main = "Pareto derivative plot")
  par(mfrow = c(1,1))
}

pareto_plots <- function(y_var, u){
  par(mfrow= c(1,3))
  pareto_claims = sort(y_var)
  q = c(1:length(pareto_claims))/(length(pareto_claims)+1)
  
  plot(log(1-q[which(q>u)]) ~ log(pareto_claims[which(q>u)]), pch = 20,
       main = "Log tail of distribution", xlab = "Log severity", ylab = "", cex = 1.25, col = "blue")
  
  theta_index = log(length(pareto_claims)+1-c(2:length(pareto_claims)))/log(length(pareto_claims)) 
  hill_claims = sapply(round(u*length(pareto_claims)):length(pareto_claims), function(i){
    length(pareto_claims[i:length(pareto_claims)])/sum(log(pareto_claims[i:length(pareto_claims)]/pareto_claims[i-1]))})
  
  plot(hill_claims ~ theta_index[round(u*length(pareto_claims)):length(pareto_claims)], pch = 20, 
       main = "Hill estimator of tail index alpha", xlab = "Theta", ylab = "", cex = 1.25, col = "blue")
  
  excess_claims = sapply(round(u*length(pareto_claims)):length(pareto_claims), function(i){
    mean(pareto_claims[i:length(pareto_claims)])-pareto_claims[i-1]})
  
  plot(excess_claims ~ pareto_claims[(round(u*length(pareto_claims))-1):(length(pareto_claims)-1)], pch = 20,
       main = "Mean excess plot", xlab = "Threshold", ylab = "", cex = 1.25, col = "blue")
  par(mfrow= c(1,1))
}

fit_gpd <- function(x) {
  alpha = c()
  theta = c()
  sort_x = sort(x)
  model_ml = fit_pareto(x, "ML", "Complete", 0)
  
  init_1 = 1 / model_ml$alpha
  init_2 = model_ml$lambda * init_1 / mean(x)
  
  for (i in 1:(length(x) - 10)) {
    y = sort_x[which(sort_x > sort_x[i])] - sort_x[i]
    loglik_function <- function(p) {
      z = -gpareto_loglik_function(p, y / mean(y))
      return(z)
    }
    
    opt_result <- try(optim(p = c(init_1, init_2), loglik_function), silent = TRUE)
    if (inherits(opt_result, "try-error")) {
      break
    } else {
      model_ml <- opt_result
    }
    
    init_1 = model_ml$par[1]
    init_2 = model_ml$par[2]
    alpha = c(alpha, 1 / init_1)
    theta = c(theta, log(length(x) - i) / log(length(x)))
  }
  
  plot(alpha ~ theta, pch = 20, main = "GPD estimator", xlab = "Theta", ylab = "Tail index alpha", cex = 1.25, col = "blue")
}

fit_frechet <- function(x) {
  p_val = quantile(x, 0.25, names = FALSE)
  q_val = quantile(x, 0.75, names = FALSE)
  g_val = log(log(4)) / log(log(4/3))
  theta_init = exp((g_val * log(q_val) - log(p_val)) / (g_val - 1))
  alpha_init = log(log(4)) / (log(theta_init) - log(p_val))
  
  loglik_function <- function(param) {
    alpha <- param[1]
    theta <- param[2]
    if(alpha <= 0 || theta <= 0) return(1e10)
    log_f = log(alpha) + alpha * log(theta) - (alpha + 1) * log(x) - (theta / x)^alpha
    return(-sum(log_f))
  }
  
  model_ml = optim(p = c(alpha_init, theta_init), fn = loglik_function, hessian = TRUE)
  alpha = model_ml$par[1]
  theta = model_ml$par[2]
  i_fisher = tryCatch(solve(model_ml$hessian), error = function(e) matrix(NA, 2, 2))
  
  get_moment <- function(k) {
    if (alpha > k) return(theta^k * base::gamma(1 - k/alpha))
    else return(Inf)
  }
  
  raw_moments = c(get_moment(1), get_moment(2), get_moment(3), get_moment(4))
  variance = if(alpha > 2) (raw_moments[2] - raw_moments[1]^2) else Inf
  std_dev = if(variance != Inf && variance > 0) sqrt(variance) else Inf
  
  moments = c(raw_moments[1], std_dev, if(alpha > 3) my_skewness(raw_moments) else Inf, if(alpha > 4) my_kurtosis(raw_moments) else Inf)
  
  probs = c(0.01, 0.05, c(1:9)/10, 0.95, 0.99)
  quantiles = theta * (-log(probs))^(-1/alpha)
  distribution_values = exp(-(theta / sort(x))^alpha)
  
  z = list("alpha" = alpha, "theta" = theta, "moments" = moments, "quantiles" = quantiles, "distribution_values" = distribution_values, "cov" = i_fisher)
  return(z)
}

fit_loglogistic <- function(x) {
  p = quantile(x, 0.25, names = FALSE)
  q = quantile(x, 0.75, names = FALSE)
  gamma_init = 2 * log(3) / (log(q) - log(p))
  theta_init = exp((log(q) + log(p)) / 2)
  
  loglik_function <- function(param) {
    gamma <- param[1]
    theta <- param[2]
    if(gamma <= 0 || theta <= 0) return(1e10)
    log_f = log(gamma) + gamma * (log(x) - log(theta)) - log(x) - 2 * log(1 + (x/theta)^gamma)
    return(-sum(log_f))
  }
  
  model_ml = optim(p = c(gamma_init, theta_init), fn = loglik_function, hessian = TRUE)
  gamma = model_ml$par[1]
  theta = model_ml$par[2]
  i_fisher = tryCatch(solve(model_ml$hessian), error = function(e) matrix(NA, 2, 2))
  
  get_moment <- function(k) {
    if (gamma > k) return(theta^k * base::gamma(1 + k/gamma) * base::gamma(1 - k/gamma))
    else return(Inf)
  }
  raw_moments = c(get_moment(1), get_moment(2), get_moment(3), get_moment(4))
  variance = if(gamma > 2) (raw_moments[2] - raw_moments[1]^2) else Inf
  std_dev = if(variance != Inf && variance > 0) sqrt(variance) else Inf
  
  moments = c(raw_moments[1], std_dev, if(gamma > 3) my_skewness(raw_moments) else Inf, if(gamma > 4) my_kurtosis(raw_moments) else Inf)
  
  probs = c(0.01, 0.05, c(1:9)/10, 0.95, 0.99)
  quantiles = theta * (probs^(-1) - 1)^(-1/gamma)
  u_val = (sort(x)/theta)^gamma
  distribution_values = u_val / (1 + u_val)
  
  z = list("gamma" = gamma, "theta" = theta, "moments" = moments, "quantiles" = quantiles, "distribution_values" = distribution_values, "cov" = i_fisher)
  return(z)
}

fit_spliced_frechet_loglogis <- function(x, u = 40000) {
  x_small = x[x <= u]
  x_large = x[x > u] - u
  c_weight = length(x_small) / length(x)
  
  # Frechet Estimation
  ll_body <- function(p) {
    alpha = p[1]; theta = p[2]
    if(alpha <= 0 || theta <= 0) return(1e10)
    log_f = log(alpha) + alpha * log(theta) - (alpha + 1) * log(x_small) - (theta / x_small)^alpha
    log_F_u = -(theta / u)^alpha
    return(-sum(log_f - log_F_u))
  }
  opt_b = optim(c(1, mean(x_small)), ll_body)
  alpha_b = opt_b$par[1]
  theta_b = opt_b$par[2]
  
  # Loglogistic Estimation
  ll_tail <- function(p) {
    gamma = p[1]; theta_t = p[2]
    if(gamma <= 0 || theta_t <= 0) return(1e10)
    log_f = log(gamma) + gamma * (log(x_large) - log(theta_t)) - log(x_large) - 2 * log(1 + (x_large/theta_t)^gamma)
    return(-sum(log_f))
  }
  opt_t = optim(c(1, mean(x_large)), ll_tail)
  gamma_t = opt_t$par[1]
  theta_t = opt_t$par[2]
  
  # CDF
  F_u_val = exp(-(theta_b / u)^alpha_b)
  dist_vals = numeric(length(x))
  sorted_x = sort(x)
  
  for(i in 1:length(sorted_x)) {
    v = sorted_x[i]
    if(v <= u) {
      dist_vals[i] = c_weight * exp(-(theta_b / v)^alpha_b) / F_u_val
    } else {
      u_val = ((v - u) / theta_t)^gamma_t
      dist_vals[i] = c_weight + (1 - c_weight) * (u_val / (1 + u_val))
    }
  }
  
  # Quantiles
  probs = c(0.01, 0.05, c(1:9)/10, 0.95, 0.99)
  quants = numeric(length(probs))
  for(i in 1:length(probs)) {
    p = probs[i]
    if(p <= c_weight) {
      p_adj = p * F_u_val / c_weight
      quants[i] = theta_b * (-log(p_adj))^(-1/alpha_b)
    } else {
      p_tail = (p - c_weight) / (1 - c_weight)
      quants[i] = u + theta_t * (p_tail^(-1) - 1)^(-1/gamma_t)
    }
  }
  
  # Moments
  f_spliced <- function(v) {
    ifelse(v <= u, 
           c_weight * (alpha_b * theta_b^alpha_b * v^(-alpha_b - 1) * exp(-(theta_b/v)^alpha_b)) / F_u_val,
           (1 - c_weight) * (gamma_t * ((v-u)/theta_t)^gamma_t) / ((v-u) * (1 + ((v-u)/theta_t)^gamma_t)^2))
  }
  
  get_moment <- function(k) {
    res = tryCatch(integrate(function(v) (v^k) * f_spliced(v), 0, Inf, subdivisions=2000)$value, error=function(e) Inf)
    return(res)
  }
  
  raw_moments = c(get_moment(1), get_moment(2), get_moment(3), get_moment(4))
  var_val = raw_moments[2] - raw_moments[1]^2
  sd_val = if(!is.na(var_val) && var_val > 0 && var_val != Inf) sqrt(var_val) else Inf
  
  skew_val = if(raw_moments[3] != Inf && sd_val != Inf) my_skewness(raw_moments) else Inf
  kurt_val = if(raw_moments[4] != Inf && sd_val != Inf) my_kurtosis(raw_moments) else Inf
  
  moments = c(raw_moments[1], sd_val, skew_val, kurt_val)
  
  z = list("alpha_b" = alpha_b, "theta_b" = theta_b, "gamma_t" = gamma_t, "theta_t" = theta_t, "c_weight" = c_weight, "u" = u, "moments" = moments, "quantiles" = quants, "distribution_values" = dist_vals)
  return(z)
}

fit_spliced_frechet_pareto <- function(x, u = 40000) {
  x_small = x[x <= u]
  x_large = x[x > u] - u
  c_weight = length(x_small) / length(x)
  
  # Frechet Estimation
  ll_body <- function(p) {
    alpha = p[1]; theta = p[2]
    if(alpha <= 0 || theta <= 0) return(1e10)
    log_f = log(alpha) + alpha * log(theta) - (alpha + 1) * log(x_small) - (theta / x_small)^alpha
    log_F_u = -(theta / u)^alpha
    return(-sum(log_f - log_F_u))
  }
  opt_b = optim(c(1, mean(x_small)), ll_body)
  alpha_b = opt_b$par[1]
  theta_b = opt_b$par[2]
  
  # Pareto Estimation
  ll_tail <- function(p) {
    alpha_t = p[1]; lambda_t = p[2]
    if(alpha_t <= 0 || lambda_t <= 0) return(1e10)
    log_f = log(alpha_t) + alpha_t * log(lambda_t) - (alpha_t + 1) * log(x_large + lambda_t)
    return(-sum(log_f))
  }
  opt_t = optim(c(2, mean(x_large)), ll_tail)
  alpha_t = opt_t$par[1]
  lambda_t = opt_t$par[2]
  
  # CDF
  F_u_val = exp(-(theta_b / u)^alpha_b)
  dist_vals = numeric(length(x))
  sorted_x = sort(x)
  
  for(i in 1:length(sorted_x)) {
    v = sorted_x[i]
    if(v <= u) {
      dist_vals[i] = c_weight * exp(-(theta_b / v)^alpha_b) / F_u_val
    } else {
      dist_vals[i] = c_weight + (1 - c_weight) * (1 - (lambda_t / (v - u + lambda_t))^alpha_t)
    }
  }
  
  # Quantiles
  probs = c(0.01, 0.05, c(1:9)/10, 0.95, 0.99)
  quants = numeric(length(probs))
  for(i in 1:length(probs)) {
    p = probs[i]
    if(p <= c_weight) {
      p_adj = p * F_u_val / c_weight
      quants[i] = theta_b * (-log(p_adj))^(-1/alpha_b)
    } else {
      p_tail = (p - c_weight) / (1 - c_weight)
      quants[i] = u + lambda_t * ((1 - p_tail)^(-1/alpha_t) - 1)
    }
  }
  
  # Moments
  f_spliced <- function(v) {
    ifelse(v <= u, 
           c_weight * (alpha_b * theta_b^alpha_b * v^(-alpha_b - 1) * exp(-(theta_b/v)^alpha_b)) / F_u_val,
           (1 - c_weight) * (alpha_t * lambda_t^alpha_t) / ((v - u + lambda_t)^(alpha_t + 1)))
  }
  
  get_moment <- function(k) {
    res = tryCatch(integrate(function(v) (v^k) * f_spliced(v), 0, Inf, subdivisions=2000)$value, error=function(e) Inf)
    return(res)
  }
  
  raw_moments = c(get_moment(1), get_moment(2), get_moment(3), get_moment(4))
  var_val = raw_moments[2] - raw_moments[1]^2
  sd_val = if(!is.na(var_val) && var_val > 0 && var_val != Inf) sqrt(var_val) else Inf
  
  skew_val = if(raw_moments[3] != Inf && sd_val != Inf) my_skewness(raw_moments) else Inf
  kurt_val = if(raw_moments[4] != Inf && sd_val != Inf) my_kurtosis(raw_moments) else Inf
  
  moments = c(raw_moments[1], sd_val, skew_val, kurt_val)
  
  z = list("alpha_b" = alpha_b, "theta_b" = theta_b, "alpha_t" = alpha_t, "lambda_t" = lambda_t, "c_weight" = c_weight, "u" = u, "moments" = moments, "quantiles" = quants, "distribution_values" = dist_vals)
  return(z)
}

fit_em_exp <- function(x, u = 40000) {
  # Data preparation
  x_obs = x[x <= u]
  n = length(x_obs)
  sum_x = sum(x_obs)
  
  # Initial value for EM
  theta_current = mean(x_obs)
  
  # Expectation-Maximization loop
  for (iter in 1:100) {
    # EXPECTATION STEP
    F_u = pexp(u, rate = 1 / theta_current)
    
    # number of missing observations (M)
    M = n * ((1 - F_u) / F_u)
    
    # expected value from memoryless property
    E_X_tail = u + theta_current
    
    # MAXIMIZATION STEP
    theta_next = (sum_x + M * E_X_tail) / (n + M)
    
    # Convergence check
    if (abs(theta_next - theta_current) < 1e-6) {
      theta_current = theta_next
      cat("EM convergence reached in", iter, "iterations.\n")
      break
    }
    theta_current = theta_next
  }
  
  theta = theta_current
  
  # Estimator variance
  var_theta = (theta^2) / (n + M)
  cov_mat = matrix(var_theta, 1, 1)
  
  # Theoretical moments of exponential distribution
  raw_moments = c(theta, 2*theta^2, 6*theta^3, 24*theta^4)
  std_dev = theta
  
  moments = c(raw_moments[1], std_dev, my_skewness(raw_moments), my_kurtosis(raw_moments))
  
  # Quantile vector (matched to summary_fit table)
  probs = c(0.01, 0.05, c(1:9)/10, 0.95, 0.99)
  quantiles = qexp(probs, rate = 1/theta)
  
  # Theoretical CDF for the entire input vector x
  distribution_values = pexp(sort(x), rate = 1/theta)
  
  z = list("theta" = theta, "moments" = moments, "quantiles" = quantiles, "distribution_values" = distribution_values, "cov" = cov_mat)
  return(z)
}

plot_pit <- function(x, p_func, model_name = "Model") {
  u <- p_func(x)
  hist(u, breaks = 20, freq = FALSE, 
       main = paste("PIT -", model_name), 
       xlab = "F(x)", ylab = "Density", 
       col = "lightblue", border = "white")
  abline(h = 1, col = "red", lwd = 2, lty = 2)
}

calc_ad <- function(x, p_func) {
  n <- length(x)
  x_sort <- sort(x)
  
  # CDF function application
  u <- p_func(x_sort)
  
  # Numerical safeguard
  u <- pmax(pmin(u, 1 - 1e-10), 1e-10)
  
  # A^2 statistic calculation
  i <- 1:n
  S <- sum((2 * i - 1) * (log(u) + log(1 - u[n - i + 1])))
  A2 <- -n - S/n
  
  # Classic critical value for 0.05 significance level
  critical_value <- 2.492
  
  # Statistical decision
  reject_H0 <- A2 > critical_value
  decision_msg <- ifelse(reject_H0, 
                         "Reject H0 at 0.05 level", 
                         "Do not reject H0 at 0.05 level")
  
  return(list(
    statistic = A2, 
    critical_value = critical_value,
    reject_H0 = reject_H0, 
    decision = decision_msg
  ))
}

plot_exposure <- function(x, p_func, lev_func = NULL, d_max = quantile(x, 0.95), model_name = "Model") {
  d_seq <- seq(0, d_max, length.out = 100)
  
  # 1. Empirical Curve
  emp_EX <- mean(x)
  emp_G <- sapply(d_seq, function(d) mean(pmin(x, d)) / emp_EX)
  
  # 2. Theoretical Curve
  if (!is.null(lev_func)) {
    theo_EX <- lev_func(Inf)
    theo_G <- sapply(d_seq, function(d) lev_func(d) / theo_EX)
  } else {
    theo_G <- sapply(d_seq, function(d) {
      num <- integrate(function(t) 1 - p_func(t), lower = 0, upper = d, subdivisions=2000, stop.on.error = FALSE)$value
      den <- integrate(function(t) 1 - p_func(t), lower = 0, upper = Inf, subdivisions=2000, stop.on.error = FALSE)$value
      return(num / den)
    })
  }
  
  plot(d_seq, emp_G, type = "l", col = "black", lwd = 2, 
       xlab = "Contract limit (d)", ylab = "G(d)", main = paste("Exposure Curve -", model_name))
  lines(d_seq, theo_G, col = "red", lwd = 2, lty = 2)
  legend("bottomright", legend = c("Empirical", "Theoretical"), 
         col = c("black", "red"), lty = c(1, 2), lwd = 2)
}

calc_aic <- function(x, d_func, k) {
  log_pdf <- log(d_func(x))
  log_pdf[!is.finite(log_pdf)] <- -1e10 
  
  logL <- sum(log_pdf)
  aic <- 2 * k - 2 * logL
  return(c(LogLikelihood = logL, AIC = aic))
}