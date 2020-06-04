generatePoissonData = function(psi, num_reps, beta, gamma) {
  y1 = rpois(num_reps, lambda=gamma*psi + beta)
  y2 = rpois(num_reps, lambda=beta*t)
  y3 = rpois(num_reps, lambda=gamma*u)
  Y = matrix(data=c(y1, y2, y3), nrow=3, ncol=length(y1))
  return(Y)
}
Y = generatePoissonData(3, 2, 0.321, 0.175)
print(Y)
for (i in 1:ncol(Y)) {
  print(Y[1, i])
}

power_lrt = function(psis, beta, gamma, alpha, num_reps) {
  res = c()
  for (psi in psis) {
    Y = generatePoissonData(psi, num_reps, beta, gamma)
    
    num_rejections = 0
    for (i in 1:ncol(Y)) {
      y = Y[,i]
      beta_p = beta_p(psi, y)
      gamma_p = gamma_p(psi, y)
  
      chi_alpha = qchisq(1-alpha, 1)
      q = (chi_alpha/2 + gamma_p*psi) * (1/log((gamma_p*psi/beta_p)+1))
      
      if (!is.nan(q) && !is.na(q) && y[1] > q) {
        num_rejections = num_rejections + 1
      }
    }
    res = append(res, num_rejections / num_reps)
  }
  return(res)
}
power_lrt(psis=c(0.004), beta=0.291, gamma=0.175, alpha=0.1, num_reps=100)

power_r = function(psis, beta, gamma, alpha, num_reps) {
  res = c()
  for (psi in psis) {
    Y = generatePoissonData(psi, num_reps, beta, gamma)
    
    num_rejections = 0
    for (i in 1:ncol(Y)) {
      y = Y[,i]
      beta_p = beta_p(psi, y)
      gamma_p = gamma_p(psi, y)
      psi_global_mle = getGlobalMLE(psi, y)
      
      rate_mle = gamma_p*psi_global_mle + beta_p
      psi_0 = beta_p
      
      z_alpha = qnorm(1-alpha)
      
      q = ((z_alpha^2)/2 + gamma_p*psi_global_mle) * (1/log(rate_mle/psi_0))
      if (!is.nan(q) && !is.na(q) && y[1] > q) {
        num_rejections = num_rejections + 1
      }
    }
    print(psi_global_mle)
    
    res = append(res, num_rejections / num_reps)
    
  }
  return(res)
}
power_r(psis=c(4), beta=0.296, gamma=0.175, alpha=0.1, num_reps=100)

getCriticalValueViaBoot = function(psi, psi_global_mle, Y,
                                   alpha, num_boot_samples) {
  y1_boot = sample(Y[1,], num_boot_samples, replace=TRUE)
  y2_boot = sample(Y[2,], num_boot_samples, replace=TRUE)
  y3_boot = sample(Y[3,], num_boot_samples, replace=TRUE)
  r_vals=c()
  for (i in 1:num_boot_samples) {
    y_boot = c(y1_boot[i], y2_boot[i], y3_boot[i])
    r_vals = append(r_vals, r_psi(psi, psi_global_mle, y_boot))
  }
  r_vals = sort(r_vals)
  
  r_ecdf = 1:length(r_vals)/length(r_vals)

  crit = r_vals[which(r_ecdf >= 1-alpha)[1]]
}

power_r_boot = function(psis, beta, gamma, alpha, 
                        num_reps, num_boot_samples) {
  res = c()
  for (psi in psis) {
    Y = generatePoissonData(psi, num_reps, beta, gamma)
    
    num_rejections = 0
    crit = getCriticalValueViaBoot(psi, 4.021, Y,
                                   alpha, num_boot_samples)
    num_rejected = 0
    for (i in 1:ncol(Y)) {
      y = Y[,i]
      beta_p = beta_p(psi, y)
      gamma_p = gamma_p(psi, y)
      psi_global_mle = getGlobalMLE(psi, y)
      
      rate_mle = gamma_p*psi_global_mle + beta_p
      psi_0 = beta_p
    
      r_val = r_psi(psi, psi_global_mle, y)
      if (!is.nan(r_val) && !is.na(r_val) && r_val > crit) {
        num_rejected = num_rejected + 1
      }
    }
    res = append(res, num_rejected/num_reps)
    
  }
  return(res)
}
power_r_boot(psi=c(0.5), beta=0.296, gamma=0.175, alpha=0.1, 
             num_reps=200, num_boot_samples=500)


power_r_star = function(psi, alpha, num_boot_samples) {
    
}
power_r_star(0.0005, 0.1, 100)



###########################################
#### Plotting the various power functions
###########################################



psis = seq(0.001, 10, 0.5)
b = 0.296
g = 0.175
a = 0.1
psi_string = expression(psi) 
title = paste("Power fns; beta=", b, ", gamma=", g, ", alph=", a, sep="")
plot(psis, power_lrt(psis, b, g, a, num_reps=500), type="l", col="blue", 
     xlab=psi_string, ylab="Power", main=title,
     ylim=c(0,1.2))
legend(x=6, y=0.3, 
       legend=c("LRT", "r", "r.boot"), 
       col=c("blue", "red", "purple"),
       lty=1:1, cex=0.8
)

lines(psis, power_r(psis, b, g, a, num_reps=300), col="red")

lines(psis, power_r_boot(psis, b, g, a, 
                         num_reps=200, 
                         num_boot_samples = 500), col="purple")



#### Zoom in closer to zero
psis = seq(0.0001, 2.0001, 0.1)
title = paste("Power fns; beta=", b, ", gamma=", g, ", alph=", a, sep="")
plot(psis, power_lrt(psis, b, g, a, num_reps=500), 
     type="l", col="blue", 
     ylim=c(0, 0.9),
     xlab=psi_string, ylab="Power", main=title)
lines(psis, power_r(psis, b, g, a, num_reps=300), col="red")
lines(psis, power_r_boot(psis, b, g, a, 
                         num_reps=200, 
                         num_boot_samples = 500), col="purple")
legend(x=1.0, y=0.2, 
       legend=c("LRT", "r", "r.boot", "Feldman-Cousins"), 
       col=c("blue", "red", "purple", "dark orange"),
       lty=1:1, cex=0.8
)
