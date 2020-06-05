generatePoissonData = function(psi, num_reps, beta, gamma) {
  y1 = rpois(num_reps, lambda=gamma*psi + beta)
  y2 = rpois(num_reps, lambda=beta*t)
  y3 = rpois(num_reps, lambda=gamma*u)
  Y = matrix(data=t(c(y1, y2, y3)), nrow=3, ncol=length(y1), byrow=TRUE)
  return(Y)
}

# 
# power() is a helper function - it is passed the relevant parameters,
# the TestStat (such as r(psi)) and the criticalValue, and it 
# computes the power of the test statistic for the given value of psi
power = function(psi, psi_0, beta, gamma, TestStat,
                 num_reps, criticalValue) {
  Y = generatePoissonData(psi, num_reps, beta, gamma)
  num_rejected = 0
  for (i in 1:ncol(Y)) {
    y = Y[,i]
    testStatValue = TestStat(psi_0, y)
    
    if (!is.nan(testStatValue) && 
        !is.na(testStatValue) && 
        testStatValue > criticalValue) {
      num_rejected = num_rejected + 1
    }
  }
  return(num_rejected/num_reps)
}

# calculates power of the LRT for all values in psis
power_lrt = function(psis, psi_0, beta, gamma, alpha, num_reps) {
  res = c()
  chisq_alpha = qchisq(1-alpha, 1)
  
  for (psi in psis) {
    Y = generatePoissonData(psi, num_reps, beta, gamma)
    res = append(res, power(psi, psi_0, beta, gamma, 
                            lrt, num_reps, chisq_alpha))
  }
  return(res)
}

#quick test of power_lrt
power_lrt(psis=c(3), psi_0=0, 
          beta=0.291, gamma=0.175, alpha=0.1, 
          num_reps=1000)

power_r = function(psis, psi_0, beta, gamma, alpha, num_reps) {
  res = c()
  z_alpha = qnorm(1-alpha)
  
  for (psi in psis) {
    Y = generatePoissonData(psi, num_reps, beta, gamma)
    res = append(res, power(psi, psi_0, beta, gamma, 
                            r_psi, num_reps, z_alpha))
  }
  return(res)
}

#quick test of power_r
power_r(psis=c(0.05), psi_0=0, 
        beta=0.296, gamma=0.175, alpha=0.1, num_reps=200)

getCriticalValueViaBoot = function(psi_0, Y, alpha, 
                                   num_boot_samples) {
  y1_boot = sample(Y[1,], num_boot_samples, replace=TRUE)
  y2_boot = sample(Y[2,], num_boot_samples, replace=TRUE)
  y3_boot = sample(Y[3,], num_boot_samples, replace=TRUE)

  r_vals=c()
  for (i in 1:num_boot_samples) {
    y_boot = c(y1_boot[i], y2_boot[i], y3_boot[i])
    r_vals = append(r_vals, r_psi(psi_0, y_boot))
  }
  r_vals = sort(r_vals)
  
  r_ecdf = 1:length(r_vals)/length(r_vals)

  crit = r_vals[which(r_ecdf >= 1-alpha)[1]]
  return(crit)
}

power_r_boot = function(psis, psi_0, beta, gamma, alpha,
                        num_reps, num_boot_samples) {
  res = c()
  for (psi in psis) {
    Y = generatePoissonData(psi_0, num_reps=num_reps, beta, gamma)
    crit = getCriticalValueViaBoot(psi_0, Y,
                                   alpha, num_boot_samples)
    res = append(res, power(psi, psi_0, beta, gamma, 
                            r_psi, num_reps, crit))
  }
  return(res)
}
# quick test of bootstrap r(psi) power function
power_r_boot(psis=c(5), psi_0=0, beta=0.296, gamma=0.175, alpha=0.1, 
               num_reps=100, num_boot_samples=200)

power_r_star = function(psi, alpha, num_boot_samples) {
    ##### TODO
}

###########################################
#### Plotting the various power functions
###########################################

psis = seq(0.001, 10, 0.5)
psi_0 = 0
b = 0.296
g = 0.175
a = 0.1
psi_string = expression(psi) 
title = paste("Power fns; beta=", b, ", gamma=", g, ", alph=", a, sep="")
plot(psis, power_lrt(psis, psi_0, b, g, a, num_reps=200), type="l", col="blue", 
     xlab=psi_string, ylab="Power", main=title,
     ylim=c(0,1.2))
legend(x=6, y=0.3, 
       legend=c("LRT", "r", "r.boot"), 
       col=c("blue", "red", "purple"),
       lty=1:1, cex=0.8
)

lines(psis, power_r(psis, psi_0, b, g, a, num_reps=300), col="red")

lines(psis, power_r_boot(psis, psi_0, b, g, a, 
                         num_reps=100, 
                         num_boot_samples = 200), col="purple")



###########################################
#### Plot that is closer to zero
###########################################

psis = seq(0, 1.0001, 0.05)
title = paste("Power fns; beta=", b, ", gamma=", g, ", alph=", a, sep="")
plot(psis, power_lrt(psis, psi_0, b, g, a, num_reps=200), 
     type="l", col="blue", 
     ylim=c(0, 0.4),
     xlab=psi_string, ylab="Power", main=title)
lines(psis, power_r(psis, psi_0, b, g, a, num_reps=300), col="red")
lines(psis, power_r_boot(psis, psi_0, b, g, a, 
                         num_reps=200, 
                         num_boot_samples = 300), col="purple")
legend(x=0, y=0.4, 
       legend=c("LRT", "r", "r.boot"), 
       col=c("blue", "red", "purple"),
       lty=1:1, cex=0.8
)
