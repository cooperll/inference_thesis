generatePoissonData = function(psi, num_reps, beta, gamma) {
  y1 = rpois(num_reps, lambda=gamma*psi + beta)
  y2 = rpois(num_reps, lambda=beta*t)
  y3 = rpois(num_reps, lambda=gamma*u)
  Y = matrix(data=t(c(y1, y2, y3)), nrow=3, ncol=length(y1), byrow=TRUE)
  return(Y)
}
getBootstrapPValue = function(psi_0, B, r_obs) {
  count = 0
  for (i in 1:ncol(B)) {
    boot = B[,i]
    r_boot = r_psi(psi_0, boot)
    print("BOOT VALS")
    print(r_boot)
    print(r_obs)
    
    if (!is.nan(r_boot) &&
        !is.na(r_boot) &&
        r_boot >= r_obs) {
      count = count + 1
    }
  }
  #r_vals = sort(r_vals)
  #r_ecdf = 1:length(r_vals)/length(r_vals)
  #crit = r_vals[which(r_ecdf >= 1-alpha)[1]]
  return(count/length(B))
}

#
# power() is a helper function - it is passed the relevant parameters,
# the TestStat (such as r(psi)) and the criticalValue, and it
# computes the power of the test statistic for the given value of psi
power = function(psi, psi_0, beta, gamma, TestStat,
                 num_reps, criticalValue, ...) {
  Y = generatePoissonData(psi, num_reps, beta, gamma)
  num_rejected = 0
  for (i in 1:ncol(Y)) {
    y = Y[,i]
    testStatValue = TestStat(psi_0, y, ...)

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
    res = append(res, power(psi, psi_0, beta, gamma,
                            lrt, num_reps, chisq_alpha))
  }
  return(res)
}

power_r = function(psis, psi_0, beta, gamma, alpha, num_reps) {
  res = c()
  z_alpha = qnorm(1-alpha)

  for (psi in psis) {
    res = append(res, power(psi, psi_0, beta, gamma,
                            r_psi, num_reps, z_alpha))
  }
  return(res)
}

power_c = function(psis, psi_0, beta, gamma, alpha, 
                   num_reps, num_boot_samples, nth) {
  res = c()
  for (psi in psis) {
    Y = generatePoissonData(psi_0, num_reps=num_reps, beta, gamma)
    crit = getCriticalValueViaBoot(psi_0, Y,
                                   alpha, num_boot_samples)
    print(crit)
    res = append(res, power(psi, psi_0, beta, gamma,
                            c_psi, num_reps, crit, nth))
  }
  return(res)
}

power_r_boot = function(psis, psi_0, beta, gamma, alpha,
                        num_reps, num_boot_samples) {
  beta_c = beta_p_0(y)
  gamma_c = gamma_p_0(y)
  
  res = c()
  for (psi in psis) {
    Y = generatePoissonData(psi, num_reps, beta, gamma)
    
    num_rejected = 0
    for (i in 1:ncol(Y)) {
      B = generatePoissonData(psi_0, num_boot_samples, beta_c, gamma_c)
      
      y = Y[,i]
      r_obs = r_psi(psi_0, y)
      print("R_OBS")
      print(r_obs)
      pValue = getBootstrapPValue(psi_0, B, r_obs)

      if (!is.nan(pValue) &&
          !is.na(pValue) &&
          pValue < alpha) {
        num_rejected = num_rejected + 1
      }
    }
    res = append(res, num_rejected/num_reps)  
  }
  return(res)
}
testY = generatePoissonData(1, 10, 0.321, 0.175)
testY[,4]
# quick test of bootstrap r(psi) power function
power_r_boot(psis=c(2), psi_0=0, beta=0.296, gamma=0.175, alpha=0.1,
             num_reps=5, num_boot_samples=5)

power_r_star = function(psis, psi_0, beta, gamma, alpha, 
                        num_reps) {
  res = c()
  z_alpha = qnorm(1-alpha)
  
  for (psi in psis) {
    Y = generatePoissonData(psi, num_reps, beta, gamma)
    res = append(res, power(psi, psi_0, beta, gamma,
                            r_star_psi, num_reps, z_alpha))
  }
  return(res)
}

power_r_star_boot = function(psis, psi_0, beta, gamma, alpha,
                        num_reps, num_boot_samples) {
  res = c()
  for (psi in psis) {
    Y = generatePoissonData(psi_0, num_reps=num_reps, beta, gamma)
    crit = getCriticalValueViaBoot(psi_0, Y, alpha, num_boot_samples)
    res = append(res, power(psi, psi_0, beta, gamma,
                            r_star_psi, num_reps, crit))
  }
  return(res)
}

power_fc = function(psis, psi_0, beta, gamma, alpha, 
                        num_reps) {
  res = c()
  z_alpha = qnorm(1-alpha)
  
  for (psi in psis) {
    Y = generatePoissonData(psi, num_reps, beta, gamma)
    res = append(res, power(psi, psi_0, beta, gamma,
                            fc_psi, num_reps, z_alpha))
  }
  return(res)
}


#####################################
##### Tests of power functions
#####################################

#quick test of power_lrt
power_lrt(psis=c(3), psi_0=0, beta=0.291, gamma=0.175, alpha=0.1,
          num_reps=1000)

#quick test of power_r
power_r(psis=c(0.05), psi_0=0, beta=0.296, gamma=0.175, alpha=0.1, 
        num_reps=200)

#quick test of power_r
power_c(psis=c(0, 0.05, 1, 5), psi_0=0,
        beta=0.296, gamma=0.175, alpha=0.1, 
        num_reps=10000, num_boot_samples=10000, nth=15)

# quick test of bootstrap r(psi) power function
power_r_boot(psis=c(2), psi_0=0, beta=0.296, gamma=0.175, alpha=0.1,
             num_reps=100, num_boot_samples=200)

# quick test of r*(psi) power function
power_r_star(psis=c(5), psi_0=0, beta=0.296, gamma=0.175, alpha=0.1,
             num_reps=100)

# quick test of bootstrap r*(psi) power function
power_r_star_boot(psis=c(5), psi_0=0, beta=0.296, gamma=0.175, alpha=0.1,
             num_reps=100, num_boot_samples=200)

# quick test of Feldman-Cousins power function
power_fc(psis=c(5), psi_0=0, beta=0.296, gamma=0.175, alpha=0.1,
             num_reps=100)

###########################################
#### Plotting the various power functions
###########################################

psis = seq(0, 10, 0.5)
psi_0 = 0
b = 0.296
g = 0.175
a = 0.1
psi_string = expression(psi)
title = paste("Power fns; beta=", b, ", gamma=", g, ", alph=", a, sep="")
plot(psis, power_lrt(psis, psi_0, b, g, a, num_reps=10000), type="l", col="blue",
     xlab=psi_string, ylab="Power", main=title,
     ylim=c(0,1.2))
legend(x=0, y=1.2,
       legend=c("LRT", "r", "r.boot", "r*"),
       col=c("blue", "red", "purple", "dark green"),
       lty=1:1, cex=0.8
)

lines(psis, power_r(psis, psi_0, b, g, a, num_reps=5000), col="red")

lines(psis, power_r_boot(psis, psi_0, b, g, a,
                         num_reps=5000,
                         num_boot_samples = 5000), col="purple")

lines(psis, power_r_star(psis, psi_0, b, g, a,
                         num_reps=5000), col="dark green")

lines(psis, power_fc(psis, psi_0, b, g, a,
                         num_reps=5000), col="dark orange")


###########################################
#### Plot that is closer to zero
###########################################

psis = seq(0, 1.0001, 0.05)
title = paste("Power fns; beta=", b, ", gamma=", g, ", alph=", a, sep="")
plot(psis, power_lrt(psis, psi_0, b, g, a, num_reps=10000),
     type="l", col="blue",
     ylim=c(0, 0.4),
     xlab=psi_string, ylab="Power", main=title)
lines(psis, power_r(psis, psi_0, b, g, a, num_reps=5000), col="red")
lines(psis, power_r_boot(psis, psi_0, b, g, a,
                         num_reps=2500,
                         num_boot_samples = 5000), col="purple")
legend(x=0, y=0.4,
       legend=c("LRT", "r", "r.boot", "r*", "Feldman-Cousins"),
       col=c("blue", "red", "purple", "dark green", "dark orange"),
       lty=1:1, cex=0.8
      )
