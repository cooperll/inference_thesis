### getBootstrapPValue 
###
### psi_0 is the parameter of interest assuming the null hypothesis
### bootData is the data set produced via parametric bootstrap
### obsValue is the observed value of the test statistic 
###    produced via the original data
### TestStat is a *function*. It is the test statistic (i.e. one
###    of the other functions listed in this file)
getBootstrapPValue = function(psi_0, bootData, obsValue, TestStat) {
  count = 0
  total = 0
  for (i in 1:ncol(bootData)) {
    boot = bootData[,i]
    bootValue = TestStat(psi_0, boot)
    
    ### TODO: remove all these obscene checks 
    ### for bad values coming from the test statistic.
    ### At the moment, r*(psi) is to blame. 
    if (is.infinite(bootValue[1]) || 
        is.nan(bootValue[1]) || 
        is.infinite(obsValue[1]) ||
        is.nan(obsValue[1])) {
      next
    }
    total = total + 1
    
    if (bootValue[1] >= obsValue[1]) {
      count = count + 1
    }
  }
  return(count/total)
}

# Likelihood Ratio Test
lrt = function(psi_0, y, ...) {
  beta_p = beta_p_0(y)
  gamma_p = gamma_p_0(y)
  psi_global_mle = getGlobalMLE(y)[1]

  val = 2 * (l(psi_global_mle, beta_p, gamma_p, y) - l(psi_0, beta_p, gamma_p, y))
  return(c(val, 1-pchisq(val, df=1)))
}

# Signed root likelihood statistic
r_psi = function(psi_0, y, ...) {
  psi_global_mle = getGlobalMLE(y)[1]
  diff = l_p(psi_global_mle, y) - l_p(psi_0, y)
  val = sign(psi_global_mle - psi_0) * sqrt(2 * diff)
  return(c(val, 1-pnorm(val)))
}

# Bootstrapped version of the 
# signed root likelihood statistic
r_psi_boot = function(psi_0, y, BootSamples) {
  r = r_psi(psi_0, y)
  return(c(r, getBootstrapPValue(psi_0, BootSamples, r, r_psi)))
}

# Barndorff-Nielsen's r* statistic
r_star_psi = function(psi_0, y, ...) {
  MLEs = getGlobalMLE(y)
  beta_p = beta_p_0(y)
  gamma_p = gamma_p_0(y)
  
  dl_dtheta_hat_MLE = dl_dtheta_hat(MLEs[1], MLEs[2], MLEs[3], y)
  dl_dtheta_hat_constrained = dl_dtheta_hat(psi_0, beta_p, gamma_p, y)
  
  first_row = dl_dtheta_hat_MLE - dl_dtheta_hat_constrained
  second_row = dl_dlambda_dtheta_hat(psi_0, beta_p, gamma_p, y)
  m = rbind(dl_dtheta_hat_MLE - dl_dtheta_hat_constrained, 
            dl_dlambda_dtheta_hat(psi_0, beta_p, gamma_p, y))
  det_dl_term = det(m)
  
  j_lambda = j(psi_0, beta_p, gamma_p, y)[2:3, 2:3]
  det_lambda = det(j_lambda)
  
  j = j(MLEs[1], MLEs[2], MLEs[3], y)
  det_j = det(j)
  
  if (!is.nan(det_j) && det_j <= 0) {
    j_term = NaN
  } else {
    j_term = sqrt(det_lambda) * sqrt(det_j)
  }
  
  v_p = det_dl_term / j_term
  r_psi = r_psi(psi_0, y)
  
  quotient = v_p/r_psi[1]
  
  if (is.nan(quotient)) {
    ### This case occurs when y_2 or y_3 is 0
    ### (and thus, beta's and/or gamma's MLE is 0)
    quotient = 0.00001
  } else if (is.infinite(quotient) && quotient == Inf) {
    ### 
    ### TODO: This appears to the big remaining
    ### problem case (along with the -Inf case below).
    
    ### There appears to be some issues coming from the determinant
    ### of the matrix m above. I would assume it's supposed to be 
    ### positive definite, but I am not always getting that. 
    
    ### Setting quotient to a small number in this case to 
    ### minimize the contribution of this case to increased rejection 
    ### rate in repeated sampling
    quotient = 0.00001
  } else if (is.infinite(quotient) && quotient == -Inf) {
    ### This case occurs when y_1 is 0 (and thus the psi MLE is 0, and 
    ### r(psi) is 0.
    ### For some reason, the determinant of m is always negative in this case
    ### as well... This seems fishy. 
    quotient = 0.00001
  }
  
  val = r_psi[1] + ((log(quotient)) / r_psi[1])
  return(c(val, 1-pnorm(val)))
}

# Barndorff-Nielsen's r* statistic with a 
# bootstrapped p-value
r_star_boot = function(psi_0, y, BootSamples) {
  r_star = r_star_psi(psi_0, y)
  
  ### TODO: remove after underlying numeric 
  ### issues with r*(psi) are 100% resolved. 
  ### Mostly unneeded now, but this makes for easier debugging. 
  if (is.nan(r_star[1])) {
    return(c(r_star, NULL))
  }
  
  return(c(r_star, getBootstrapPValue(psi_0, BootSamples, r_star, r_star_psi)))
}
