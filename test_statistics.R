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
  if (is.nan(diff) || is.na(diff)) {
    browser()
  }
  if (diff < 0) {
    print("r_psi : diff is negative")
  }
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
  m = matrix(data=c(
    dl_dtheta_hat_MLE - dl_dtheta_hat_constrained,
    dl_dlambda_dtheta_hat(psi_0, beta_p, gamma_p, y)
  ), nrow=3, ncol=3)
  det_dl_term = det(m)
  
  j_lambda = j(psi_0, beta_p, gamma_p, y)[2:3, 2:3]
  det_lambda = det(j_lambda)
  
  det_j = det(j(MLEs[1], MLEs[2], MLEs[3], y))
  j_term = sqrt(det_lambda) * sqrt(det_j)
  
  v_p = det_dl_term / j_term
  r_psi = r_psi(psi_0, y)
  
  quotient = v_p/r_psi[1]
  
  ### TODO: These checks make me nervous. 
  ### And I don't understand why det_dl_term is so small
  ### so often. It's always right around 0, with 
  ### positive or negative sign with about equal frequency.
  ### I ought to triple-check these calculations. 
  #if (is.nan(quotient) || 
  #    is.infinite(quotient) ||
  #    quotient <= 0) {
  #  quotient = 0.00001
  #}
  ### TODO: If we set the global MLE for psi to 0
  ### when we would otherwise have a negative
  ### global MLE for psi, r(psi) will have a value
  ### of 0 (and then a p-value of 0.5). This causes
  ### quotient to be Inf, and the pValue to be 0, 
  ### which should be fine I think? Maybe I should remove
  ### this check.
  #if (quotient == Inf) {
  #  quotient = 1000000000
  #}
  val = r_psi[1] + ((log(quotient)) / r_psi[1])
  return(c(val, 1-pnorm(val)))
}

# Barndorff-Nielsen's r* statistic
r_star_psi_2 = function(psi_0, y, ...) {
  MLEs = getGlobalMLE(y)
  beta_p = beta_p_0(y)
  gamma_p = gamma_p_0(y)
  
  dl_dtheta_hat_MLE = dl_dtheta_hat(MLEs[1], MLEs[2], MLEs[3], y)
  dl_dtheta_hat_constrained = dl_dtheta_hat(psi_0, beta_p, gamma_p, y)
  m = matrix(data=c(
    dl_dtheta_hat_MLE - dl_dtheta_hat_constrained,
    dl_dlambda_dtheta_hat(psi_0, beta_p, gamma_p, y)
  ), nrow=3, ncol=3)
  det_dl_term = det(m)
  
  j_lambda = j(psi_0, beta_p, gamma_p, y)[2:3, 2:3]
  det_lambda = det(j_lambda)
  
  det_j = det(j(MLEs[1], MLEs[2], MLEs[3], y))
  j_term = sqrt(det_lambda) * sqrt(det_j)
  
  v_p = det_dl_term / j_term
  r_psi = r_psi(psi_0, y)
  
  quotient = v_p/r_psi[1]
  
  ### TODO: These checks make me nervous. 
  ### And I don't understand why det_dl_term is so small
  ### so often. It's always right around 0, with 
  ### positive or negative sign with about equal frequency.
  ### I ought to triple-check these calculations. 
  if (is.nan(quotient) || 
    is.infinite(quotient) ||
      quotient <= 0) {
    quotient = 0.00001
  }
  ### TODO: If we set the global MLE for psi to 0
  ### when we would otherwise have a negative
  ### global MLE for psi, r(psi) will have a value
  ### of 0 (and then a p-value of 0.5). This causes
  ### quotient to be Inf, and the pValue to be 0, 
  ### which should be fine I think? Maybe I should remove
  ### this check.
  if (quotient == Inf) {
    quotient = 1000000000
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

# Barndorff-Nielsen's r* statistic with a 
# bootstrapped p-value
r_star_boot_2 = function(psi_0, y, BootSamples) {
  r_star = r_star_psi_2(psi_0, y)
  
  ### TODO: remove after underlying numeric 
  ### issues with r*(psi) are 100% resolved. 
  ### Mostly unneeded now, but this makes for easier debugging. 
  if (is.nan(r_star[1])) {
    return(c(r_star, NULL))
  }
  
  return(c(r_star, getBootstrapPValue(psi_0, BootSamples, r_star, r_star_psi_2)))
}