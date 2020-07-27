# Likelihood Ratio Test
lrt = function(psi_0, psi_mle, L, Y, ...) {
  beta_constr = beta_c(Y)
  gamma_constr = gamma_c(Y)

  val = 2 * (l(psi_mle, beta_constr, gamma_constr, Y) - 
             l(psi_0, beta_constr, gamma_constr, Y))
  return(c(val, 1-pchisq(val, df=1)))
}

# Signed root likelihood statistic
r_psi = function(psi_0, Y, ...) {
  psi_global_mle = getGlobalMLE(Y)[1]
  diff = l_p(psi_global_mle, Y) - l_p(psi_0, Y)
  val = sign(psi_global_mle - psi_0) * sqrt(2 * diff)
  return(c(val, 1-pnorm(val)))
}

# Same as r_psi, but we utilize the bootstrap to
# standardize r_psi, in hopes that the standardized
# statistic has a distribution closer to N(0,1)
r_std = function(psi_0, psi_mle, L, Y, ...) {
  sampleSize = nrow(Y)
  bootSize = 25
  
  r_val = r(psi_0, psi_mle, L, Y)[1]
  mu_boot = boot(mean, r, psi_0, L, Y, 
                 sampleSize, bootSize)
  sd_boot = boot(sd, r, psi_0, L, Y, 
                 sampleSize, bootSize)
  
  val = (r_val - mu_boot) / sd_boot
  pval = 1-pnorm(val)
  return(c(val, 1-pnorm(val)))
}

# Signed root likelihood statistic (general function)
r = function(theta_0, theta_mle, L, X, ...) {
  diff = log(L(theta_mle, X) / L(theta_0, X))
  if (is.nan(diff)) {
    print("NaN")
    return(NaN)
  }
  if (is.infinite(diff)) {
    print("Infinite")
  }
  val = sign(theta_mle - theta_0) * sqrt(2 * diff)
  return(c(val, 1-pnorm(val)))
}

# Bootstrapped version of the 
# signed root likelihood statistic
r_psi_boot = function(psi_0, Y, BootSamples) {
  r = r_psi(psi_0, Y)
  return(c(r, getBootstrapPValue(psi_0, BootSamples, r, r_psi)))
}

# Barndorff-Nielsen's r* statistic
r_star = function(psi_0, psi_mle, L, Y, ...) {
  MLEs = getGlobalMLE(Y)
  beta_constr = beta_c(Y)
  gamma_constr = gamma_c(Y)
  
  dl_dtheta_hat_MLE = dl_dtheta_hat(MLEs[1], MLEs[2], MLEs[3], Y)
  dl_dtheta_hat_constrained = 
    dl_dtheta_hat(psi_0, beta_constr, gamma_constr, Y)
  
  first_row = dl_dtheta_hat_MLE - dl_dtheta_hat_constrained
  second_row = dl_dlambda_dtheta_hat(psi_0, beta_constr, gamma_constr, Y)
  m = rbind(dl_dtheta_hat_MLE - dl_dtheta_hat_constrained, 
            dl_dlambda_dtheta_hat(psi_0, beta_constr, gamma_constr, Y))
  det_dl_term = det(m)
  
  j_lambda = j(psi_0, beta_constr, gamma_constr, Y)[2:3, 2:3]
  det_lambda = det(j_lambda)
  
  j = j(MLEs[1], MLEs[2], MLEs[3], Y)
  det_j = det(j)
  
  if (!is.nan(det_j) && det_j <= 0) {
    j_term = NaN
  } else {
    j_term = sqrt(det_lambda) * sqrt(det_j)
  }
  
  v_p = det_dl_term / j_term
  r_psi = r(psi_0, psi_mle, L, Y, ...)
  
  quotient = v_p/r_psi[1]
  
  if (is.nan(quotient)) {
    ### This case occurs when y_2 or y_3 is 0
    ### (and thus, beta's and/or gamma's MLE is 0)
    #quotient = 0.00001
    return(c(NaN, NaN))
  } else if (is.infinite(quotient) && quotient == Inf) {
    ### 
    ### TODO: This appears to the big remaining
    ### problem case (along with the -Inf case below).
    
    ### There appears to be some issues coming from the determinant
    ### of the matrix m above. I would assume it's supposed to be 
    ### positive definite, but I am not always getting that. 
    
    ### It can also occur if r_psi is 0
    
    ### Setting quotient to a small number in this case to 
    ### minimize the contribution of this case to increased rejection 
    ### rate in repeated sampling
    
    #quotient = 0.00001
    return(c(NaN, NaN))
  } else if (is.infinite(quotient) && quotient == -Inf) {
    ### This case occurs when y_1 is 0 (and thus the psi MLE is 0, and 
    ### r(psi) is 0.
    ### For some reason, the determinant of m is always negative in this case
    ### as well... This seems fishy. 
    
    return(c(NaN, NaN))
    #quotient = 0.00001
  }
  
  val = r_psi[1] + ((log(quotient)) / r_psi[1])
  return(c(val, 1-pnorm(val)))
}

r_star_std = function(psi_0, psi_mle, L, Y, ...) {
  sampleSize = nrow(Y)
  bootSize = 20
  
  r_val = r_star(psi_0, psi_mle, L, Y)[1]
  mu_boot = boot(mean, r_star, psi_0, L, Y, 
                 sampleSize, bootSize)
  sd_boot = boot(sd, r_star, psi_0, L, Y, 
                 sampleSize, bootSize)
  
  val = (r_val - mu_boot) / sd_boot
  pval = 1-pnorm(val)
  return(c(val, 1-pnorm(val)))
}

# Barndorff-Nielsen's r* statistic with a 
# bootstrapped p-value
r_star_boot = function(psi_0, Y, BootSamples) {
  r_star = r_star_psi(psi_0, Y)
  
  ### TODO: remove after underlying numeric 
  ### issues with r*(psi) are 100% resolved. 
  ### Mostly unneeded now, but this makes for easier debugging. 
  if (is.nan(r_star[1])) {
    return(c(r_star, NULL))
  }
  
  return(c(r_star, 
           getBootstrapPValue(psi_0, BootSamples, r_star, r_star_psi)))
}
