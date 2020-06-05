# Likelihood Ratio Test  
lrt = function(psi_0, y) {
  beta_p = beta_p(psi_0+0.000001, y)
  gamma_p = gamma_p(psi_0, y)
  psi_global_mle = getGlobalMLE(psi_0, y)
  
  #print("PSI GLOBAL MLE")
  #print(psi_global_mle)
  
  diff = l(psi_global_mle, beta_p, gamma_p, y) - l(psi_0, beta_p, gamma_p, y)
  return(2 * diff)
}

# Signed root likelihood statistic
r_psi = function(psi_0, y) {
  psi_global_mle = getGlobalMLE(psi_0, y)
  diff = l_p(psi_global_mle, y) - l_p(psi_0, y)
  
  #print("PROFILE LIKELIHOODS")
  #print(l_p(psi_global_mle, y))
  #print(l_p(psi_0, y))
  
  val = sign(psi_global_mle - psi_0) * sqrt(2 * diff)
  return(val)
}

# Barndorff-Nielsen's r* statistic
r_star_psi = function(psi, y) {
  return(r_psi + (log(10) - log(5)) / r_psi)
}

# Helper function to the r_star statistic
q_psi = function(psi) {
  
}
