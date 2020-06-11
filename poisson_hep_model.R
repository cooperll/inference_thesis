# initial values as seen in Davison & Sartori's paper
y = matrix(c(1, 8, 14), nrow=3, ncol=1, byrow=FALSE)
t = c(27)
u = c(80)

# Likelihood function
L = function(psi, beta, gamma, y) {
  L1 = ((gamma[1]*psi + beta[1])^y[1]) * exp(-(gamma[1]*psi + beta[1]))
  L2 = ((beta[1]*t[1])^y[2]) * exp(-beta[1]*t[1])
  L3 = ((gamma[1]*u[1])^y[3]) * exp(-gamma[1]*u[1])
  return(L1 * L2 * L3)
  #return(exp(l(psi, beta, gamma, y)))
}

# Log-likelihood function
l = function(psi, beta, gamma, y) {
  L1 = y[1] * log(gamma[1]*psi + beta[1]) - (gamma[1]*psi + beta[1])
  L2 = y[2] * log(beta[1]*t[1]) - beta[1]*t[1]
  L3 = y[3] * log(gamma[1]*u[1]) - gamma[1]*u[1]
  return(L1 + L2 + L3)
}

# Constrained MLE of gamma
gamma_p = function(psi, y) {
  K = psi * t[1] - u[1]

  A = K * (psi + u[1])
  B = (psi + u[1]) * (y[2] + y[3]) - K*y[1] - K*y[3]
  C = - y[3] * (y[1] + y[2] + y[3])

  numerator = -B + sqrt(B^2 - 4*A*C)
  return(numerator / (2*A))
}

# Constrained MLE of beta
beta_p = function(psi, y) {
  if (psi == 0) {
    return(0.321)
  }
  K = (psi * t[1]) - u[1]
  gamma = gamma_p(psi, y)

  numerator = psi * y[2] * gamma
  return(numerator / ((K*gamma) + y[3]))
}

beta_p_0 = function(y) {
  return( (y[1] + y[2])/(1+t) )
}
gamma_p_0 = function(y) {
  return( y[3]/u )
}

l_p = function(psi, y) {
  if (is.na(beta_p_0(y))) {
    print("UH OH")
    print(y)
  }
  val = l(psi, beta_p_0(y), gamma_p_0(y), y)
  return(val)
}

l_p(2, y) - l_p(5, y)



# Observed Fisher information
j = function(psi, beta, gamma, y) {
  theta_mle = getGlobalMLE(y)
  beta_mle = theta_mle[2]
  gamma_mle = theta_mle[3]
  y1_expected = gamma*psi + beta

  # remember: y[1] = gamma_mle*psi_mle + beta_mle
  dl_dpsi2 = -y[1]*(gamma^2)/(y1_expected^2)
  dl_dpsi_dbeta = -y[1]*gamma/(y1_expected^2)
  dl_dpsi_dgamma = -( (y[1]*psi*gamma/(y1_expected^2)) + 1)
    
  dl_dbeta2 = -( (y[1]/(y1_expected^2)) + (beta_mle*t/(beta^2)))
  dl_dbeta_dgamma = -y[1]*psi / (y1_expected^2)
    
  dl_dgamma2 = -( (y[1]*(psi^2)/(y1_expected^2)) + (gamma_mle*u/(gamma^2)))
  
  return(matrix(c(
      c(-dl_dpsi2,       -dl_dpsi_dbeta,   -dl_dpsi_dgamma),
      c(-dl_dpsi_dbeta,  -dl_dbeta2,       -dl_dbeta_dgamma),
      c(-dl_dpsi_dgamma, -dl_dbeta_dgamma, -dl_dgamma2)
    ), nrow=3, ncol=3)
  )
}

j(4.029, 0.321, 0.175, y)

# derivatives of the log-likelihood w.r.t. psi, beta, and gamma
dl = function(psi, beta, gamma, y) {
  y1_expected = gamma*psi + beta
  dl_dp = (y[1]*gamma/(y1_expected)) 
  dl_db = (y[1]/y1_expected) + (y[2]/beta) - (1+t)
  dl_dg = (y[1]*psi/(y1_expected)) + (y[3]/gamma) - (psi + u)
  return(c(dl_dp, dl_db, dl_dg))
}
dl(4.029, 0.321, 0.175, y)

# derivatives of the log-likelihood w.r.t. the MLEs
dl_dtheta_hat = function(psi, beta, gamma, y) {
  theta_mle = getGlobalMLE(y)
  psi_mle = theta_mle[1]
  gamma_mle = theta_mle[3]
  y1_expected = gamma*psi + beta
  
  y1_expected = gamma*psi + beta
  dl_dp = gamma_mle*log(y1_expected)
  dl_db = log(y1_expected) + t*log(beta*t)
  dl_dg = psi_mle*log(y1_expected) + u*log(gamma*u)
  return(c(dl_dp, dl_db, dl_dg))
}
dl_dtheta_hat(4.029, 0.296, 0.175, y)

# second derivatives of the log-likelihood w.r.t. 
# first: the MLEs and second: the nuisance parameters 
# (i.e. ∂l/∂λdθ_MLE )
dl_dlambda_dtheta_hat = function(psi, beta, gamma, y) {
  theta_mle = getGlobalMLE(y)
  psi_mle = theta_mle[1]
  beta_mle = theta_mle[2]
  gamma_mle = theta_mle[3]
  
  y1_expected = gamma*psi + beta
  dl_dbeta_dpsi_hat = gamma_mle/y1_expected
  dl_dbeta_dbeta_hat = (t/beta) + (1/y1_expected) 
  dl_dbeta_dgamma_hat = (psi_mle/y1_expected) 
  
  dl_dgamma_dpsi_hat = (gamma_mle*psi)/(y1_expected)
  dl_dgamma_dbeta_hat = psi/y1_expected
  dl_dgamma_dgamma_hat = (u/gamma) + ((psi_mle^2)/y1_expected)
  
  return(
    matrix(data=c(
            c(dl_dbeta_dpsi_hat, dl_dbeta_dbeta_hat, dl_dbeta_dgamma_hat),
            c(dl_dgamma_dpsi_hat, dl_dgamma_dbeta_hat, dl_dgamma_dgamma_hat)
          ), nrow=2, ncol=3)
  )
}
dl_dlambda_dtheta_hat(0, 0.296, 0.175, y)

# Barndorff-Nielsen's p* formula for the density of the MLE
p_star = function(psi_mle, psi_true, y) {
  beta_p = beta_p_0(y)
  gamma_p = gamma_p_0(y)
  
  # Corner of the Observed Information corresponding to
  # the nuisance parameters
  j_lambda = j(psi_mle, beta_p, gamma_p, y)[2:3, 2:3]
  j_term = det(j(psi_mle, beta_p, gamma_p, y)) / det(j_lambda)

  return(sqrt(j_term) * exp(-l_p(psi_true, y) + l_p(psi_mle, y)))
}

###################################################
####### Numerical solution helper functions
###################################################

# a wrapper around the likelihood to be used by optim()
optim_likelihood = function(param) {
  -l(param[1], param[2], param[3], y)
}

# a wrapper around l_p to be used by optim()
optim_l_p = function(psi, y) {
  -l_p(psi, y)
}

# Obtains the global MLE
getGlobalMLE = function(y) {
  beta = y[2]/t
  gamma = y[3]/u
  psi = (y[1] - beta) / gamma
  return(c(psi, beta, gamma))
}
