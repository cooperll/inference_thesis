# Obtains the global MLE
getGlobalMLE = function(y) {
  beta = y[2]/t
  gamma = y[3]/u
  
  # TODO: This seems extremely dodgy, and I dislike having to do this.
  # Are other inferences even valid with this? Thankfully, gamma's
  # MLE being zero is a rare event thanks to u = 80, but 
  # I'm forced to leave it here until I find an analytic solution. 
  # Maybe a transformation..? 
  if (gamma == 0) {
    psi = y[1] - beta
  } else {
    psi = (y[1] - beta) / gamma
  }
  
  ### Make psi obey the non-negativity constraint. 
  if (psi < 0) {
    psi = 0
  }
  return(c(psi, beta, gamma))
}

# Generates num_reps samples using the HEP Poisson model
generatePoissonData = function(psi, num_reps, beta, gamma) {
  y1 = rpois(num_reps, lambda=gamma*psi + beta)
  y2 = rpois(num_reps, lambda=beta*t)
  y3 = rpois(num_reps, lambda=gamma*u)
  Y = matrix(data=t(c(y1, y2, y3)), nrow=3, ncol=length(y1), byrow=TRUE)
  return(Y)
}

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
  v = gamma[1]*psi + beta[1]
  w = beta[1]*t[1]
  u = gamma[1]*u[1]
  
  if (is.nan(v) || is.na(v) || v <= 0) {
    v = 0.0000001
  }
  if (is.nan(w) || is.na(w) || w <= 0) {
    w = 0.0000001
  }
  if (is.nan(u) || is.na(u) || u <= 0) {
    u = 0.0000001
  }
  
  L1 = y[1] * log(v) - (gamma[1]*psi + beta[1])
  L2 = y[2] * log(w) - beta[1]*t[1]
  L3 = y[3] * log(u) - gamma[1]*u[1]
  return(L1 + L2 + L3)
}

# Constrained MLE of beta
beta_p_0 = function(y) {
  return( (y[1] + y[2])/(1+t) )
}

# Constrained MLE of gamma
gamma_p_0 = function(y) {
  return( y[3]/u )
}

# Profile likelihood 
L_p = function(psi, y) {
  beta_p = beta_p_0(y)
  gamma_p = gamma_p_0(y)
  return(L(psi, beta_p, gamma_p, y))
}

# Derivative of the profile likelihood
dL_p = function(psi, y) {
  beta_p = beta_p_0(y)
  gamma_p = gamma_p_0(y)
  
  y_fact = factorial(y[1]) * factorial(y[2]) * factorial(y[3])
  Q = ((beta_p*t)^y[2]) * ((gamma_p*u)^y[3])
  if (y[1] == 0) {
    psi_term_1 = gamma_p
  } else {
    psi_term_1 = y[1] * gamma_p * ((gamma_p * psi + beta_p)^(y[1]-1))
  }
  psi_term_2 =  gamma_p * exp(-((gamma_p * psi + beta_p) + beta_p*t + gamma_p*u))
  return(y_fact * Q * psi_term_1 * psi_term_2)
}

# Profile log likelihood
l_p = function(psi, y) {
  if (is.na(beta_p_0(y))) {
    print("UH OH")
    print(y)
  }
  val = l(psi, beta_p_0(y), gamma_p_0(y), y)
  return(val)
}

# Observed Fisher information
j = function(psi, beta, gamma, y) {
  theta_mle = getGlobalMLE(y)
  psi_mle = theta_mle[1]
  beta_mle = theta_mle[2]
  gamma_mle = theta_mle[3]
  
  y1_expected = gamma * psi + beta

  # I've already added the negative sign for j into 
  # all of the 2nd derivatives here
  dl_dpsi2 = y[1] * (gamma^2)/(y1_expected^2)
  dl_dpsi_dbeta = y[1] * gamma/(y1_expected^2)
  dl_dpsi_dgamma = (y[1] * psi * gamma /(y1_expected^2)) + 1
    
  dl_dbeta2 = (y[1]/(y1_expected^2)) + (y[2]/(beta^2))
  dl_dbeta_dgamma = y[1]  * psi / (y1_expected^2)
    
  dl_dgamma2 = (y[1]*(psi^2)/(y1_expected^2)) + (y[3]/(gamma^2))
  
  return(rbind(
      c(dl_dpsi2,       dl_dpsi_dbeta,   dl_dpsi_dgamma),
      c(dl_dpsi_dbeta,  dl_dbeta2,       dl_dbeta_dgamma),
      c(dl_dpsi_dgamma, dl_dbeta_dgamma, dl_dgamma2)
    )
  )
}


# Example of issue. I thought the observed information
# was positive determinant and thus j should have a
# positive determinant...
j1 = j(21, 0.2962963, 0.125, ytest)
j1
det(j1)

getGlobalMLE(c(2, 4, 15))
j2 = j(9.8765432, 0.1481481, 0.1875, c(2, 4, 15))
det(j2)

# derivatives of the log-likelihood w.r.t. psi, beta, and gamma
dl = function(psi, beta, gamma, y) {
  y1_expected = gamma*psi + beta

  dl_dp = (y[1]*gamma/(y1_expected)) - gamma
  dl_db = (y[1]/y1_expected) + (y[2]/beta) - (1+t)
  dl_dg = (y[1]*psi/(y1_expected)) + (y[3]/gamma) - (psi + u)
  return(c(dl_dp, dl_db, dl_dg))
}

# derivatives of the log-likelihood w.r.t. the MLEs
dl_dtheta_hat = function(psi, beta, gamma, y) {
  theta_mle = getGlobalMLE(y)
  psi_mle = theta_mle[1]
  gamma_mle = theta_mle[3]
  y1_expected = gamma*psi + beta
  if (y1_expected <= 0) {
    y1_expected = 0.0000001
  }
  if (beta <= 0) {
    beta = 0.0000001
  }
  if (gamma <= 0) {
    gamma = 0.0000001
  }
  dl_dp = gamma_mle*log(y1_expected)
  dl_db = log(y1_expected) + t*log(beta*t)
  dl_dg = psi_mle*log(y1_expected) + u*log(gamma*u)
  return(c(dl_dp, dl_db, dl_dg))
}

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
  dl_dgamma_dgamma_hat = (u/gamma) + (psi*psi_mle/y1_expected)
  
  return(
    rbind(c(dl_dbeta_dpsi_hat, dl_dbeta_dbeta_hat, dl_dbeta_dgamma_hat),
          c(dl_dgamma_dpsi_hat, dl_dgamma_dbeta_hat, dl_dgamma_dgamma_hat))
  )
}

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