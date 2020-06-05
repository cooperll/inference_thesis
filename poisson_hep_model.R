# initial values as seen in Davison & Sartori's paper
y = matrix(c(1, 8, 14), nrow=3, ncol=1, byrow=FALSE)
t = c(27)
u = c(80)

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

l_p = function(psi, y) {
  val = l(psi, beta_p(psi, y), gamma_p(psi, y), y)
  return(val)
}


##########################################
####### Numerical solution helper functions
##########################################

# a wrapper around the likelihood to be used by optim()
optim_likelihood = function(param) {
  -l(param[1], param[2], param[3], y)
}

# a wrapper around l_p to be used by optim()
optim_l_p = function(psi, y) {
  -l_p(psi, y)
}

# Obtains a global MLE by maximizing the
# profile likelihood
getGlobalMLE = function(psi_init, y) {
  res = optim(par=c(psi_init),
              fn=optim_l_p, gr=NULL, y,
              method="Brent", lower=-15, upper=150)
  theta_global_mle = res$par
  return(theta_global_mle)
}
