# Obtains the global MLE
# Input: Y - an N x 3 matrix of the N samples
getGlobalMLE = function(Y) {
  N = nrow(Y)
  
  y1_bar = mean(Y[,1])
  y2_bar = mean(Y[,2])
  y3_bar = mean(Y[,3])
  beta = y2_bar/t[1]
  gamma = y3_bar/u[1]

  # TODO: This seems extremely dodgy, and I dislike having to do this.
  # Are other inferences even valid with this? Thankfully, gamma's
  # MLE being zero is a rare event thanks to u = 80, but
  # I'm forced to leave it here until I find an analytic solution.
  # Maybe a transformation..?
  if (gamma == 0) {
    psi = y1_bar - beta
  } else {
    psi = (y1_bar - beta) / gamma
  }

  ### Make psi obey the non-negativity constraint.
  #if (psi < 0) {
  #  psi = 0
  #}
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
L = function(psi, beta, gamma, Y) {
  N = nrow(Y)
  likelihood = 1
  for (i in 1:N) {
    L1 = ((gamma[1]*psi + beta[1])^Y[i, 1]) * exp(-(gamma[1]*psi + beta[1]))
    L2 = ((beta[1]*t[1])^Y[i, 2]) * exp(-beta[1]*t[1])
    L3 = ((gamma[1]*u[1])^Y[i, 3]) * exp(-gamma[1]*u[1])
    factorial = factorial(Y[i, 1]) * factorial(Y[i, 2]) * factorial(Y[i, 3])
    likelihood = likelihood * (L1 * L2 * L3 / factorial)
  }
  return(likelihood)
  #return(exp(l(psi, beta, gamma, y)))
}

# Log-likelihood function
l = function(psi, beta, gamma, Y) {
  N = nrow(Y)
  likelihood = 0

  for (i in 1:N) {
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

    L1 = Y[i, 1] * log(v) - (gamma[1]*psi + beta[1])
    L2 = Y[i, 2] * log(w) - beta[1]*t[1]
    L3 = Y[i, 3] * log(u) - gamma[1]*u[1]
    likelihood = likelihood + (L1 + L2 + L3)
  }
  return(likelihood)
}

# Constrained MLEs of beta
beta_c = function(Y) {
  y1_bar = mean(Y[,1])
  y2_bar = mean(Y[,2])
  beta_constr = (y1_bar + y2_bar)/(1 + t[1])
  return(beta_constr)
}

# Constrained MLEs of gamma
gamma_c = function(Y) {
  y3_bar = mean(Y[,3])
  gamma_constr = y3_bar/u[1]
  return(gamma_constr)
}

# Profile likelihood
L_p = function(psi, Y) {
  beta_constr = beta_c(Y)
  gamma_constr = gamma_c(Y)
  return(L(psi, beta_constr, gamma_constr, Y))
}

# Derivative of the profile likelihood
dL_p = function(psi, Y) {
  beta_constr = beta_c(Y)
  gamma_constr = gamma_c(Y)

  y_fact = factorial(Y[1]) * factorial(Y[2]) * factorial(Y[3])
  Q = ((beta_constr*t)^Y[2]) * ((gamma_constr*u)^Y[3])
  if (Y[1] == 0) {
    psi_term_1 = gamma_constr
  } else {
    psi_term_1 = Y[1] * gamma_constr * ((gamma_constr * psi + beta_constr)^(Y[1]-1))
  }
  psi_term_2 =  gamma_constr * exp(-((gamma_constr * psi + beta_constr) + beta_constr*t + gamma_constr*u))
  return(y_fact * Q * psi_term_1 * psi_term_2)
}

# Profile log likelihood
l_p = function(psi, Y) {
  val = l(psi, beta_c(Y), gamma_c(Y), Y)
  return(val)
}

# Observed Fisher information
j = function(psi, beta, gamma, Y) {
  N = 1
  mles = getGlobalMLE(Y)
  psi_mle = mles[1]
  beta_mle = mles[2]
  gamma_mle = mles[3]

  P = gamma * psi + beta
  Phat = gamma_mle * psi_mle + beta_mle

  # I've already added the negative sign for j into
  # all of the 2nd derivatives here
  dl_dpsi2 = N * Phat * (gamma^2)/(P^2)
  dl_dpsi_dbeta = N * Phat * gamma/(P^2)
  dl_dpsi_dgamma = (N* Phat * psi * gamma /(P^2)) + 1

  dl_dbeta2 = (N * Phat/(P^2)) + (N*beta_mle*t[1]/(beta^2))
  dl_dbeta_dgamma = N * Phat  * psi / (P^2)

  dl_dgamma2 = (N * Phat*(psi^2)/(P^2)) + (N*gamma_mle*u[1]/(gamma^2))

  return(rbind(
      c(dl_dpsi2,       dl_dpsi_dbeta,   dl_dpsi_dgamma),
      c(dl_dpsi_dbeta,  dl_dbeta2,       dl_dbeta_dgamma),
      c(dl_dpsi_dgamma, dl_dbeta_dgamma, dl_dgamma2)
    )
  )
}

# derivatives of the log-likelihood w.r.t. psi, beta, and gamma
dl = function(psi, beta, gamma, Y) {
  y1_expected = gamma*psi + beta

  dl_dp = (Y[1]*gamma/(y1_expected)) - gamma
  dl_db = (Y[1]/y1_expected) + (Y[2]/beta) - (1+t[1])
  dl_dg = (Y[1]*psi/(y1_expected)) + (Y[3]/gamma) - (psi + u[1])
  return(c(dl_dp, dl_db, dl_dg))
}

# derivatives of the log-likelihood w.r.t. the MLEs
dl_dtheta_hat = function(psi, beta, gamma, Y) {
  N = 1
  mles = getGlobalMLE(Y)
  psi_mle = mles[1]
  gamma_mle = mles[3]
  P = gamma*psi + beta
  if (P <= 0) {
    P = 0.0000001
  }
  if (beta <= 0) {
    beta = 0.0000001
  }
  if (gamma <= 0) {
    gamma = 0.0000001
  }
  dl_dp = N*gamma_mle*log(P)
  dl_db = N*(log(P) + t[1]*log(beta*t[1]))
  dl_dg = N*(psi_mle*log(P) + u[1]*log(gamma*u[1]))
  return(c(dl_dp, dl_db, dl_dg))
}

# second derivatives of the log-likelihood w.r.t.
# first: the MLEs and second: the nuisance parameters
# (i.e. ∂l/∂λdθ_MLE )
dl_dlambda_dtheta_hat = function(psi, beta, gamma, Y) {
  N = 1
  mles = getGlobalMLE(Y)
  psi_mle = mles[1]
  beta_mle = mles[2]
  gamma_mle = mles[3]

  P = gamma*psi + beta
  
  dl_dbeta_dpsi_hat = N*gamma_mle/P
  dl_dbeta_dbeta_hat = N*((t[1]/beta) + (1/P))
  dl_dbeta_dgamma_hat = N*(psi_mle/P)

  dl_dgamma_dpsi_hat = N*(gamma_mle*psi)/(P)
  dl_dgamma_dbeta_hat = N*psi/P
  dl_dgamma_dgamma_hat = N*((u[1]/gamma) + (psi*psi_mle/P))

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
