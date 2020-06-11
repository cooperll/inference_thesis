# Likelihood Ratio Test
lrt = function(psi_0, y) {
  beta_p = beta_p_0(y)
  gamma_p = gamma_p_0(y)
  psi_global_mle = getGlobalMLE(y)[1]

  diff = l(psi_global_mle, beta_p, gamma_p, y) - l(psi_0, beta_p, gamma_p, y)
  return(2 * diff)
}

# Signed root likelihood statistic
r_psi = function(psi_0, y) {
  psi_global_mle = getGlobalMLE(y)[1]
  diff = l_p(psi_global_mle, y) - l_p(psi_0, y)
  val = sign(psi_global_mle - psi_0) * sqrt(2 * diff)
  return(val)
}

# Barndorff-Nielsen's r* statistic
r_star_psi = function(psi_0, y) {
  MLEs = getGlobalMLE(y)
  beta_p = beta_p_0(y)
  gamma_p = gamma_p_0(y)
  
  r_psi = r_psi(psi_0, y)
  
  dl_dtheta_hat_MLE = dl_dtheta_hat(MLEs[1], MLEs[2], MLEs[3], y)
  dl_dtheta_hat_constrained = dl_dtheta_hat(psi_0, beta_p, gamma_p, y)
  m = matrix(data=c(
    dl_dtheta_hat_MLE - dl_dtheta_hat_constrained,
    dl_dlambda_dtheta_hat(psi_0, beta_p, gamma_p, y)
  ), nrow=3, ncol=3)
  det_dl_term = det(m)
  
  j_lambda = j(psi_0, beta_p, gamma_p, y)[2:3, 2:3]
  j_term = sqrt(det(j_lambda)) * sqrt(det(j(MLEs[1], MLEs[2], MLEs[3], y)))
  
  v_p = det_dl_term / j_term
  return(r_psi + (log(v_p/r_psi)) / r_psi)
}
1 - pnorm(r_psi(-1, y))

Y = seq(-5, 15, 1)
plot(Y, r_star_psi(0.5, Y), 
     type="l", col="blue",
     xlab=expression(psi), ylab="Value", 
     main="Verification of r*")


R = function(psi_0, y) {
  return(p_star(Y, Y, y)/p_star(Y, 0, y))
}

fc_psi = function(psi_0, y) {
  
}


psis = seq(-1.5, 15, 0.01)
plot(Y, p_star(Y, 0, y), 
     type="l", col="blue",
     xlab=expression(psi), ylab="Value", 
     main="p*")
lines(Y, p_star(Y, Y, y), 
      type="l", col="blue",
      xlab=expression(psi), ylab="Value", 
      main="p*")
lines(Y, p_star(Y, 4, y), 
      type="l", col="blue",
      xlab=expression(psi), ylab="Value", 
      main="p*")

plot(Y, p_star(Y, Y, y)/p_star(Y, 0, y), 
     type="l", col="blue",
     xlab=expression(psi), ylab="Value", 
     main="p*")
lines(Y, 1 - 4000*p_star(Y, 4, y), 
      type="l", col="red",
      xlab=expression(psi), ylab="Value", 
      main="p*")

p_star(1, 4, y)
4000*p_star(Y, 4, y)



psis = seq(-1.5, 15, 0.01)
res = c()
for (psi in psis) {
  res = append(res, p_star(psi, 0, y))
}
plot(psis, res, 
     type="l", col="blue",
     xlab=expression(psi), ylab="Value", 
     main="p*")
p_star(3, 5, y)
