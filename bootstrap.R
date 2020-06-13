# rv is the random variable whose critical value at 
# the alpha level we desire
estimateCritValue = function(rv, num_boot_samples, alpha, y) {
  psi_mle = getGlobalMLE(y)[1]
  rvs = rv(num_boot_samples, psi_mle, y)
  rvs = sort(rvs)
  rv_ecdf = 1:length(rvs)/length(rvs)
  crit = rvs[which(rv_ecdf >= 1-alpha)[1]]
  return(crit)
}

#estimateCritValue(rStrangeDist, 1, 0.1, y)
#psi_mle = getGlobalMLE(y)[1]
#rvs = rStrangeDist(10, psi_mle, y)
#rvs = sort(rvs)
#rv_ecdf = 1:length(rvs)/length(rvs)
#crit = rvs[which(rv_ecdf >= 0.9)[1]]

getBootstrapPValue = function(psi_0, B, r_obs) {
  count = 0
  total = 0
  for (i in 1:ncol(B)) {
    boot = B[,i]
    r_boot = r_psi(psi_0, boot)
    if (is.nan(r_boot) || is.na(r_boot) || 
        is.nan(r_obs) || is.na(r_obs)) {
      #next
      print("BOOT VALS")
      #browser()
      value = r_psi(psi_0, boot)
      print(r_boot)
      print(r_obs)
    }
    total = total + 1
    
    if (!is.nan(r_boot) &&
        !is.na(r_boot) &&
        r_boot >= r_obs) {
      count = count + 1
    }
  }
  #r_vals = sort(r_vals)
  #r_ecdf = 1:length(r_vals)/length(r_vals)
  #crit = r_vals[which(r_ecdf >= 1-alpha)[1]]
  return(count/total)
}
