#
# power() is a helper function - it is passed the relevant parameters,
# the TestStat (such as r(psi)) and the criticalValue, and it
# computes the power of the test statistic for the given value of psi
#
# If the test statistic is going to utilize bootstrap for estimation
# of p values, then the additional parameters origData and 
# num_boot_samples *cannot* be NULL. 
power = function(psis, psi_0, beta, gamma, TestStat,
                 num_reps, alpha, isBoot=FALSE, 
                 origData=NULL, num_boot_samples=NULL) {
  res = c()
  
  ### These are set only for the bootstrapped statistics.
  ### They get used by TestStat, the test statistic whose power
  ### function we desire. See test_statistics.R for individual
  ### test statistic implementations. 
  bootData = NULL
  if (isBoot) {
    beta_c = beta_p_0(origData)
    gamma_c = gamma_p_0(origData)
    bootData = generatePoissonData(psi_0, num_boot_samples, beta_c, gamma_c)
  }
  
  total = 0
  for (psi in psis) {
    Y = generatePoissonData(psi, num_reps, beta, gamma)
    num_rejected = 0
    for (i in 1:ncol(Y)) {
      y = Y[,i]
      testStatRes = TestStat(psi_0, y, bootData)

      ### TODO: remove after stability is guaranteed.
      ### Right now, just makes sure that power doesn't
      ### fail on certain cases, like when y_2 or y_3 are 0.
      if (is.nan(testStatRes[2])) {
        next
      }
      total = total + 1
      
      if (testStatRes[2] < alpha) {
        num_rejected = num_rejected + 1
      }
    }
    res = append(res, num_rejected/num_reps)
  }
  return(res)
}
