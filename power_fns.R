#
# power() is a helper function - it is passed the relevant parameters,
# the TestStat (such as r(psi)) and the criticalValue, and it
# computes the power of the test statistic for the given value of psi
#
# If the test statistic is going to utilize bootstrap for estimation
# of p values, then the additional parameters origData and 
# num_boot_samples *cannot* be NULL. 
power = function(psis, psi_0, beta, gamma, TestStat, L,
                 num_reps, alpha, sampleSize, isOneSided=FALSE, 
                 isBoot=FALSE, origData=NULL, num_boot_samples=NULL) {
  res = c()
  
  ### These are set only for the bootstrapped statistics.
  ### They get used by TestStat, the test statistic whose power
  ### function we desire. See test_statistics.R for individual
  ### test statistic implementations. 
  bootData = NULL
  if (isBoot) {
    beta_constr = beta_c(origData)
    gamma_constr = gamma_c(origData)
    bootData = generatePoissonData(psi_0, num_boot_samples, 
                                   beta_constr, gamma_constr)
  }
  
  for (psi in psis) {
    data = t(generatePoissonData(psi, sampleSize*num_reps, beta, gamma))
    num_rejected = 0
    total = 0
    
    for (i in 1:num_reps) {
      Y = matrix(nrow=sampleSize, ncol=ncol(data))
      idx = ((i-1)*sampleSize) + 1
      
      if (sampleSize == 1) {
        Y[1,] = data[idx,]
      } else {
        Y = data[idx:(idx+(sampleSize-1)),]
      }
      
      psi_mle = getGlobalMLE(Y)[1]
      testStatRes = TestStat(psi_0, psi_mle, L, Y, bootData)

      ### TODO: remove after stability is guaranteed.
      ### Right now, just makes sure that power doesn't
      ### fail on certain cases, like when y_2 or y_3 are 0.
      if (is.nan(testStatRes[2]) || is.na(testStatRes[2])) {
        test = TestStat(psi_0, psi_mle, L, Y, bootData)
        #print("SKIP")
        next
      }
      total = total + 1
      
      if (isOneSided) {
        if (testStatRes[2] < alpha) {
          num_rejected = num_rejected + 1
        }
      } else {
        if (testStatRes[2] < alpha/2 || testStatRes[2] > 1-(alpha/2)) {
          num_rejected = num_rejected + 1
        }
      }
    }
    res = append(res, num_rejected/num_reps)
  }
  return(res)
}
