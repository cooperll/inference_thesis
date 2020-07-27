plotPowerFns = function(psis, psi_0, beta, gamma, alpha, 
                        sampleSize, num_reps, num_boot_samples, 
                        origData, ylim, title) {
  plot(psis, power(psis, psi_0, beta, gamma, lrt, L_p,
                   num_reps=num_reps, alpha, sampleSize, isOneSided=TRUE), 
       type="l", col="blue", xlab=psi_string, ylab="Power", 
       main=title, ylim=ylim)
  
  lines(psis, power(psis, psi_0, beta, gamma, r_star, L_p,
                    num_reps=num_reps, alpha, sampleSize),
        col="dark green")

  lines(psis, power(psis, psi_0, beta, gamma, r,  L_p,
                    num_reps=num_reps, alpha, sampleSize), 
        col="red")
  
  lines(psis, power(psis, psi_0, beta, gamma, r_std, L_p,
                    num_reps=num_reps, alpha, sampleSize), 
        col="turquoise")
  
  #lines(psis, power(psis, psi_0, beta, gamma, r_star_std, L_p,
  #                  num_reps=num_reps, alpha, sampleSize),
  #      col="black")
  
  #lines(psis, power(psis, psi_0, beta, gamma, r_psi_boot, num_reps=num_reps/15, alpha, 
  #                  sampleSize, isBoot=TRUE, 
  #                  origData=origData, num_boot_samples=num_boot_samples), 
  #      col="purple")
  
  #lines(psis, power(psis, psi_0, beta, gamma, r_star_boot, num_reps=num_reps/15, alpha, 
  #                  isBoot=TRUE, 
  #                  sampleSize, origData=origData, num_boot_samples=num_boot_samples), 
  #      col="dark orange")
}
