plotPowerFns = function(psis, psi_0, beta, gamma, alpha, 
                        num_reps, num_boot_samples, 
                        origData, ylim, title) {
  plot(psis, power(psis, psi_0, beta, gamma, lrt, num_reps=num_reps, alpha), 
       type="l", col="blue", xlab=psi_string, ylab="Power", 
       main=title, ylim=ylim)
  
  lines(psis, power(psis, psi_0, beta, gamma, r_star_psi, num_reps=num_reps, alpha), 
        col="dark green")
  
  lines(psis, power(psis, psi_0, beta, gamma, r_psi, num_reps=num_reps, alpha), 
        col="red")
  
  lines(psis, power(psis, psi_0, beta, gamma, r_psi_boot, num_reps=num_reps/15, alpha, 
                    isBoot=TRUE, num_boot_samples=num_boot_samples, origData=y), 
        col="purple")
  
  lines(psis, power(psis, psi_0, beta, gamma, r_star_boot, num_reps=num_reps/15, alpha, 
                    isBoot=TRUE, num_boot_samples=num_boot_samples, origData=y), 
        col="dark orange")
}
