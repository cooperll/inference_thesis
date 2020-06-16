source('~/Developer/inference_thesis/poisson_hep_model.R')
source('~/Developer/inference_thesis/test_statistics.R')
source('~/Developer/inference_thesis/power_fns.R')

# initial values as seen in Davison & Sartori's paper
y = matrix(c(1, 8, 14), nrow=3, ncol=1, byrow=FALSE)
t = c(27)
u = c(80)

# Values for plots
psis = seq(0, 10, 2)
psi_0 = 0
beta = 0.296
gamma = 0.175
alpha = 0.1
ylim = c(0, 1.2)

num_reps = 1000

num_reps_smaller = 100
num_boot_samples = 100

psi_string = expression(psi)
title = paste("Power fns; beta=", beta, 
              ", gamma=", gamma, 
              ", alph=", alpha, sep="")

plot(psis, power(psis, psi_0, beta, gamma, lrt, num_reps, alpha), 
     type="l", col="blue", xlab=psi_string, ylab="Power", 
     main=title, ylim=ylim)

legend(x=0, y=1.2,
       legend=c("LRT", "r", "r.boot", "r*_v1", 
                "r*_v2", "r*_v1 boot", "r*_v2 boot"),
       col=c("blue", "red", "purple", "dark green", 
             "dark red", "dark orange", "turquoise"),
       lty=1:1, cex=0.8)

lines(psis, power(psis, psi_0, beta, gamma, r_psi, num_reps=num_reps, alpha), 
      col="red")

lines(psis, power(psis, psi_0, beta, gamma, r_psi_boot, 
                  num_reps_smaller, alpha, 
                  isBoot=TRUE, num_boot_samples, origData=y), 
      col="purple")

lines(psis, power(psis, psi_0, beta, gamma, r_star_psi, num_reps, alpha), 
      col="dark green")

lines(psis, power(psis, psi_0, beta, gamma, r_star_psi_2, num_reps, alpha), 
      col="dark red")

lines(psis, power(psis, psi_0, beta, gamma, r_star_boot, 
                  num_reps_smaller, alpha, 
                  isBoot=TRUE, num_boot_samples, origData=y), 
      col="dark orange")

lines(psis, power(psis, psi_0, beta, gamma, r_star_boot_2, 
                  num_reps_smaller, alpha, 
                  isBoot=TRUE, num_boot_samples, origData=y), 
      col="turquoise")
