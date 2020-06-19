### NOTE: Change this variable to the local directory that contains
###       these files
sourceDir = '~/Developer/inference_thesis/'

source(paste(sourceDir, 'poisson_hep_model.R', sep=""))
source(paste(sourceDir, 'test_statistics.R', sep=""))
source(paste(sourceDir, 'power_fns.R', sep=""))
source(paste(sourceDir, 'plotting.R', sep=""))

# initial values as seen in Davison & Sartori's paper
y = matrix(c(1, 8, 14), nrow=3, ncol=1, byrow=FALSE)
t = c(27)
u = c(80)

# Values for plots
psis = seq(0, 10, 2)
psi_0 = 0
beta = 0.296
gamma = 0.175
alpha = 0.10
ylim = c(0, 0.8)

num_reps = 2000

psi_string = expression(psi)
title = paste("Power fns; beta=", beta, 
              ", gamma=", gamma, 
              ", alph=", alpha, sep="")

plotPowerFns(psis, psi_0, beta, gamma, alpha,
             num_reps, num_boot_samples, origData=y,
             ylim, title)

plot(psis, power(psis, psi_0, beta, gamma, lrt, num_reps=num_reps, alpha), 
     type="l", col="blue", xlab=psi_string, ylab="Power", 
     main=title, ylim=ylim)

lines(psis, power(psis, psi_0, beta, gamma, r_star_psi, num_reps=num_reps, alpha), 
      col="dark green")

legend(x=0, y=0.8,
       legend=c("LRT", "r", "r.boot", "r*", "r* boot"),
       col=c("blue", "red", "purple", "dark green", 
             "dark orange"),
       lty=1:1, cex=0.8)

###################################
##### Plot that's closer to zero
###################################


psis = seq(0, 0.04, 0.002)
psi_0 = 0
beta = 2000
gamma = 2000
alpha = 0.10
ylim = c(0, 1)

num_reps = 2000
num_boot_samples = 100

psi_string = expression(psi)
title = paste("Power fns; beta=", beta, 
              ", gamma=", gamma, 
              ", alph=", alpha, sep="")

plotPowerFns(psis, psi_0, beta, gamma, alpha,
             num_reps, num_boot_samples, origData=y,
             ylim, title)

legend(x=0, y=1.0,
       legend=c("LRT", "r", "r.boot", "r*", "r* boot"),
       col=c("blue", "red", "purple", "dark green", 
             "dark orange"),
       lty=1:1, cex=0.8)
