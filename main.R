### NOTE: Change this variable to the local directory that contains
###       these files
sourceDir = '~/Developer/inference_thesis/'

source(paste(sourceDir, 'poisson_hep_model.R', sep=""))
source(paste(sourceDir, 'test_statistics.R', sep=""))
source(paste(sourceDir, 'power_fns.R', sep=""))
source(paste(sourceDir, 'plotting.R', sep=""))
source(paste(sourceDir, 'bootstrap.R', sep=""))

# Example of current issue. I thought the observed information
# was positive determinant and thus j should have a
# positive determinant...
t = c(27, 27)
u = c(80, 80)

#MLEs = getGlobalMLE(c(4, 8, 10))
#j1 = j(MLEs[1], MLEs[2], MLEs[3], y=c(4, 8, 10))
#j1
#det(j1)

#differentMLEs = getGlobalMLE(c(2, 4, 15))
#j2 = j(differentMLEs[1], differentMLEs[2], differentMLEs[3], c(2, 4, 15))
#det(j2)



########################################
##### Some plots of power functions
########################################

# Values for plots
psis = seq(0, 10, 2)
psi_0 = 0
beta = 0.296
gamma = 0.175
alpha = 0.10
ylim = c(0, 1)

sampleSize = 3
num_reps = 1000
num_boot_samples = num_reps/12
  
psi_string = expression(psi)
title = paste("Power fns; beta=", beta, 
              ", gamma=", gamma, 
              ", alph=", alpha, sep="")
Y = matrix(nrow = 1, ncol = 3)
Y[1,] = c(1, 8, 14)

plotPowerFns(psis, psi_0, beta, gamma, alpha, sampleSize,
             num_reps, num_boot_samples, origData=Y,
             ylim, title)

legend(x=0, y=0.8,
       legend=c("LRT", "r", "r*", "r - standardized"),
       col=c("blue", "red", "dark green", "turquoise"),
       lty=1:1, cex=0.8)




###################################
##### Plot that's closer to zero
###################################


psis = seq(0, 1, 0.2)
psi_0 = 0
beta = 0.296
gamma = 0.175
alpha = 0.10
ylim = c(0, 0.5)

sampleSize = 5
num_reps = 500
num_boot_samples = 20

psi_string = expression(psi)
title = paste("Power fns; beta=", beta, 
              ", gamma=", gamma, 
              ", alph=", alpha, sep="")

plotPowerFns(psis, psi_0, beta, gamma, alpha, sampleSize,
             num_reps, num_boot_samples, origData=Y,
             ylim, title)

legend(x=0, y=0.5,
       legend=c("LRT", "r", "r*", "r - standardized"),
       col=c("blue", "red", "dark green", "turquoise"),
       lty=1:1, cex=0.8)
