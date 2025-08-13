# Running this code should produce a vignette that fits the SIRS model
# to simulations using the ensemble MCMC method from Goodman and Weare (2010).

# See the readme for more detail, but note that model_functions.R involves
# a compilation of the differential equations in C. The dyn.load command therein
# is currently set to work for Linux/macOS. For Windows, users will need
# to modify the dyn.load command from dyn.load(sir.so) to dyn.load(sir.dll) 
# (or something like that) 

# By default, running this file should invoke the ensemble MCMC fitting procedure for 100
#	iterations (and these should print in the console)

########
## Load in all of the things we need:
#######

#load the ODE solver package:
require(deSolve)

# load in the model and the log-likelihood:
source('model_functions.R')

# load in the ensemble MCMC sampler (the function getEnsemble):
source('ensemble_sampler.R')

# load in some functions that are convenient for plotting the profile likelihoods:
source('plotfunctions.R')

# load in some default parameters:
source('default_parameters.R')


######
##  try running some stuff:
#####

## Make a simulated data set:
simdat <- getSimDat(theta_init)

## Calculate exected test results for comparison:
simmean <- getSimMean(theta_init)

## Calculate trajectories in the overall population:
simprev <- getSimPrev(theta_init)

## for plotting purposes, need to add a time column into simdat:
simdat_plot = cbind(simprev$t, simdat)
simmean_plot = cbind(simprev$t, simmean)
colnames(simdat_plot) <- c('t','y','x')
colnames(simmean_plot) <- c('t','y','x')


par(mfrow=c(1,2))
# Plot the simulated data:
plot(x~t,simdat_plot, xlab='Time [wks]', ylab='No. of detections',
	main='Simulation Data')
points(y~t,simdat_plot, col='red')
lines(x~t,simmean_plot,col='black', lwd=3)
lines(y~t,simmean_plot,col='black', lwd=3)
legend('topleft', col=c('black','red','black'),
	lty=c(NA,NA,1),	pch=c(19,19,NA),
	legend=c('Negative','Positive','Mean'))

# plot the trajectories in the overall population:
plot(S~t,simprev,type='l',col='blue',lwd=2,
	ylim=c(0, 327e6), ylab='No. of Population',
	main='Population Trajectories',xlab='Time [wks]')
lines(R~t,simprev,col='darkgreen',lwd=2)
lines(I~t,simprev,col='darkgrey',lwd=2)
legend('topright',legend=c('S','I','R'),lwd=2,
	col=c('blue','darkgrey','darkgreen'),
	lty=1)

######
##  try fitting the model:
#####


## specify a log-likelihood function
loglik <- function(theta) getlogl(simdat,theta)

## try fitting the model to the simulated data:
# running for dimTime=100 from the true values appears sufficient to resolve the likelihood surface somewhat:
# (try changing the number of timesteps (dimTime) and the ensemble size (dimEnsemble)): 

## run the ensemble MCMC:
ens <- getEnsemble(loglik, theta_init, dimTime = 100, dimEnsemble = 50)

## plot the profile likelihoods:
# we plot epidemiological parameters in a separate call from the testing model parameters

## plot Epidemiological parameter estimates:
plotEpiPar(ens)

## plot testing model parameter estimates:
plotPsiPar(ens)

print('True parameter values:')
print(c('R0'=R01, 'epsilon' = Epsilon, 'tau_beta' = (phase1+26)/52, 
	'gamma_inv' = 1/Gam1, 'rho_inv' = 1/Rho1, 'phi'=Phi1))


#save(ens, file='ensemble_mcmc_output.Rdata')
