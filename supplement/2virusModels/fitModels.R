rm(list=ls())

load('2virussims/df1.Rdata')
load('2virussims/df2.Rdata')
load('2virussims/df3.Rdata')

load('2virussims/simdat1.Rdata')
load('2virussims/simdat2.Rdata')
load('2virussims/simdat3.Rdata')


## load in the single-virus SIRS model definitiions:
source('model_functions.R')
source('plotfunctions.R')

# here are the parameter vectors close to the ones we estimated in the paper:
# feel free to try fitting from different guesses:
theta1hat <- c(log(1.36), -log(1.97), -log(14.94), -15.58, -16.56, qlogis(0.33), qlogis(0.3), qlogis(14/52), qlogis(0.51), -16.79)
theta2hat <- c(log(1.73), -log(1.54), -log(25.14), -15.60, -14.65, qlogis(0.66), qlogis(0.3), qlogis(14/52), qlogis(0.50), -16.78) 
theta3hat <- c(log(1.46), -log(1.38), -log(16.52), -15.57, -15.24, qlogis(0.58), qlogis(0.3), qlogis(14/52), qlogis(0.51), -16.81) 

loglik1 <- function(theta) getlogl(df1,theta)
loglik2 <- function(theta) getlogl(df2,theta)
loglik3 <- function(theta) getlogl(df3,theta)

# uncomment this to try fitting these:
source('ensemble_sampler.R')
ens1 <- getEnsemble(loglik1, theta1hat, 100, 30)
ens2 <- getEnsemble(loglik2, theta2hat, 100, 30)
ens3 <- getEnsemble(loglik3, theta3hat, 100, 30)

thetahat1 <- getmcmcmle(ens1)
# check if it looks like things converged:
plotEpiPar(ens1)
plotPsiPar(ens1)

# try for the others:
thetahat2 <- getmcmcmle(ens2)
plotEpiPar(ens2)
plotPsiPar(ens2)

thetahat3 <- getmcmcmle(ens3)
plotEpiPar(ens3)
plotPsiPar(ens3)



par(mfrow=c(2,2),mar=c(4,4,2,2))
plot(simdat1[,2],ylim=c(0,200),ylab='Cases',xlab='Time [weeks]',type='l',
	main='2 different viruses')
lines(simdat1[,3],col='red')
lines(simdat1[,4],col='blue')
lines(simdat1[,1],col='purple')
plot(df1[,2],ylim=c(0,200),ylab='Cases',xlab='Time [weeks]',type='l',
	main='Aggregating detections')
lines(df1[,1])
plotFitMean(df1,thetahat1,add=T)


par(mfrow=c(2,2),mar=c(4,4,2,2))
plot(simdat2[,2],ylim=c(0,200),ylab='Cases',xlab='Time [weeks]',type='l',
	main='2 different viruses')
lines(simdat2[,3],col='red')
lines(simdat2[,4],col='blue')
lines(simdat2[,1],col='purple')
plot(df2[,2],ylim=c(0,200),ylab='Cases',xlab='Time [weeks]',type='l',
	main='Aggregating detections')
lines(df2[,1])
plotFitMean(df2,thetahat2,add=T)

par(mfrow=c(2,2),mar=c(4,4,2,2))
plot(simdat3[,2],ylim=c(0,200),ylab='Cases',xlab='Time [weeks]',type='l',
	main='2 different viruses')
lines(simdat3[,3],col='red')
lines(simdat3[,4],col='blue')
lines(simdat3[,1],col='purple')
plot(df3[,2],ylim=c(0,200),ylab='Cases',xlab='Time [weeks]',type='l',
	main='Aggregating detections')
lines(df3[,1])
plotFitMean(df3,thetahat3,add=T)




