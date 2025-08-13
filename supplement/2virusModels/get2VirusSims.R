source('getModelEstimates.R')

require(xtable)
# set parameters for the 2-virus models:
# This generates the 2-virus simulations for Figures S25-27 of the Supplement.
# Use this code to generate the simulations from 2-virus models, and then try fitting single-virus SIRS models
#  to the aggregated detections

# 2virusmod1:
theta1 = c(log(2.15), log(2.13), -1, -1, -log(104), -log(104), -15.58, -16, -16, 0, 0, -3.46e-1, -4.42, -4.13, -1, 2.94e-2, -16.80)

# 2virusmod2
theta2 = c(log(2.5), log(2.5), -0.3, -0.3, -log(52), -log(52), -15.58, -15, -15, 0, 0, -3.46e-1, -4.42, -4.13, -1, 2.94e-2, -16.80)

# 2virusmod3
theta3 = c(log(2.5), log(1.5), -0.3, -0.3, -log(52), -log(52), -15.58, -15, -15, 0, 0, -3.46e-1, -4.42, -4.13, -1, 2.94e-2, -16.80)

# Initial condition that we use in the paper:
Y <- c(0.2134986, 0.006026265, 0.0001458091, 4.115645e-06, 0, 0.1751144,
       0.0002251385, 0.004942823, 0.2703876)

simdat1 <- getSimDat(theta1, Y=Y)
simdat2 <- getSimDat(theta2, Y=Y)
simdat3 <- getSimDat(theta3, Y=Y)

df1 = simdat1
df1$x1 = df1$x11+df1$x10+df1$x01
df1$x0 = df1$x00
df1 = df1[,c('x1','x0')]

df2 = simdat2
df2$x1 = df2$x11+df2$x10+df2$x01
df2$x0 = df2$x00
df2 = df2[,c('x1','x0')]

df3 = simdat3
df3$x1 = df3$x11+df3$x10+df3$x01
df3$x0 = df3$x00
df3 = df3[,c('x1','x0')]



# make a convenient table of the parameter values:
modpars = rbind(theta1,theta2,theta3)
rownames(modpars) = paste0('Model',1:3)
modpars[,1:2] = exp(modpars[,1:2])
modpars[,3:6] = exp(-modpars[,3:6])
modpars[,12] = plogis(modpars[,12])
modpars[,13] = 52*plogis(modpars[,13])-26
modpars[,14] = 52*plogis(modpars[,14])-26
modpars[,15] = 52*plogis(modpars[,15])+18
modpars[,16] = plogis(modpars[,16])
modpars = modpars[,-c(10,11)]
colnames(modpars) <- c('r01','r02','gam1','gam2','rho1','rho2','psi0','phi1','phi2','eps','phase1','phase2',
	'phasepsi','epspsi','psi1')
modpars = modpars[,c(1,2,3,4,5,6,10,11,12,8,9,7,15,13,14)]
xtable(modpars)

# Initial condition that we use in the paper:
Y <- c(0.2134986, 0.006026265, 0.0001458091, 4.115645e-06, 0, 0.1751144,
       0.0002251385, 0.004942823, 0.2703876)

simdat1 <- getSimDat(theta1, Y=Y)
simdat2 <- getSimDat(theta2, Y=Y)
simdat3 <- getSimDat(theta3, Y=Y)

df1 = simdat1
df1$x1 = df1$x11+df1$x10+df1$x01
df1$x0 = df1$x00
df1 = df1[,c('x1','x0')]

df2 = simdat2
df2$x1 = df2$x11+df2$x10+df2$x01
df2$x0 = df2$x00
df2 = df2[,c('x1','x0')]

df3 = simdat3
df3$x1 = df3$x11+df3$x10+df3$x01
df3$x0 = df3$x00
df3 = df3[,c('x1','x0')]

save(simdat1, file='2virussims/simdat1.Rdata')
save(simdat2, file='2virussims/simdat2.Rdata')
save(simdat3, file='2virussims/simdat3.Rdata')
save(df1, file='2virussims/df1.Rdata')
save(df2, file='2virussims/df2.Rdata')
save(df3, file='2virussims/df3.Rdata')

