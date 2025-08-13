# set initial parameters. Specifing parameters for 2-virus models (but can use them
#  for 1-virus models as well)
R01 = 3
R02 = 2.5
Gam1 = 1
Gam2 = 1
Rho1 = 1/52
Rho2 = 1/52
Psi0 = plogis(-16)
Phi1 = plogis(-15)
Phi2 = plogis(-15)
Sigma = 0.8
Chi = 1
Epsilon = 0.15
phase1=4
phase2=52*plogis(0)-26
phasepsi=52*plogis(0)+18
epspsi=0.75
deltapsi=1e-8

# assign parameter values for the single-pathogen SIRS model
R0 <- R01
Gam <- Gam1
Rho <- Rho1
Phi <- Phi1
tau_beta <- phase1

# place parameters in convenient vectors
# we transform parameters to R^N for MCMC/optimizing loglikelihood:
theta_coinf <- c(
		log(R01),
		log(R02),
		log(Gam1),
		log(Gam2),
		log(Rho1),
		log(Rho2),
		qlogis(Psi0),
		qlogis(Phi1),
		qlogis(Phi2),
		log(Sigma),
		log(Chi),
		qlogis(Epsilon),
		qlogis(.5),
		qlogis(.5),
		qlogis(.5),
		qlogis(.75))

theta_init <- c(
		log(R01),
		log(Gam1),
		log(Rho1),
		qlogis(Psi0),
		qlogis(Phi1),
		qlogis(Epsilon),
		qlogis((phase1+26)/52),
		qlogis((phasepsi-18)/52),
		qlogis(epspsi),
		qlogis(deltapsi))


