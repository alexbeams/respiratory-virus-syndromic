# This code contains all of the functions that define the SIR model:
#	1. getSim: user supplies parameters and initial conditions, and this returns
#		a dataframe of the solution to the ODE
#	2. getSimDat: user supplies the same parameters and intial conditions as getSim,
#		as well as parameters for the testing model. Returns a dataframe of 
#		positive detections (y) and negative detections (x).
#		By default, this runs getSim for 8000 weeks to eliminate transients,
#		and returns the last 1399 days (b/c that's how much data we use for 
#		our main analysis). For our analysis, we are trying to fit the attractor
#		of the ODE model
#	3. getSimMean: similar to getSimDat, but returns expected values of positive
#		and negative detections (x, y respectively)
#	4. getSimPrev: similar to getSim, but runs for the same amount of time as
#		getSimDat/getSimMean to visualize the attractor of the ODEs
#	5. getlogl: defines a loglikelihood function. User supplies data and parameters.
#		The data (x) need to have the same form as the output of getSimDat. 
#		By default, runs from random intial conditions, and runs ODEs for 8000
#		weeks
#
#	We optimize parameters in R^N; all functions transform parameters back before
#		running calculations.

#load the C file and compile the ODEs 
system("R CMD SHLIB sirmod/sir.c")

#load the compiled object
dyn.load("sirmod/sir.so")
# ^^ this step is different between Linux/macOS and Windows. This should work
# on Linux/macOS as-is. For Windows, the compiled object has extension .dll (I think)

#ODE model; accepts parameter values, outputs a big dataframe with state variables' timeseries.
# Able to run with default values.
getSim <- function(
	r0 = 3, #baseline R0's for pathogens
	eps =  0.1, #amplitude of seasonal forcing term
	phase=0, 	#phase shift(s) of seasonal forcing term for each pathogen
	tInit = 0,	#simulation start time	
	tEnd = 8000,	#simulation end time 
	DT = 1/7,	#simulate for each day; rates are all per week
	wndw=52*6,	#how many weeks to include in the final output?
	Npop = 327e6,	#fixed population size
	gam = 1,	#recovery rates, per week
	rho = 1/52,	#rates of immune loss, per week
	Y=runif(3),	#initial condition
	T_pd=52,	#period of seasonal transmission forcing term
	beta = r0*gam #transmission rates, parametrized in terms of R0
){
	parms <- c(
		eps=eps,
		phase=phase,
		beta=beta,
		gam=gam,
		rho=rho,
		T_pd=T_pd,
		Npop=Npop)

	#set up the initial conditions
	Y <- Y/sum(Y) * Npop
	# we replace dI/dt with d/dt log(I) b/c numbers get close to 0:	
	Y[2] = log(Y[2])	
	#specify the time at which we will have solutions
	times <- seq(tInit, tEnd, by = DT/1)
	#simulate the ODE
	out <- ode(Y, times, func = "derivs", parms = parms,
		jacfunc = "jac", dllname = "sir",
		initfunc = "initmod", nout = 1, outnames = "Sum")
	#convert output to an R data frame
	out <- as.data.frame(out)
	out <- out[seq(1,dim(out)[1], by=1),]

	#name the columns of the dataframe
	names(out) <- c('t','S','I','R','Forcing')
	# transform log(I) back
	out$I = exp(out$I)
	dat=out
	return(dat)
}


#define the log-likelihood 
getlogl <- function(x,theta,Y=runif(3),tEnd=8000){

	Npop = 327e6
	r0  = exp(theta[1])
	gam = exp(theta[2])
	rho = exp(theta[3])
	psi0 = plogis(theta[4])
	phi = plogis(theta[5])
	eps=plogis(theta[6])
	phase=52*plogis(theta[7])-26
	phasepsi=52*plogis(theta[8])+18
	epspsi=plogis(theta[9])	
	psi1=plogis(theta[10])

	# initial data Y get normalized and transformed to log(I) inside getSim():	
	sim = getSim(r0=r0,gam=gam,rho=rho,eps=eps,Npop=Npop,phase=phase,Y=Y,tEnd=tEnd) 
	sim = tail(sim, n=dim(x)[1])
	
	#aggregate by infected and uninfected: 
	Q = cbind(c(sim$I),c(sim$S+sim$R))
	
	#independent testing probability (symptom-driven, but unrelated to the 
	#  focal pathogen of interest); seasonal with a linear increase:
	psi = (psi0+psi1*(sim$t-sim$t[1])/52 )  * (1 + epspsi * sin(2*pi*(sim$t-phasepsi)/52)) 

	# testing probabilities for each state:
	f1 = 1-(1-psi)*(1-phi)
	f0 = psi 
	
	#state-specific and time-dependent testing probabilities:
	f = data.frame( cbind(f1,f0)) 
	Q = as.matrix(Q)


	# positive detections, y, are realizations of Y ~ Poisson( (phi+psi)*I),
	# negative detections, x, are realizations of X ~ Poission( psi*(S+R) ):
	logl <- sum( 
			apply( 
				cbind(x,f*Q), 1, function(x) sum(dpois(x[1:2],lambda=x[3:4],log=T)) 
			)
		)

	return(logl)	
}


# produce a simulated dataset - calls getSim:

getSimDat <- function(theta,DT=1/7,wndw=52*6,tEnd=8000, Y=runif(3)){
	r0  = exp(theta[1])
	gam = exp(theta[2])
	rho = exp(theta[3])
	psi0 = plogis(theta[4])
	phi = plogis(theta[5])
	eps=plogis(theta[6])
	phase=52*plogis(theta[7])-26
	phasepsi=52*plogis(theta[8])+18
	epspsi=plogis(theta[9])	
	psi1=plogis(theta[10])

	Npop=327e6

	#simulate the ODEs 
	sim = getSim(r0=r0,gam=gam,rho=rho,
		eps=eps,phase=phase,DT=DT,wndw=wndw,tEnd=tEnd, 
		Y=Y)
	sim=tail(sim,n=1399)
	#aggregate in terms of the "contingency table" variables
	Q = cbind(c(sim$I),c(sim$S+sim$R))

	#testing probabilities for each state:
	psi = (psi0+psi1*(sim$t-sim$t[1])/52 )  * (1 + epspsi * sin(2*pi*(sim$t-phasepsi)/52)) 

	f1 = 1-(1-psi)*(1-phi)
	f0 = psi 
	
	#state-specific and time-dependent testing probabilities:
	f = data.frame( cbind(f1,f0) )
	Q = as.matrix(Q)

	#average number of tests each day:
	lambda = rowSums( f*Q )

	#sample the number of tests at each time:
	n = sapply(lambda,function(x) rpois(1,x))
	#n = sapply(lambda, function(x) rnbinom(1,mu=x,size=100) ) #no. of tests has more variance than poisson
	
	#allocate tests among different types of test results:
	x = t(apply(cbind(n,f*Q),1,function(y) rmultinom(1,unlist(y[1]),unlist(y[2:3])) ))
	x = as.data.frame(x)
	names(x) <- c('y','x')
	
	return(x)	
}

# similar to getSimDat, but just produces mean counts 

getSimMean <- function(theta,DT=1/7,wndw=52*6,tEnd=8000, Y=runif(3),nrows=1399,Npop=327e6){
	r0  = exp(theta[1])
	gam = exp(theta[2])
	rho = exp(theta[3])
	psi0 = plogis(theta[4])
	phi = plogis(theta[5])
	eps=plogis(theta[6])
	phase=52*plogis(theta[7])-26
	phasepsi=52*plogis(theta[8])+18
	epspsi=plogis(theta[9])	
	psi1=plogis(theta[10])

	#Npop=327e6
	#simulate the ODEs	
	sim = getSim(r0=r0,gam=gam,rho=rho,
		eps=eps,phase=phase,DT=DT,wndw=wndw,tEnd=tEnd, Y=Y,Npop=Npop)
	sim=tail(sim,n=nrows)
	#aggregate in terms of the "contingency table" variables
	Q = cbind(c(sim$I),c(sim$S+sim$R))

	#testing probabilities for each state:
	psi = (psi0+psi1*(sim$t-sim$t[1])/52 )  * (1 + epspsi * sin(2*pi*(sim$t-phasepsi)/52)) 

	f1 = 1-(1-psi)*(1-phi)
	f0 = psi 
	
	#state-specific and time-dependent testing probabilities:
	f = data.frame( cbind(f1,f0) )
	Q = as.matrix(Q)

	#average number of tests each day:
	lambda = rowSums( f*Q )

	x <- f*Q
	names(x) <- c('y','x')

	return(x)	
}

# similar to getSimDat and getSimMean, but produces the population trajecotories
#	that are not observed

getSimPrev <- function(theta,DT=1/7,wndw=52*6,tInit=0,tEnd=8000, Y=runif(3),nrows=1399){
	r0  = exp(theta[1])
	gam = exp(theta[2])
	rho = exp(theta[3])
	psi0 = plogis(theta[4])
	phi = plogis(theta[5])
	eps=plogis(theta[6])
	phase=52*plogis(theta[7])-26
	phasepsi=52*plogis(theta[8])+18
	epspsi=plogis(theta[9])	
	psi1=plogis(theta[10])

	Npop=327e6
	
	#simulate the ODEs - focal pathogen is the 1st, runs independent of the 2nd	
	sim = getSim(r0=r0,gam=gam,rho=rho,
		eps=eps,phase=phase,DT=DT,
		wndw=wndw,tInit=tInit,
		tEnd=tEnd, Y=Y)
	sim=tail(sim,n=nrows)
	return(sim)	
}


