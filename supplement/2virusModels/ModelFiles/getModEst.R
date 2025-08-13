
#load the ODE solver
require(deSolve)


#load the C file and compiler to run the ODEs at lightspeed
system("R CMD SHLIB ModelFiles/sirx2.c")
dyn.load("ModelFiles/sirx2.so")
# ^^ change this to sirx2.dll if you are Windows (I think)


# set initial parameter values
#R01 = 3
#R02 = 2.5
#Gam1 = 1
#Gam2 = 1
#Rho1 = 1/52
#Rho2 = 1/52
#Psi0 = plogis(-16)
#Phi1 = plogis(-15)
#Phi2 = plogis(-15)
#Sigma = 0.8
#Chi = 1
#Epsilon = 0.15
#phase1=52*plogis(0)-26
#phase2=52*plogis(0)-26
#phasepsi=52*plogis(0)+18
#epspsi=0.75

#transform parameters to real numbers for the MCMC
#theta_init <- c(
#		log(R01),
#		log(R02),
#		log(Gam1),
#		log(Gam2),
#		log(Rho1),
#		log(Rho2),
#		qlogis(Psi0),
#		qlogis(Phi1),
#		qlogis(Phi2),
#		log(Sigma),
#		log(Chi),
#		qlogis(Epsilon),
#		qlogis(.5),
#		qlogis(.5),
#		qlogis(.5),
#		qlogis(.75))


#ODE model; accepts parameter values, outputs a big dataframe with state variables' timeseries
getSim <- function(
	r01 = 3, #baseline R0's for pathogens
	r02 = 2,
	eps =  0.1, #amplitude of seasonal forcing term
	phase1=0, 	#phase shift(s) of seasonal forcing term for each pathogen
	phase2=0,
	tInit = 0,	#simulation start time	
	tEnd = 8000,	#simulation end time 
	DT = 1/7,	#simulate for each day; rates are all per week
	wndw=52*6,	#how many weeks to include in the final output?
	Npop = 327e6,	#fixed population size
	chi=1,	#symmetric cross-immunity term
	chi1 = chi, #asymmetric cross-immunity from virus 1 against virus 2
	chi2 = chi, #asymmetric cross-immunity from virus 2 against virus 1
	gam1 = 1,	#recovery rates, per week
	gam2 = 1,
	rho1 = 1/52,	#rates of immune loss, per week
	rho2 = 1/52,
	Y=runif(9),	#initial condition
	T_pd=52,	#period of seasonal transmission forcing term
	eta = 1, 	#altered coinfection duration
	gam1p = eta*gam1,
	gam2p = eta*gam2,
	nu = 1,		#altered infection duration from cross-immunity
	gam1tilde = nu*gam1,
	gam2tilde = nu*gam2,
	beta1 = r01*gam1, #transmission rates, parametrized in terms of R0
	beta2 = r02*gam2,
	sigma=1,		#altered suscepibility to coinfection
	sigma1=sigma,		#altered susceptibility to virus 1 if infected by virus 2
	beta1p = sigma1*beta1,
	beta1tilde = chi2*beta1,
	sigma2=sigma,		#altered susceptibility to virus 2 if infected by virus 1
	beta2p = sigma2*beta2,
	beta2tilde = chi1*beta2,
	rho1p=rho1,	#rate at which immunity to pathogen 1 is lost, given infection with pathogen 2
	rho1tilde=rho1,	#rate at which immunity to pathogen 1 is lost, given immunity to pathogen 2
	rho2p=rho2,	#rate at which immunity to pathogen 2 is lost, given infection with pathogen 1
	rho2tilde=rho2	#rate at which immunity to pathogen 2 is lost, given immunity to pathogen 1
){
	parms <- c(
		eps=eps,
		phase1=phase1,
		phase2=phase2,
		beta1=beta1,
		beta1p=beta1p,
		beta1tilde=beta1tilde,
		beta2=beta2,
		beta2p=beta2p,
		beta2tilde=beta2tilde,
		gam1=gam1,
		gam2=gam2,
		gam1p=gam1p,
		gam2p=gam2p,
		gam1tilde=gam1tilde,
		gam2tilde=gam2tilde,
		rho1=rho1,
		rho2=rho2,
		rho1p=rho1p,
		rho2p=rho2p,
		rho1tilde=rho1tilde,
		rho2tilde=rho2tilde,
		T_pd=T_pd,
		Npop=Npop)

	#set up the initial conditions
	Y <- Y/sum(Y) * Npop
	#specify the time at which we will have solutions
	times <- seq(tInit, tEnd, by = DT/1)
	#simulate the ODE
	out <- ode(Y, times, func = "derivs", parms = parms,
		jacfunc = "jac", dllname = "sirx2",
		initfunc = "initmod", nout = 1, outnames = "Sum")
	#convert output to a data frame
	out <- as.data.frame(out)
	out <- out[seq(1,dim(out)[1], by=1),]

	#name the columns of the dataframe
	names(out) <- c('t','S','I1','I2','I12','R1','R2','R1I2','R2I1','R1R2','Forcing')
	#define aggregate variables-
	#out$s1 <- out$S + out$I2 + out$R2
	#out$i1 <- out$I1 + out$I12 + out$R2I1
	#out$r1 <- out$R1+out$R1R2 + out$R1I2
	#out$s2 <- out$S + out$I1 + out$R1
	#out$i2 <- out$I2 + out$I12 + out$R1I2
	#out$r2 <- out$R2+out$R1R2 + out$R2I1
	#define ``contingency table'' variables
	out$q11 <- out$I12
	out$q00 <- (out$S + out$R1 + out$R2 + out$R1R2)
	out$q10 <- (out$I1 + out$R2I1)
	out$q01 <- (out$I2 + out$R1I2)
	#calculate the actual odds ratio
	#out$omega <- out$q11*out$q00/out$q10/out$q01
	#just return the last portion of the simulation to elminate transient dynamics
	#dat=out[out$t > max(times)-wndw,]
	dat=out
	return(dat)
}

getSimPred <- function(theta){
	r01  = exp(theta[1])
	r02  = exp(theta[2])
	gam1 = exp(theta[3])
	gam2 = exp(theta[4])
	rho1 = exp(theta[5])
	rho2 = exp(theta[6])	
	psi0 = plogis(theta[7])
	phi1 = plogis(theta[8])
	phi2 = plogis(theta[9])
	sigma= exp(theta[10])
	chi=exp(theta[11])
	eps=plogis(theta[12])
	phase1=52*plogis(theta[13])-26
	phase2=52*plogis(theta[14])-26
	phasepsi=52*plogis(theta[15])+18
	epspsi=plogis(theta[16])	
	
	sim = getSim(r01=r01,r02=r02,gam1=gam1,gam2=gam2,rho1=rho1,rho2=rho2,
		sigma=sigma,chi=chi,eps=eps,phase1=phase1,phase2=phase2) #assume altered susceptibility to coinfection
	return(sim)
}

plotSimPred <- function(theta,path1nm='Pathogen 1', path2nm='Pathogen 2', nm=NULL,save=F){
	sim=getSimPred(theta)
	sim=as.data.frame(sim)

	if(save==T){jpeg(file=nm,height=4,width=4,units='in',res=300)}
	plot((I1+I12+R2I1)~t,sim,type='l',col='blue',ylim=c(0,3e7))
	lines((I2+I12+R1I2)~t,sim,col='red')
	legend('topright',legend=c(path1nm,path2nm),lty=1,col=c('red','blue'))
	if(save==T){dev.off()}
}


getSimDat <- function(theta,DT=1/7,wndw=52*6,tEnd=8000, Y=runif(9)){
	r01  = exp(theta[1])
	r02  = exp(theta[2])
	gam1 = exp(theta[3])
	gam2 = exp(theta[4])
	rho1 = exp(theta[5])
	rho2 = exp(theta[6])	
	psi0 = plogis(theta[7])
	phi1 = plogis(theta[8])
	phi2 = plogis(theta[9])
	sigma= exp(theta[10])
	chi=exp(theta[11])
	eps=plogis(theta[12])
	phase1=52*plogis(theta[13])-26
	phase2=52*plogis(theta[14])-26
	phasepsi=52*plogis(theta[15])+18
	epspsi=plogis(theta[16])
	psi1=plogis(theta[17])


	Npop=327e6
	#simulate the ODEs	
	sim = getSim(r01=r01,r02=r02,gam1=gam1,gam2=gam2,rho1=rho1,rho2=rho2,
		sigma=sigma,chi=chi,eps=eps,phase1=phase1,phase2=phase2, #assume altered susceptibility to coinfection
		DT=DT, wndw=wndw, tEnd=tEnd,Y=Y)
	sim=tail(sim,n=1399)

	#aggregate in terms of the "contingency table" variables
	Q = sim[,c('q11','q00','q10','q01')]

	#testing probabilities for each state:
	#psi = psi0  * (1 + epspsi * sin(2*pi*(sim$t-phasepsi)/52)) 
	psi = (psi0+psi1*(sim$t-sim$t[1])/52 )  * (1 + epspsi * sin(2*pi*(sim$t-phasepsi)/52)) 


	f11 = 1-(1-psi)*(1-phi1)*(1-phi2)
	f00 = psi 
	f10 = 1-(1-psi)*(1-phi1)
	f01 = 1-(1-psi)*(1-phi2)
	
	#state-specific and time-dependent testing probabilities:
	f = data.frame( cbind(f11,f00,f10,f01) )
	Q = as.matrix(Q)

	#average number of tests each day:
	lambda = rowSums( f*Q )

	#sample the number of tests at each time:
	n = sapply(lambda,function(x) rpois(1,x))
	#n = sapply(lambda, function(x) rnbinom(1,mu=x,size=100) ) #no. of tests has more variance than poisson
	
	#allocate tests among different types of test results:
	x = t(apply(cbind(n,f*Q),1,function(y) rmultinom(1,unlist(y[1]),unlist(y[2:5])) ))
	x = as.data.frame(x)
	x[,5] <- n
	names(x) <- c('x11','x00','x10','x01','n')

	return(x)	
}

#define the log-likelihood
getlogl <- function(x,theta){

	Npop = 327e6
	r01  = exp(theta[1])
	r02  = exp(theta[2])
	gam1 = exp(theta[3])
	gam2 = exp(theta[4])
	rho1 = exp(theta[5])
	rho2 = exp(theta[6])	
	psi0 = plogis(theta[7])
	phi1 = plogis(theta[8])
	phi2 = plogis(theta[9])
	sigma= exp(theta[10])
	chi=exp(theta[11])
	eps=plogis(theta[12])
	phase1=52*plogis(theta[13])-26
	phase2=52*plogis(theta[14])-26
	phasepsi=52*plogis(theta[15])+18
	epspsi=plogis(theta[16])	
	
	#simulate the ODEs	
	sim = getSim(chi=chi,sigma=sigma,r01=r01,r02=r02,gam1=gam1,gam2=gam2,rho1=rho1,rho2=rho2,eps=eps,Npop=Npop,phase1=phase1,phase2=phase2) 
	sim = tail(sim, n=dim(x)[1])
	
	#aggregate in terms of the "contingency table" variables
	Q = sim[,c('q11','q00','q10','q01')]

	#testing probabilities for each state:
	psi = psi0  * (1 + epspsi * sin(2*pi*(sim$t-phasepsi)/52)) 

	f11 = 1-(1-psi)*(1-phi1)*(1-phi2)
	f00 = psi 
	f10 = 1-(1-psi)*(1-phi1)
	f01 = 1-(1-psi)*(1-phi2)
	
	#state-specific and time-dependent testing probabilities:
	f = data.frame( cbind(f11,f00,f10,f01) )
	Q = as.matrix(Q)

	logl <- sum( 
			apply( 
				cbind(x,f*Q), 1, function(x) sum(dpois(x[1:4],lambda=x[5:8],log=T)) 
			)
		)

	return(logl)	
}


#MCMC scheme
getEnsemble <- function(likelihood,thetaCurrent,dimTime,dimEnsemble){
        dimTheta = length(thetaCurrent)

        ensemble <- array( dim=c(dimTheta+1,dimTime,dimEnsemble) )

        #initialize the ensemble - Gaussian around the initial point
        for(i in 1:dimEnsemble){
                likvl <- -Inf

                while(likvl == -Inf ){
                prmvl <- rnorm(dimTheta,thetaCurrent,1e-2)
                likvl <- likelihood(prmvl)
        }

        ensemble[1:dimTheta,1,i] <- prmvl
        ensemble[dimTheta+1,1,i] <- likvl
        }
        #need to not have likelihoods = -Inf to start


        #generate all the independent stretches initially
        y <- runif(2*dimTime*dimEnsemble)
        a <- 10 #tuning parameter; stretches at most by a factor of a or compresses down by a factor of 1/a
        Z <- (y*(a-1)+1)^2/a #normalized the g(z) from Goodman and Weare (2010), generated unif(0,1)'s and inverted to get stretches
        Zind <- 1

        for(i in 1:(dimTime-1)){
         print(i)
         for(k in 1:dimEnsemble){
          Xk = ensemble[1:dimTheta,i,k]
          likelihood_Current = ensemble[1+dimTheta,i,k]
          j = sample(setdiff(1:dimEnsemble,k),1)
          Xj = ensemble[1:dimTheta,i,j]
          Y <<- Xj + (Xk-Xj)*Z[Zind]
	  likelihood_Proposed = likelihood(Y)
          if(is.na(likelihood_Proposed)){likelihood_Proposed=-Inf}
          proposal_Probability = Z[Zind]*exp( likelihood_Proposed - likelihood_Current ) #metropolis-hastings acceptance probability from Goodman and Weare (2010)
          if(runif(1)<proposal_Probability){
           ensemble[,i+1,k] <- c(Y,likelihood_Proposed)
          }else{
           ensemble[,i+1,k] <- c(Xk,likelihood_Current)
          }
          Zind <- Zind + 1
         }
        }
        return(ensemble)
}

#simdat <- getSimDat()[,-5]

#loglik <- function(theta) getlogl(simdat,theta)

#ens <- getEnsemble(loglik, theta_true , 2000, 25);save(ens,file='ModEstOutput16.Rdata')

#define the log-likelihood for just a single pathogen
getlogl1 <- function(x,theta,Y=runif(9),tEnd=8000){

	Npop = 327e6
	r01  = exp(theta[1])
	gam1 = exp(theta[2])
	rho1 = exp(theta[3])
	psi0 = plogis(theta[4])
	phi1 = plogis(theta[5])
	eps=plogis(theta[6])
	phase1=52*plogis(theta[7])-26
	phasepsi=52*plogis(theta[8])+18
	epspsi=plogis(theta[9])	
	psi1=plogis(theta[10])

	#normalize initial data
	Y = Y/sum(Y)

	#simulate the ODEs - arbitrary parameter values for the second pathogen which is independent	
	sim = getSim(r01=r01,r02=2,gam1=gam1,gam2=1,rho1=rho1,rho2=1/52,eps=eps,Npop=Npop,phase1=phase1,Y=Y, chi=1,sigma=1,nu=1,eta=1,tEnd=tEnd) 
	sim = tail(sim, n=dim(x)[1])
	
	#aggregate in terms of the "contingency table" variables
	Q = cbind(c(sim$q11+sim$q10),c(sim$q00+sim$q01))
	
	#independent testing probability (symptom-driven, but unrelated to the focal pathogen of interest); seasonal with a linear increase:
	psi = (psi0+psi1*(sim$t-sim$t[1])/52 )  * (1 + epspsi * sin(2*pi*(sim$t-phasepsi)/52)) 

	# testing probabilities for each state:
	f1 = 1-(1-psi)*(1-phi1)
	f0 = psi 
	
	#state-specific and time-dependent testing probabilities:
	f = data.frame( cbind(f1,f0)) 
	Q = as.matrix(Q)

	logl <- sum( 
			apply( 
				cbind(x,f*Q), 1, function(x) sum(dpois(x[1:2],lambda=x[3:4],log=T)) 
			)
		)

	return(logl)	
}


getSimDat1 <- function(theta,DT=1/7,wndw=52*6,tEnd=8000, Y=runif(9)){
	r01  = exp(theta[1])
	gam1 = exp(theta[2])
	rho1 = exp(theta[3])
	psi0 = plogis(theta[4])
	phi1 = plogis(theta[5])
	eps=plogis(theta[6])
	phase1=52*plogis(theta[7])-26
	phasepsi=52*plogis(theta[8])+18
	epspsi=plogis(theta[9])	
	psi1=plogis(theta[10])

	Npop=327e6
	# normalize the initial condition
	Y = Y/sum(Y)
	
	#simulate the ODEs - focal pathogen is the 1st, runs independent of the 2nd	
	sim = getSim(r01=r01,gam1=gam1,rho1=rho1,
		eps=eps,phase1=phase1,DT=DT,wndw=wndw,tEnd=tEnd, Y=Y,chi=1,sigma=1,eta=1,nu=1)
	sim=tail(sim,n=1399)
	#aggregate in terms of the "contingency table" variables
	Q = cbind(c(sim$q11+sim$q10),c(sim$q00+sim$q01))

	#testing probabilities for each state:
	psi = (psi0+psi1*(sim$t-sim$t[1])/52 )  * (1 + epspsi * sin(2*pi*(sim$t-phasepsi)/52)) 

	f1 = 1-(1-psi)*(1-phi1)
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
	x[,3] <- n
	names(x) <- c('x1','x0','n')
	
	#print(Y)
	return(x)	
}

getSimMean1 <- function(theta,DT=1/7,wndw=52*6,tEnd=8000, Y=runif(9),nrows=1399,Npop=327e6){
	r01  = exp(theta[1])
	gam1 = exp(theta[2])
	rho1 = exp(theta[3])
	psi0 = plogis(theta[4])
	phi1 = plogis(theta[5])
	eps=plogis(theta[6])
	phase1=52*plogis(theta[7])-26
	phasepsi=52*plogis(theta[8])+18
	epspsi=plogis(theta[9])	
	psi1=plogis(theta[10])

	#Npop=327e6
	#simulate the ODEs	
	sim = getSim(r01=r01,gam1=gam1,rho1=rho1,
		eps=eps,phase1=phase1,DT=DT,wndw=wndw,tEnd=tEnd, Y=Y,chi=1,sigma=1,eta=1,nu=1,Npop=Npop)
	sim=tail(sim,n=nrows)
	#aggregate in terms of the "contingency table" variables
	Q = cbind(c(sim$q11+sim$q10),c(sim$q00+sim$q01))

	#testing probabilities for each state:
	psi = (psi0+psi1*(sim$t-sim$t[1])/52 )  * (1 + epspsi * sin(2*pi*(sim$t-phasepsi)/52)) 

	f1 = 1-(1-psi)*(1-phi1)
	f0 = psi 
	
	#state-specific and time-dependent testing probabilities:
	f = data.frame( cbind(f1,f0) )
	Q = as.matrix(Q)

	#average number of tests each day:
	lambda = rowSums( f*Q )

	x <- f*Q
	x <- cbind(x, lambda)
	names(x) <- c('x1','x0','n')

	return(x)	
}

getSimPrev1 <- function(theta,DT=1/7,wndw=52*6,tInit=0,tEnd=8000, Y=runif(9),nrows=1399){
	r01  = exp(theta[1])
	gam1 = exp(theta[2])
	rho1 = exp(theta[3])
	psi0 = plogis(theta[4])
	phi1 = plogis(theta[5])
	eps=plogis(theta[6])
	phase1=52*plogis(theta[7])-26
	phasepsi=52*plogis(theta[8])+18
	epspsi=plogis(theta[9])	
	psi1=plogis(theta[10])

	Npop=327e6
	# normalize the initial condition
	Y = Y/sum(Y)
	
	#simulate the ODEs - focal pathogen is the 1st, runs independent of the 2nd	
	sim = getSim(r01=r01,gam1=gam1,rho1=rho1,
		eps=eps,phase1=phase1,DT=DT,wndw=wndw,tInit=tInit,tEnd=tEnd, Y=Y,chi=1,sigma=1,eta=1,nu=1)
	sim=tail(sim,n=nrows)
	sim$i1 <- sim$I1 + sim$I12 + sim$R2I1
	sim$s1 <- sim$S + sim$I2 + sim$R2
	sim$r1 <- sim$R1 + sim$R1R2 + sim$R1I2	

	return(sim)	
}

