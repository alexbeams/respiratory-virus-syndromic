# load this package to run computations in parallel:
require(parallel)

# load the package to solve ODEs:
require(deSolve)

#######
## Load in the other necessary files:
#######

# load in the model definitions and the loglikelihood function:
source('model_functions.R')

# load in the ensemble MCMC sampler:
source('ensemble_sampler.R')


# parameter combinations for benchmarking simulations
getSimSet <- function(){
        r0 <- log(c(2,4,8))
        gam <- log(1)
        rho <- log(c(1/52,1/104))
        psi0 <- -16
        phi <- -15
        eps <- qlogis(0.15)
	phase <- 0
	phasepsi <- 0
	epspsi <- qlogis(0.75)
	psi1 <- -Inf

        pars <- expand.grid(r0,gam,rho,psi0,phi,eps,phase,phasepsi,epspsi,psi1)
        return(pars)

}

# generate the parameter combinations:
simpars <- getSimSet()
simpars <- t(apply(simpars, 1, as.numeric))
colnames(simpars) <- c('r0','gam','rho','psi0','phi','eps','phase','phasepsi',
	'epspsi','psi1')

# generate nsims simulations for each parameter vector:
nsims = 30
simpars <- simpars[rep(seq_len(nrow(simpars)), each = nsims), ]

# generate the simulations and store them in the simdats folder
# they are indexed as simdatk.Rdata, where k is the row of simpars just above
getSimDatSet <- function(){
        for(i in 1:dim(simpars)[1]){
		simdat=getSimDat(unlist(simpars[i,]))
                save(simdat,file=paste0('Figure3/simdats/','simdat',i,'.Rdata'))
        }
}

Nsims =  dim(simpars)[1]

# make a convenient list of names of all of the simulations:
simdatlist <- paste0('simdat',1:Nsims,'.Rdata')

# assign these names to the rows of the dataframe of parameter values:
rownames(simpars) = simdatlist

# generate the simulations and store them in the simdats folder:
getSimDatSet()
