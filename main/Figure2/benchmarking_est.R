#  In the paper, we set numiters = 2000, nparticles = 30 to generate 
#   Figure 3 of the main text.
#  Initializing from the true values, numiters = 100 may be sufficient
#  This defualt is set to a small value to ensure the code runs properly:

numiters = 5


# this function fits the model using the ensemble MCMC method
#  to the simulation in the simdats folder (with the name simdatnm)
#  by default this uses 25 particles for 5 timesteps - this probably

getEst <- function(simdatnm,niters = numiters,nparticles=30){
        load(paste0('Figure3/simdats/',simdatnm))
        loglik <- function(theta) getlogl(simdat,theta)
        theta_true <- simpars[simdatnm,]
        ens <- getEnsemble(loglik, theta_true , niters, nparticles)
        save(ens,file=paste0('Figure3/simdatfits/',simdatnm))
}

# run this line to fit the models to the simulations with the 
#  ensemble MCMC method:
mclapply(simdatlist, getEst, mc.cores=detectCores())

getbnchres <- function(simnm){
        # move into the simdatfits working directory
        load(paste0('Figure3/simdatfits/',simnm))
        simind = which(simdatlist == simnm)
        #all of these have loglik in ens[11,,]
        loglikmax = max(ens[11,,])
        #filter out values less than 15 loglikelihoods from loglikmax
        ensinds = which(ens[11,,]>loglikmax-15)
        maxind = which(ens[11,,]==loglikmax)[1] #there might be multiple
        mles = sapply(1:10, function(x) ens[x,,][maxind])
        lower_95 = sapply(1:10, function(x) quantile(ens[x,,][ensinds],.025) )
        upper_95 = sapply(1:10, function(x) quantile(ens[x,,][ensinds],.975) )
        true_values = simpars[simnm,]
        sigmas = sapply(1:10, function(x) sd(ens[x,,][ensinds])   )
        zval = (mles-true_values)/sigmas

        res <- rbind(lower_95,true_values,upper_95,mles,zval )
        colnames(res) <- colnames(simpars)
        return(res)
}



simnms = list.files(path='Figure3/simdatfits')
## some of these may have run into errors and stopped. check the dimensions and toss out
## any runs which ran into errors
## change the 2000 to however many timesteps you use in the fitting procedure:
#goodsims <- sapply(simnms, function(x){load(paste0('Figure3/simdatfits/',x));dim(ens)[2]==2000})
#goodsims <- simnms[goodsims]


getbnchreslist <- function(){
        bnchreslist <- lapply(simdatlist, getbnchres)
        return(bnchreslist)
}


res <- getbnchreslist()

## mmake some plots of estimates vs true values
r0 = sapply(res, function(x) exp(x['true_values','r0']))
r0hat = sapply(res, function(x) exp(x['mles','r0']))

gam = sapply(res, function(x) exp(x['true_values','gam']))
gamhat = sapply(res, function(x) exp(x['mles','gam']))

rho = sapply(res, function(x) exp(x['true_values','rho']))
rhohat = sapply(res, function(x) exp(x['mles','rho']))

psi0 = sapply(res, function(x) plogis(x['true_values','psi0']))
psi0hat = sapply(res, function(x) plogis(x['mles','psi0']))

phi = sapply(res, function(x) plogis(x['true_values','phi']))
phihat = sapply(res, function(x) plogis(x['mles','phi']))

eps = sapply(res, function(x) plogis(x['true_values','eps']))
epshat = sapply(res, function(x) plogis(x['mles','eps']))

phase = sapply(res, function(x) 52*plogis(x['true_values','phase'])-26 )
phasehat = sapply(res, function(x) 52* plogis(x['mles','phase'])-26)

phasepsi = sapply(res, function(x) 52*plogis(x['true_values','phasepsi'])+18 )
phasepsihat = sapply(res, function(x) 52* plogis(x['mles','phasepsi'])+18)

epspsi = sapply(res, function(x) plogis(x['true_values','epspsi']) )
epspsihat = sapply(res, function(x)  plogis(x['mles','epspsi']))


# reproduce Figure 3 of the main text (boxplots and histograms of the benchmarking results):
#pdf(file='Figure3.pdf')
par(mfrow=c(3,3),mar=c(5,5,4,2))
boxplot(r0hat~r0,xlab=bquote('Simulated '~R[0]),
	ylab= bquote('Estimated '~R[0]),
	main='Basic reproduction number',cex.lab=2,cex.axis=1.3)
abline(h=unique(r0),lty='dotted')

hist(c(1/gamhat),xlab=bquote('Estimated '~gamma^{-1}),freq=F,
	main='Mean duration of infection',cex.lab=2,cex.axis=1.3)

boxplot(c(1/rhohat)~c(1/rho),xlab=bquote('Simulated '~rho^{-1}),
	ylab=bquote('Estimated '~rho^{-1}),
	main='Mean duration of immunity',cex.lab=2,cex.axis=1.3)
abline(h=unique(1/rho),lty='dotted')


# just did a single estimation for a bunch of others... probably we need to run more. Damn.
hist(phasehat,xlab=bquote('Estimated '~tau[beta]),freq=F,main='Phase of maximum transmission',
	cex.lab=2,cex.axis=1.3)
abline(v=unique(phase))
hist(epshat,xlab=bquote('Estimated '~epsilon[beta]),freq=F,main='Amplitude of seasonality',
	cex.lab=2,cex.axis=1.3)
abline(v=unique(eps))
hist(phihat,xlab=bquote('Estimated '~phi),main='Virus-attributable testing rate',freq=F,
	cex.lab=2,cex.axis=1.3)
abline(v=unique(phi))


hist(psi0hat,xlab=bquote('Estimated '~psi[0]),main='Average baseline testing rate',freq=F,
	cex.lab=2,cex.axis=1.3)
hist(epspsihat,xlab=bquote('Estimated '~epsilon[psi]),main='Amplitude in baseline testing rate',freq=F,
	cex.lab=2,cex.axis=1.3)
hist(phasepsihat,xlab=bquote('Estimated '~tau[psi]),main='Phase of max baseline testing rate',freq=F,
	cex.lab=2,cex.axis=1.3)

#dev.off()
