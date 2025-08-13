plotFitSim <- function(dat,theta,DT=1/7,wndw=52*6,tEnd=8000,Y=runif(3),save=F,filnm=F,pathnm='Pathogen',start=16801,end=18199){
	tms= seq(start,end,by=1) #list the dates on the x-axis
	simdat=getSimDat(theta,DT=DT,wndw=wndw,tEnd=tEnd,Y=Y)	
	if(save==T){pdf(file=filnm,height=8,width=8)}
	par(mfrow=c(2,1))
	plot(as.Date(tms),dat[,1],ylim=c(0,1.3*max(dat[,1])),col='gray',main=paste0(pathnm, ' Positive'),ylab='No. of tests',
		xlab='Time [years]',type='l')
	lines(as.Date(tms),simdat$y,col='red',lwd=.2)
	legend('topright',lty=1,col=c('gray','red'),legend=c('Data','Model'))	
	plot(as.Date(tms),dat[,2],col='gray',main=paste0(pathnm,' Negative'),ylab='No. of tests',xlab='Time [years]',type='l')	
	lines(as.Date(tms),simdat$x,col='red',lwd=.2)
	#mtext(side=3,outer=TRUE,line=-1,text=pathnm)	
	if(save==T){dev.off()}
}

plotFitMean <- function(dat,theta,
	DT=1/7,
	wndw=52*6,
	tEnd=8000,
	Y=runif(3),
	save=F,
	filnm=F,
	pathnm='Pathogen',
	start=16801,
	end=18199,
	Npop=327e6,
	tExtra=0,
	phase=0,
	add=F,
	fittm=18199){
	# this is code is idiotically written. tEnd is measured in weeks, but start and end are julian days.
	# the point is that we simulated the ode for 8000 weeks to make sure we're at an attractor,
	# and the fitting procedure involved choosing parameters so that the last 1399 days (199.8571 weeks) 
	# of that 8000-week long simulation match the data. So, with the default tEnd=8000, day 1 of the data 
	# fitting period coincides with 7800.143 weeks' worth of simulation from the randomly chosen initial 
	# condition.
	# tExtra is a variable that allows us to generate predictions for tExtra days after the end of the data 
	# fitting period, which ends on julian=18199, or Oct 31, 2019 by default.	
	# phase=0 is the default; for some of the out-of-phase solutions (229E and PIV2), set phase=52

	end = end + tExtra
	tEnd = tEnd + tExtra/7 + phase
	tmslbl= seq(start,end,by=1) #list the dates on the x-axis
	simdat=getSimMean(theta,DT=DT,wndw=wndw,tEnd=tEnd,Y=Y,nrows=dim(dat)[1],Npop=Npop)	
	if(save==T){pdf(file=filnm,height=8,width=8)}
	if(add==F){par(mfrow=c(2,1))}
	plot(as.Date(tmslbl),dat[,1],ylim=c(0,1.3*max(dat[,1])),col='gray',main=paste0(pathnm,' Positive'),ylab='No. of tests',
		xlab='Time [years]',type='p')
	abline(v=as.Date(fittm))
	abline(v=as.Date(18335))
	lines(as.Date(tmslbl),simdat$y,col='black',lwd=1)	
	legend('topleft',lty=1,col=c('gray','black'),legend=c('Data','Model'))	
	plot(as.Date(tmslbl), dat[,2],col='gray',main=paste0(pathnm,' Negative'),ylab='No. of tests',xlab='Time [years]',type='p')	
	abline(v=as.Date(fittm))
	abline(v=as.Date(18335))
	lines(as.Date(tmslbl),simdat$x,col='black',lwd=1)
	if(save==T){dev.off()}

	# calculate RSS from the model, and from supsmu:
	rss_mod = list(pos=sum((dat[,1]-simdat$y)^2), neg = sum((dat[,2]-simdat$x)^2) )
	# fit supsmu:
	s_pos = supsmu(as.Date(tmslbl),dat[,1])
	s_neg = supsmu(as.Date(tmslbl),dat[,2])
	rss_supsmu = list(pos = sum((s_pos[[2]]-dat[,1])^2), neg = sum((s_neg[[2]]-dat[,2])^2) )
	rss_mod = unlist(rss_mod)
	rss_supsmu = unlist(rss_supsmu)	
	x = c(rss_mod[1],rss_supsmu[1],rss_mod[1]/rss_supsmu[1], 
		rss_mod[2],rss_supsmu[2],rss_mod[2]/rss_supsmu[2])
	x = round(x,2)
	return(x)
}

plotFitMeanPos <- function(dat,theta,
	DT=1/7,
	wndw=52*6,
	tEnd=8000,
	Y=runif(3),
	save=F,
	filnm=F,
	pathnm='Pathogen',
	start=16801,
	end=18199,
	Npop=327e6,
	tExtra=0,
	phase=0,
	pathcol='black',
	lin1=18199){
	# this is code is idiotically written. tEnd is measured in weeks, but start and end are julian days.
	# the point is that we simulated the ode for 8000 weeks to make sure we're at an attractor,
	# and the fitting procedure involved choosing parameters so that the last 1399 days (199.8571 weeks) 
	# of that 8000-week long simulation match the data. So, with the default tEnd=8000, day 1 of the data 
	# fitting period coincides with 7800.143 weeks' worth of simulation from the randomly chosen initial 
	# condition.
	# tExtra is a variable that allows us to generate predictions for tExtra days after the end of the data 
	# fitting period, which ends on julian=18199, or Oct 31, 2019.	
	# phase=0 is the default; for some of the out-of-phase solutions (229E and PIV2), set phase=52

	end = end + tExtra
	tEnd = tEnd + tExtra/7 + phase
	tmslbl= seq(start,end,by=1) #list the dates on the x-axis
	simdat=getSimMean(theta,DT=DT,wndw=wndw,tEnd=tEnd,Y=Y,nrows=dim(dat)[1],Npop=Npop)	
	if(save==T){pdf(file=filnm,height=8,width=8)}
	par(mfrow=c(1,1))
	plot(as.Date(tmslbl),dat[,1],ylim=c(0,1.3*max(dat[,1])),col=pathcol,main=paste0(pathnm),ylab='No. of tests',
		xlab='Time [years]',type='p', cex.lab=2.8, cex.axis=2.8, cex.main=3)
	abline(v=as.Date(lin1))
	abline(v=as.Date(18335))
	lines(as.Date(tmslbl),simdat$y,lwd=1,col='black')	
	#legend('topleft',lty=1,col=c('gray','black'),legend=c('Data','Model'))	
	

	if(save==T){dev.off()}
}


plotFitResiduals <- function(dat,theta,DT=1/7,wndw=52*6,tEnd=8000,Y=runif(3),save=F,filnm=F,pathnm='Pathogen'){
	simdat=getSimMean(theta)	
	if(save==T){pdf(file=filnm,height=8,width=8)}
	par(mfrow=c(2,2))
	plot(dat[,1]-simdat$y,
	     #ylim=c(0,1.3*max(dat[,1])),
	     col='gray',main='Positive',ylab='No. of tests - Predicted No. of tests',
		xlab='Time [days]',type='l')
	hist(dat[,1]-simdat$y,main='Positive',xlab='No. of tests - Predicted No.')
	plot(dat[,2]-simdat$x,col='gray',main='Negative',
	     ylab='No. of tests - Predicted No. of tests',xlab='Time [days]',type='l')	
	hist(dat[,2]-simdat$x,main='Negative',xlab='No of tests - Predicted No.')
	if(save==T){dev.off()}
}

getMLE <- function(data){
	loglik= function(theta){
		getlogl1(getinfdat(data),theta)
	}
	opt=optim(theta_init,loglik,control=list(fnscale=-1,maxit=5000))
	return(opt)

}

# modify getlogl1 to handle errors produced by deSolve; if it throws an error, 
# return -Inf for the log-likelihood

trylogl <- function(x,theta,Y=runif(3), tEnd=8000){
	logl1 <- try({
		getlogl1(x,theta,Y=Y,tEnd=tEnd)
		},silent = TRUE)
	if(class(logl1)=='try-error'){logl1<- -Inf}
	return(logl1)
}

# for models with multistability, use the max likelihood across several randomly selected
# initial conditions

# we can probably use the same initial conditions every parameter update
initList <- list()
for(i in 1:20) initList[[i]] <- runif(3)
rm(i)

initList <- list()
initList[[1]] <- c(.9,.01,.09)
initList[[2]] <- c(.9,.09,.01)
initList[[3]] <- c(.9,.1,0)




multilogl <- function(x,theta, initlist, tEnd=8000){
	fn <- function(init) trylogl1(x,theta,Y=init,tEnd=tEnd)	
	logliks <- mclapply(initlist, fn, mc.cores=length(initList))
	loglik <- max(unlist(logliks))
	return(loglik)
}

# mle = mclapply(goodpnms, getMLE, mc.cores=62);names(mle)<-pairnms;save(mle,file='sept_19_opts.Rdata')

# extract the MLE from the mcmc ensemble

getmcmcmle <- function(ens){
	loglikind = dim(ens)[1]
	lmax = max(ens[loglikind,,])
	ind1 = which(apply(ens[loglikind,,], 1, function(x) sum(which(x==lmax)) )>0)[1]
	ind2 = which(ens[loglikind,ind1,] == lmax)[1]
	mcmcmle = ens[,ind1,ind2]
	return(mcmcmle)
}


# write some functions to plot mcmc results

# this accepts the 11-dimensional ensemble output for the single-pathogen models
plotEpiPar <- function(ens, pathnm='Simulated Pathogen', save=F, filnm=F,wd='',t0=7800.286){
	# the 2nd dimension of the mcmcensemble array, ens, corresponds to timesteps in the MCMC
	tmind = dim(ens)[2]

	# In a real analysis we would want to discard transients from the MCMC. Do that here:	
	#ens = ens[,seq(tmind-500,tmind,by=1),]

	# transform parameters back from R^N, and identify reasonable plotting regions
	# in the vicinity of the MLE:
	lmax=max(ens[dim(ens)[1],,])
	r01 = exp(c(ens[1,,]))
	eps = plogis(c(ens[6,,]))
	phase1 = 52*plogis(c(ens[7,,]))-26
	phase1 = ((phase1-t0) %% (52))*7 + 16800 #rescale phase1 to give max transmission rate rel. to. Dec31/Jan1 2016 
	gam1 = exp(c(ens[2,,]))
	rho1 = exp(c(ens[3,,]))
	phi1 = plogis(c(ens[5,,]))
	loglik = c(ens[11,,])

	# plot the profile likelihoods, and save as a .jpeg if desired:
	if(save==T){jpeg(file=filnm,height=8,width=13,res=400, units='in')}
	par(mfrow=c(2,3))
	plot(r01,loglik,xlab=expression(R[0]), ylab='',#ylim=c(lmax-20,lmax+5),
		cex.lab=2.4, cex.axis=2, cex.main=2)
	title(ylab = 'log-likelihood', cex.lab=2)	
	plot(eps,loglik,xlab=expression(epsilon), ylab='', #ylim=c(lmax-20,lmax+5),
		cex.lab=2.4, cex.axis=2, cex.main=2)
	title(ylab = 'log-likelihood', cex.lab=2)	
	plot(as.Date(phase1),loglik,
	     xlab=expression(tau[beta]), 
	     ylab='',#ylim=c(lmax-20,lmax+5),
	     cex.lab=2.4, cex.axis=2, cex.main=2)
	title(ylab = 'log-likelihood', cex.lab=2)	
	plot(gam1^(-1),loglik,xlab=expression(gamma^{-1}), 
	     ylab='',#ylim=c(lmax-20,lmax+5),
	     cex.lab=2.4, cex.axis=2, cex.main=2)
	title(ylab = 'log-likelihood', cex.lab=2)	
	plot(rho1^(-1),loglik,xlab=expression(rho^{-1}), 
	     ylab='',#ylim=c(lmax-20,lmax+5), 
	     cex.lab=2.4, cex.axis=2, cex.main=2)
	title(ylab = 'log-likelihood', cex.lab=2)	
	plot(phi1,loglik,xlab=expression(phi), ylab='',#ylim=c(lmax-20,lmax+5),
		cex.lab=2.4, cex.axis=2, cex.main=2)
	title(ylab = 'log-likelihood', cex.lab=2)	
	mtext(side=3,outer=TRUE,line=-3,text=pathnm, cex=3)	
	if(save==T){dev.off()}
}


plotPsiPar <- function(ens, pathnm='Simulated Pathogen', save=F, filnm=F,wd='oct12/mcmcoutput/',t0=7800.286){
	# the 2nd dimension of the mcmc array ens corresponds to timesteps in the MCMC	
	tmind = dim(ens)[2]

	# in a real analysis we would want to discard MCMC transients. Do that here:
	#ens = ens[,seq(tmind-500,tmind,by=1),]

	# Transform parameters back, and identify plotting regions:
	lmax=max(ens[dim(ens)[1],,])
	psi0 = plogis(c(ens[4,,]))
	psi1 = plogis(c(ens[10,,]))
	epspsi = plogis(c(ens[9,,]))
	phasepsi = 52*plogis(ens[8,,])+18
	phasepsi = ((phasepsi-t0+13) %% (52))*7 + 16800 #rescale phasepsi to give max transmission rate rel. to. Dec31/Jan1 2016 
	loglik = c(ens[11,,])

	# plot the profile likelihoods and save the file as a jpeg, if desired
	if(save==T){jpeg(file=filnm,height=8,width=8,res=400,units='in')}
	par(mfrow=c(2,2))
	plot(psi0,loglik,xlab=expression(psi[0]), ylab='',
		cex.lab=2, cex.axis=1.5, cex.main=2)#ylim=c(lmax-20,lmax+5))
	title(ylab = 'log-likelihood', cex.lab=2)	
	plot(psi1,loglik,xlab=expression(psi[1]~' [per year]'), ylab='',
		cex.lab=2, cex.axis=1.5, cex.main=2)#ylim=c(lmax-20,lmax+5))
	title(ylab = 'log-likelihood', cex.lab=2)	
	plot(epspsi,loglik,xlab=expression(epsilon[psi]), ylab='',
		cex.lab=2, cex.axis=1.5, cex.main=2)#ylim=c(lmax-20,lmax+5))
	title(ylab = 'log-likelihood', cex.lab=2)	
	plot(as.Date(phasepsi),loglik,xlab=expression(tau[psi]), 
	     ylab='',cex.lab=2, cex.axis=1.5, cex.main=2) #,ylim=c(lmax-20,lmax+5))
	title(ylab = 'log-likelihood', cex.lab=2)	
	mtext(side=3,outer=TRUE,line=-3,text=pathnm, cex=3)	
	if(save==T){dev.off()}
}


plotSimPrev <- function(theta,pathnm,start=16801,end=18199,filnm=F,save=F){
	sim=getSimPrev(theta)
	sim$t = sim$t - sim$t[1]
	tms = seq(start,end,by=1)
	if(save==T){pdf(file=filnm,height=4,width=4)}
	par(mfrow=c(1,1))
	plot(log10(S)~as.Date(tms),sim,type='l',lty='dashed',ylim=c(4,9),
		main=paste0(pathnm,' Estimated Prevalence'),
		ylab=expression(log[10]~'Population'),
		xlab='Time [years]')
	lines(log10(R)~as.Date(tms),sim,lty='dotted')
	lines(log10(I)~as.Date(tms),sim)
	legend('bottomleft',legend=c('S','I','R'),lty=c(2,1,3))
}

#plot a collection of trajectories from the posterior, with uncertainty in inital conditions if desired:
plotSimPrevEns <- function(ensemble,pathnm,start=16801,end=18199,tEnd=8000,filnm=F,save=F,nsims=10,Y=runif(3),a=0,b=1){
	# just select the end of the mcmc:	
	theta=getmcmcmle(ensemble)
	lmax=theta[length(theta)]	
	# find the parameters that fall between 11 and 9 log-likelihoods of the MLE:	
	# changed it to between 2 and 3, to mimic a 95% confidence set:
	ci99 = ensemble[dim(ensemble)[1],,] > lmax-3 & ensemble[dim(ensemble)[1],,] < lmax-1 	
	inds = which(apply(ci99,1, sum)	> 0)
	ci99 = ensemble[,inds,]
	tms = seq(start,end,by=1)
	sim=getSimPrev(theta,tEnd=tEnd,nrows=length(tms),Y=Y)
	#tms = seq(start,end,by=1)
	#sim=tail(sim,n=length(tms))
	if(save==T){pdf(file=filnm,height=6,width=6)}
	par(mfrow=c(1,1))
	plot(log10(S)~as.Date(tms),sim,type='l',lty='dashed',ylim=c(4,9),
		main=paste0(pathnm,' Estimated Prevalence'),
		ylab=expression(log[10]~'Population'),
		xlab='Time [years]',col='blue')
	lines(log10(R)~as.Date(tms),sim,lty='dotted',col='darkgreen')
	lines(log10(I)~as.Date(tms),sim,col='red')
	legend('bottomleft',legend=c('S','I','R'),lty=c(2,1,3),col=c('blue','black','darkgreen'))

	# save the initial condition to jitter as well:
	y=sim[1,c(2,3,4)]

	for(i in 2:nsims){
		theta = ci99[,sample(1:dim(ci99)[2],1),sample(1:dim(ci99)[3],1)]	
		# comment this line out to draw parameters from the 99% confidence region boundary:	
		if(b==0){theta=getmcmcmle(ensemble)}

		yinit = as.numeric(y)
		#jitter the initial condition - comment out to just visualize uncertainty in parameters
		yinit = jitter(yinit,amount=a*yinit[2])

		#just run the simulation from Jan 1, 2016 to the end of the data fitting period (encoded in tms), 
		# which counts the number of weeks (but solver operates on timescale of days)
		sim=getSimPrev(theta,tEnd=length(tms)/7,nrows=length(tms),Y=yinit)
		sim=tail(sim,n=length(tms))	
		lines(log10(S)~as.Date(tms),sim,lty='dashed',col='blue')				
		lines(log10(R)~as.Date(tms),sim,lty='dotted',col='darkgreen')
		lines(log10(I)~as.Date(tms),sim)	
	}
	if(save==T){dev.off()}
}


# calculate the date when the transmission rate is highest
getMaxBeta <- function(mle,nrows=1399,start=16801,end=18199){
	sim=getSim()
	sim=tail(sim,n=nrows)
	tms = seq(start,end,by=1)
	phase1 = 52*plogis(mle['logitphase1'])-26
	crud = (sim$t-phase1-13) %% 52
	inds = which(diff(crud) < 0)
	max = as.Date(mean(tms[inds+1]))
	return(max)
}

# calculate the characteristic frequency of the SIRS system:
getOmega <- function(theta){
	delta = exp(c(theta[3]))
	gamma = exp(c(theta[2]))
	R0 = exp(c(theta[1]))
	omega = Im(sqrt(as.complex((delta^3 - 2*gamma*(R0 - 2)*delta^2 + gamma^2*(R0^2 - 8*R0 + 8)*delta - 4*gamma^3*(R0 - 1))*delta/(delta + gamma)^2))/2)

	return(omega)
}

getOmegaDist <- function(ens){
	# use the last 500 mcmc iterates:
	niter = dim(ens)[2]

	# in a real analysis, discard transients from the MCMC:	
	#ens = ens[,seq(niter-500,niter,by=1),]
	logr01 = c(ens[1,,])
	loggam1 = c(ens[2,,])
	logrho1 = c(ens[3,,])
	thetas = cbind(logr01,loggam1,logrho1)
	colnames(thetas) <- c('logr01','loggam1','logrho1')
	omegas = apply(thetas, 1, getOmega)
	return(omegas)
}

plotPdDist <- function(ens,save=F,filnm=F){
	omegas <- getOmegaDist(ens)
	if(save==T){pdf(filnm,height=6,width=6)}
	par(mfrow=c(1,1))	
	hist(2*pi/omegas,xlab='Characteristic period of oscillation [weeks]',
		main=pathnm,freq=F)
	if(save==T){dev.off()}

}

getParDist <- function(par,pathnms=goodpnms){
	pars = c('r0','gam','rho','psi0','phi','eps','phase1','phasepsi','epspsi','psi1')	
	ind = which(pars==par)	
	pardist = c()	
	for(nm in pathnms){
		load(paste0(wd,nm,'.Rdata'))
		pardist = cbind(pardist,c(ens[ind,c((dim(ens)[2]-100):dim(ens)[2]),]))
	
	}
	colnames(pardist) <- pathnms
	return(pardist)
}

# like simPrevEns, but then jitter the terminal solution value to model a disruption:
plotDisruption <- function(ensemble,pathnm,start=16801,end=18199,tEnd=8000,filnm=F,save=F,Y=runif(3),a=0,b=0,thresh=1,xmax=1500,simmax=1000){
	# just select the end of the mcmc:	
	theta=getmcmcmle(ensemble)
	lmax=theta[length(theta)]	
	# find the parameters that fall between 11 and 9 log-likelihoods of the MLE:	
	# changed it to between 2 and 3, to mimic a 95% confidence set:
	ci99 = ensemble[dim(ensemble)[1],,] > lmax-3 & ensemble[dim(ensemble)[1],,] < lmax-1 	
	inds = which(apply(ci99,1, sum)	> 0)
	ci99 = ensemble[,inds,]
	tms = seq(start,end,by=1)
	sim=getSimPrev(theta,tEnd=tEnd,nrows=length(tms),Y=Y)
	#tms = seq(start,end,by=1)
	#sim=tail(sim,n=length(tms))
	if(save==T){pdf(file=filnm,height=6,width=6)}
	par(mfrow=c(1,1))
	plot(log10(S)~as.Date(tms),sim,type='l',lty='dashed',ylim=c(4,9),xlim=c(as.Date(tms[1]),as.Date(tail(tms,n=1)+xmax)),
		main=paste0(pathnm,' Estimated and Projected Prevalence'),
		ylab=expression(log[10]~'Population'),
		xlab='Time [years]',col='black')
	lines(log10(R)~as.Date(tms),sim,lty='dotted',col='black')
	lines(log10(I)~as.Date(tms),sim,col='black')
	legend('bottomleft',legend=c('S','I','R'),lty=c(2,1,3),col=c('black','black','black'))

	# save the initial condition to jitter as well:
	y=sim[dim(sim)[1],]
	y=y[2:4]

	inew = unlist(y['I']) * c(1,0.5,0.2,0.05)
	pnew = unlist(y['S']/(y['S']+y['R']))
	pnew = seq(pnew,pnew+(1-pnew)/2,length=3)
	newinits = data.frame(expand.grid(inew,pnew))
	colnames(newinits) <- c('I','p')
	newinits$S <- newinits$p * (327e6 - newinits$I)
	newinits$R <- 327e6 - newinits$I - newinits$S
	newinits <- newinits[,c('S','I','R')]

	switch <- rep(0,dim(newinits)[1])
	for(i in 1:dim(newinits)[1]){
		#theta=ensemble[,sample( seq(dim(ensemble)[2]-100, dim(ensemble)[2]), 1) ,sample(dim(ensemble)[3],1)]	
		#theta = ci99[,sample(1:dim(ci99)[2],1),sample(1:dim(ci99)[3],1)]	
		# comment this line out to draw parameters from the 99% confidence region boundary:	
		if(b==0){theta=getmcmcmle(ensemble)}

		#jitter the initial condition - comment out to just visualize uncertainty in parameters
		yinit <- newinits[i,]
		
		#just run the simulation from Jan 1, 2016 to the end of the data fitting period (encoded in tms), 
		# which counts the number of weeks (but solver operates on timescale of days)
		sim1=getSimPrev(theta,tInit = tail(sim$t,n=1),tEnd=tail(sim$t,n=1)+simmax,nrows=1000000,Y=as.numeric(yinit))
		tms2 = c(end:(end+dim(sim1)[1]-1))	
		#sim=tail(sim,n=length(tms))
		# color transients according to their attractor
	
		if(i==1){att=sim1$i1}		
		dist  = sum((tail(att,n=1000)-tail(sim1$i1,n=1000))^2)/327e6^2

		if(dist>thresh){cols='blue';switch[i]=1}else{cols = 'red'}
		if(i==1){cols='black'} 
		lines(log10(S)~as.Date(tms2),sim1,lty='dashed',col=cols)				
		lines(log10(R)~as.Date(tms2),sim1,lty='dotted',col=cols)
		lines(log10(I)~as.Date(tms2),sim1,col=cols)	
	}
	if(save==T){dev.off()}
	return(list(switch,newinits))
}
