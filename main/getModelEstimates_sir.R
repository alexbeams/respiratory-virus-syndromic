rm(list=ls())

require(parallel) #to run parallel computations
require(pracma)   #to analyze ODE outputs for multistability, etc.
require(scales)   #for color transparency in plots

dat <- read.csv("CensusRegion_daily_RP_data_20211128.csv")  #n=194694 tests
#dat <- read.csv("CDCRegion_weekly_RP_data_20211128.csv")    #n=88725 tests
#dat <- read.csv("City_weekly_RP_data_20211128.csv")         #n=48537 tests; also seems to begin 2016-10-01... was this in error?

dat$julian <- as.numeric(as.Date(dat$StartTime)) #assign julian dates

# SARS.CoV.2 and B..parapertussis have a large number of NAs, but none of
# the other pathogens do, so we can just use the full data set

getinfdat <- function(nm1, lastday=18199,reg=c("Northeast","West","Midwest","South")){
	regioninds = which(dat$CensusRegionNational %in% reg)
        ddat=data.frame(cbind(dat[regioninds,c('julian',nm1,'RNA.Process.Control')]))
        ddat=ddat[ddat$julian <= lastday,]
        
	aggdat = aggregate(.~julian,ddat,sum)

	x = cbind(
		  aggdat[,nm1], 
		  aggdat$RNA.Process.Control - aggdat[,nm1] )
	
	return(cbind(x))
}

# for looking at totals of seasonal coronavirueses, influenza A's, and parainfluenzas, etc.
getinfdatgrp <- function(nms, lastday=18199){
        ddat=data.frame(cbind(dat[,c('julian',nms,'RNA.Process.Control')]))
       	pos = rowSums(ddat[,nms])
        pos[pos > 0] = 1
	ddat$pos = pos	
	ddat=ddat[ddat$julian <= lastday,]
       	ddat=ddat[,c('julian','pos','RNA.Process.Control')] 
	aggdat = aggregate(.~julian,ddat,sum)

	x = cbind(
		  aggdat[,'pos'], 
		  aggdat$RNA.Process.Control - aggdat[,'pos'] )
	
	return(cbind(x))
}

getddat <- function(nm1,nm2){
        ddat=data.frame(cbind(dat[,c('julian',nm1,nm2,'RNA.Process.Control')]))
        ddat=ddat[ddat$julian < 18200,]
        ddat$co=ddat[,nm1]+ddat[,nm2]
        ddat$co[ddat$co==1]=0
        ddat$co[ddat$co==2]=1
        ddat[,paste0(nm1,'0')]=ddat[,nm1]-ddat$co
        ddat[,paste0(nm2,'0')]=ddat[,nm2]-ddat$co

        aggdat = aggregate(.~julian,ddat,sum)

	x = cbind(
		  aggdat$co,
		  aggdat$RNA.Process.Control - aggdat$co-aggdat[,paste0(nm1,'0')]-aggdat[,paste0(nm2,'0')],
		  aggdat[,paste0(nm1,'0')],
		  aggdat[,paste0(nm2,'0')])
	
	return(cbind(x))
}

source('September/getModEst_sir.R')

pnms <- names(dat)[7:28]
goodpnms <- pnms[1:20] #B..parapertussis and SARS.CoV.2 have lots of NAs

plotFitSim1 <- function(dat,theta,DT=1/7,wndw=52*6,tEnd=8000,Y=runif(3),save=F,filnm=F,pathnm='Pathogen',start=16801,end=18199){
	tms= seq(start,end,by=1) #list the dates on the x-axis
	simdat=getSimDat1(theta,DT=DT,wndw=wndw,tEnd=tEnd,Y=Y)	
	if(save==T){pdf(file=filnm,height=8,width=8)}
	par(mfrow=c(2,1))
	plot(as.Date(tms),dat[,1],ylim=c(0,1.3*max(dat[,1])),col='gray',main=paste0(pathnm, ' Positive'),ylab='No. of tests',
		xlab='Time [years]',type='l')
	lines(as.Date(tms),simdat$x1,col='red',lwd=.2)
	legend('topright',lty=1,col=c('gray','red'),legend=c('Data','Model'))	
	plot(as.Date(tms),dat[,2],col='gray',main=paste0(pathnm,' Negative'),ylab='No. of tests',xlab='Time [years]',type='l')	
	lines(as.Date(tms),simdat$x0,col='red',lwd=.2)
	#mtext(side=3,outer=TRUE,line=-1,text=pathnm)	
	if(save==T){dev.off()}
}

plotFitMean1 <- function(dat,theta,
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
	simdat=getSimMean1(theta,DT=DT,wndw=wndw,tEnd=tEnd,Y=Y,nrows=dim(dat)[1],Npop=Npop)	
	if(save==T){pdf(file=filnm,height=8,width=8)}
	if(add==F){par(mfrow=c(2,1))}
	plot(as.Date(tmslbl),dat[,1],ylim=c(0,1.3*max(dat[,1])),col='gray',main=paste0(pathnm,' Positive'),ylab='No. of tests',
		xlab='Time [years]',type='p')
	abline(v=as.Date(fittm))
	abline(v=as.Date(18335))
	lines(as.Date(tmslbl),simdat$x1,col='black',lwd=1)	
	legend('topleft',lty=1,col=c('gray','black'),legend=c('Data','Model'))	
	plot(as.Date(tmslbl), dat[,2],col='gray',main=paste0(pathnm,' Negative'),ylab='No. of tests',xlab='Time [years]',type='p')	
	abline(v=as.Date(fittm))
	abline(v=as.Date(18335))
	lines(as.Date(tmslbl),simdat$x0,col='black',lwd=1)
	if(save==T){dev.off()}

	# calculate RSS from the model, and from supsmu:
	rss_mod = list(pos=sum((dat[,1]-simdat$x1)^2), neg = sum((dat[,2]-simdat$x0)^2) )
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
	simdat=getSimMean1(theta,DT=DT,wndw=wndw,tEnd=tEnd,Y=Y,nrows=dim(dat)[1],Npop=Npop)	
	if(save==T){pdf(file=filnm,height=8,width=8)}
	#par(mfrow=c(2,1))
	plot(as.Date(tmslbl),dat[,1],ylim=c(0,1.3*max(dat[,1])),col=pathcol,main=paste0(pathnm),ylab='No. of tests',
		xlab='Time [years]',type='p', cex.lab=2.8, cex.axis=2.8, cex.main=3)
	abline(v=as.Date(lin1))
	abline(v=as.Date(18335))
	lines(as.Date(tmslbl),simdat$x1,lwd=1,col='black')	
	#legend('topleft',lty=1,col=c('gray','black'),legend=c('Data','Model'))	
	

	if(save==T){dev.off()}
}


plotFitMean.constprevsis <- function(dat,theta) plotFitMean1(dat, c(theta[1:2],20,theta[3:4],-100,theta[5:8]))
plotFitSim.constprevsis <- function(dat,theta) plotFitSim1(dat, c(theta[1:2],20,theta[3:4],-100,theta[5:8]))


plotFitMeanReg1 <- function(pathnm,theta,DT=1/7,wndw=52*6,tEnd=8000,Y=runif(3),save=F,filnm=F,start=16801,end=18199,Npop=327e6){
	# load in region-specific data:
	dats = lapply(c('Northeast','South','Midwest','West'), function(x) getinfdat(pathnm, reg=x))
	regscls = c(0.17,0.38,0.21,0.24) #fraction of the population in Northeast, South, Midwest, and West
	
	tmslbl= seq(start,end,by=1) #list the dates on the x-axis
	if(save==T){pdf(file=filnm,height=8,width=8)}
	par(mfrow=c(2,2))	
	
	for(i in 1:4){	
		simdat=getSimMean1(theta,DT=DT,wndw=wndw,tEnd=tEnd,Y=Y,nrows=dim(dats[[i]])[1],Npop=Npop*regscls[i])

		plot(as.Date(tmslbl),dats[[i]][,1],ylim=c(0,1.3*max(dats[[i]][,1])),
		     col='black',main=c('Northeast','South','Midwest','West')[i],ylab='No. of tests',
		     xlab='Time [years]',type='l')
		lines(as.Date(tmslbl),simdat$x1,col='red',lwd=1)	
		#lines(as.Date(tmslbl), dats[[i]][,2],col='gray',lwd=1)	
		#lines(as.Date(tmslbl),simdat$x0,col='black',lwd=1)
	}	
	legend('topleft',lty=1,col=c('gray','black'),legend=c('Data','Model'))	
	mtext(pathnm,side=3,line=-2,outer=T)	
	if(save==T){dev.off()}
}

plotFitMean <- function(pathnm,wd='oct12/mcmcoutput/'){
	load(paste0(wd,pathnm,'.Rdata'))
	data=getinfdat(pathnm)
	thetahat=getmcmcmle(ens)
	plotFitMean1(data,thetahat)
}


plotFitResiduals1 <- function(dat,theta,DT=1/7,wndw=52*6,tEnd=8000,Y=runif(3),save=F,filnm=F,pathnm='Pathogen'){
	simdat=getSimMean1(theta)	
	if(save==T){pdf(file=filnm,height=8,width=8)}
	par(mfrow=c(2,2))
	plot(dat[,1]-simdat$x1,
	     #ylim=c(0,1.3*max(dat[,1])),
	     col='gray',main='Positive',ylab='No. of tests - Predicted No. of tests',
		xlab='Time [days]',type='l')
	hist(dat[,1]-simdat$x1,main='Positive',xlab='No. of tests - Predicted No.')
	plot(dat[,2]-simdat$x0,col='gray',main='Negative',
	     ylab='No. of tests - Predicted No. of tests',xlab='Time [days]',type='l')	
	hist(dat[,2]-simdat$x0,main='Negative',xlab='No of tests - Predicted No.')
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

trylogl1 <- function(x,theta,Y=runif(3), tEnd=8000){
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




multilogl1 <- function(x,theta, initlist, tEnd=8000){
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
plotEpiPar <- function(pathnm, save=F, filnm=F,wd='oct12/mcmcoutput/',t0=7800.286){
	load(paste0(wd,pathnm,'.Rdata')) #load in the mcmc output
	tmind = dim(ens)[2]
	ens = ens[,seq(tmind-500,tmind,by=1),]
	lmax=max(ens[dim(ens)[1],,])
	r01 = exp(c(ens[1,,]))
	eps = plogis(c(ens[6,,]))
	phase1 = 52*plogis(c(ens[7,,]))-26
	phase1 = ((phase1-t0+13) %% (52))*7 + 16800 #rescale phase1 to give max transmission rate rel. to. Dec31/Jan1 2016 
	gam1 = exp(c(ens[2,,]))
	rho1 = exp(c(ens[3,,]))
	phi1 = plogis(c(ens[5,,]))
	loglik = c(ens[11,,])

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

plotEpiPar.sis <- function(pathnm, save=F, filnm=F,wd='oct12/mcmcoutput/',t0=7800.286){
	load(paste0(wd,pathnm,'.Rdata')) #load in the mcmc output

	lmax=max(ens[dim(ens)[1],,])
	r01 = exp(c(ens[1,,]))
	eps = plogis(c(ens[5,,]))
	phase1 = 52*plogis(c(ens[6,,]))-26
	phase1 = ((phase1-t0+13) %% (52))*7 + 16800 #rescale phase1 to give max transmission rate rel. to. Dec31/Jan1 2016 
	gam1 = exp(c(ens[2,,]))
	#rho1 = exp(c(ens[3,,]))
	phi1 = plogis(c(ens[4,,]))
	loglik = c(ens[dim(ens)[1],,])

	if(save==T){jpeg(file=filnm,height=8,width=12,res=400, units='in')}
	par(mfrow=c(2,3))
	plot(r01,loglik,xlab=expression(R[0]), ylab='log-likelihood',ylim=c(lmax-20,lmax+5))
	plot(eps,loglik,xlab=expression(epsilon), ylab='log-likelihood',ylim=c(lmax-20,lmax+5))
	plot(as.Date(phase1),loglik,xlab=expression(tau[beta]~' [days rel. to Dec31/Jan1 2016]'), ylab='log-likelihood',ylim=c(lmax-20,lmax+5))
	plot(gam1^(-1),loglik,xlab=expression(gamma^{-1}~' [weeks]'), ylab='log-likelihood',ylim=c(lmax-20,lmax+5))
	#plot(rho1^(-1),loglik,xlab=expression(rho^{-1}~' [weeks]'), ylab='log-likelihood',ylim=c(lmax-20,lmax+5))
	plot(phi1,loglik,xlab=expression(phi), ylab='log-likelihood',ylim=c(lmax-20,lmax+5))
	mtext(side=3,outer=TRUE,line=-3,text=pathnm)	
	if(save==T){dev.off()}
}

plotEpiPar.constprevsis <- function(pathnm, save=F, filnm=F, wd='oct12/mcmcoutput/',t0=7800.286,burn=1000){
	load(paste0(wd,pathnm,'.Rdata')) #load in the mcmc output
	ens=ens[,1000:dim(ens)[2],]
	lmax=max(ens[dim(ens)[1],,])
	r01 = exp(c(ens[1,,]))
	#eps = plogis(c(ens[5,,]))
	phase1 = 52*plogis(c(ens[5,,]))-26
	phase1 = ((phase1-t0+13) %% (52))*7 + 16800 #rescale phase1 to give max transmission rate rel. to. Dec31/Jan1 2016 
	gam1 = exp(c(ens[2,,]))
	#rho1 = exp(c(ens[3,,]))
	phi1 = c(ens[4,,])
	loglik = c(ens[dim(ens)[1],,])

	if(save==T){jpeg(file=filnm,height=8,width=12,res=400, units='in')}
	par(mfrow=c(1,2))
	plot(r01,loglik,xlab=expression(R[0]), ylab='log-likelihood') #,ylim=c(lmax-20,lmax+5))
	#plot(eps,loglik,xlab=expression(epsilon), ylab='log-likelihood',ylim=c(lmax-20,lmax+5))
	#plot(as.Date(phase1),loglik,xlab=expression(tau[beta]~' [days rel. to Dec31/Jan1 2016]'), ylab='log-likelihood',ylim=c(lmax-20,lmax+5))
	#plot(gam1^(-1),loglik,xlab=expression(gamma^{-1}~' [weeks]'), ylab='log-likelihood',ylim=c(lmax-20,lmax+5))
	#plot(rho1^(-1),loglik,xlab=expression(rho^{-1}~' [weeks]'), ylab='log-likelihood',ylim=c(lmax-20,lmax+5))
	plot(phi1,loglik,xlab=expression(logit(phi)), ylab='log-likelihood') #,ylim=c(lmax-20,lmax+5))
	mtext(side=3,outer=TRUE,line=-3,text=pathnm)	
	if(save==T){dev.off()}
}

plotPsiPar <- function(pathnm, save=F, filnm=F,wd='oct12/mcmcoutput/',t0=7800.286){
	load(paste0(wd,pathnm,'.Rdata')) #load in the mcmc output
	tmind = dim(ens)[2]
	ens = ens[,seq(tmind-500,tmind,by=1),]
	lmax=max(ens[dim(ens)[1],,])
	psi0 = plogis(c(ens[4,,]))
	psi1 = plogis(c(ens[10,,]))
	epspsi = plogis(c(ens[9,,]))
	phasepsi = 52*plogis(ens[8,,])+18
	phasepsi = ((phasepsi-t0+13) %% (52))*7 + 16800 #rescale phasepsi to give max transmission rate rel. to. Dec31/Jan1 2016 
	loglik = c(ens[11,,])

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

plotPsiPar.sis <- function(pathnm, save=F, filnm=F,wd='oct12/mcmcoutput/',t0=7800.286){
	load(paste0(wd,pathnm,'.Rdata')) #load in the mcmc output

	lmax=max(ens[dim(ens)[1],,])
	psi0 = plogis(c(ens[3,,]))
	psi1 = plogis(c(ens[9,,]))
	epspsi = plogis(c(ens[8,,]))
	phasepsi = 52*plogis(ens[7,,])+18
	phasepsi = ((phasepsi-t0+13) %% (52))*7 + 16800 #rescale phasepsi to give max transmission rate rel. to. Dec31/Jan1 2016 
	loglik = c(ens[dim(ens)[1],,])

	if(save==T){jpeg(file=filnm,height=8,width=8,res=400,units='in')}
	par(mfrow=c(2,2))
	plot(psi0,loglik,xlab=expression(psi[0]), ylab='log-likelihood') #,ylim=c(lmax-20,lmax+5))
	plot(psi1,loglik,xlab=expression(psi[1]~' [per year]'), ylab='log-likelihood') #,ylim=c(lmax-20,lmax+5))
	plot(epspsi,loglik,xlab=expression(epsilon[psi]), ylab='log-likelihood') #,ylim=c(lmax-20,lmax+5))
	plot(as.Date(phasepsi),loglik,xlab=expression(tau[psi]~' [days rel. to Dec31/Jan1 2016]'), ylab='log-likelihood') #,ylim=c(lmax-20,lmax+5))
	mtext(side=3,outer=TRUE,line=-3,text=pathnm)	
	if(save==T){dev.off()}
}

plotPsiPar.constprevsis <- function(pathnm, save=F, filnm=F,wd='oct12/mcmcoutput/',t0=7800.286,burn=1000){
	load(paste0(wd,pathnm,'.Rdata')) #load in the mcmc output
	ens=ens[,1000:dim(ens)[2],]
	lmax=max(ens[dim(ens)[1],,])
	psi0 = c(ens[3,,])
	psi1 = c(ens[8,,])
	epspsi = plogis(c(ens[7,,]))
	phasepsi = 52*plogis(ens[6,,])+18
	phasepsi = ((phasepsi-t0+13) %% (52))*7 + 16800 #rescale phasepsi to give max transmission rate rel. to. Dec31/Jan1 2016 
	loglik = c(ens[dim(ens)[1],,])

	if(save==T){jpeg(file=filnm,height=8,width=8,res=400,units='in')}
	par(mfrow=c(2,2))
	plot(psi0,loglik,xlab=expression(logit(psi[0])), ylab='log-likelihood') #,ylim=c(lmax-20,lmax+5))
	plot(psi1,loglik,xlab=expression(logit(psi[1])), ylab='log-likelihood') #,ylim=c(lmax-20,lmax+5))
	plot(epspsi,loglik,xlab=expression(epsilon[psi]), ylab='log-likelihood') #,ylim=c(lmax-20,lmax+5))
	plot(as.Date(phasepsi),loglik,xlab=expression(tau[psi]), ylab='log-likelihood') #,ylim=c(lmax-20,lmax+5))
	mtext(side=3,outer=TRUE,line=-3,text=pathnm)	
	if(save==T){dev.off()}
}


plotSimPrev <- function(theta,pathnm,start=16801,end=18199,filnm=F,save=F){
	sim=getSimPrev1(theta)
	sim$t = sim$t - sim$t[1]
	tms = seq(start,end,by=1)
	if(save==T){pdf(file=filnm,height=4,width=4)}
	par(mfrow=c(1,1))
	plot(log10(s1)~as.Date(tms),sim,type='l',lty='dashed',ylim=c(4,9),
		main=paste0(pathnm,' Estimated Prevalence'),
		ylab=expression(log[10]~'Population'),
		xlab='Time [years]')
	lines(log10(r1)~as.Date(tms),sim,lty='dotted')
	lines(log10(i1)~as.Date(tms),sim)
	legend('bottomleft',legend=c('S','I','R'),lty=c(2,1,3))
}

plotSimPrev.constprevsis <- function(theta,pathnm) plotSimPrev(c(theta[1:2],20,theta[3:4],-100,theta[5:8]),pathnm)

#plot a collection of trajectories from the posterior, and uncertainty in inital conditions:
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
	sim=getSimPrev1(theta,tEnd=tEnd,nrows=length(tms),Y=Y)
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
		sim=getSimPrev1(theta,tEnd=length(tms)/7,nrows=length(tms),Y=yinit)
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
	ens = ens[,seq(niter-500,niter,by=1),]
	logr01 = c(ens[1,,])
	loggam1 = c(ens[2,,])
	logrho1 = c(ens[3,,])
	thetas = cbind(logr01,loggam1,logrho1)
	colnames(thetas) <- c('logr01','loggam1','logrho1')
	omegas = apply(thetas, 1, getOmega)
	return(omegas)
}

plotPdDist <- function(pathnm,save=F,filnm=F,wd='oct12/mcmcoutput/'){
	load(paste0(wd,pathnm,'.Rdata')) #load in the mcmc output	
	omegas <- getOmegaDist(ens)
	if(save==T){pdf(filnm,height=6,width=6)}
	par(mfrow=c(1,1))	
	hist(2*pi/omegas,xlab='Characteristic period of oscillation [weeks]',
		main=pathnm,freq=F)
	if(save==T){dev.off()}

}

goodpnms = setdiff(goodpnms, c('Influenza.A.H1','Mycoplasma.pneumoniae','Chlamydia.pneumoniae','B..pertussis')) #couldn't fit IAVH1, too few data

coronanms = c('Coronavirus.all',paste0('Coronavirus.',c('OC43','NL63','HKU1','229E'))) 
iavnms = c('Influenza.A.all','Influenza.A.H1.2009','Influenza.A.H3','Influenza.A..no.subtype','Influenza.B')
pivnms = paste0('Parainfluenza.Virus.',c(1:4))
othernms = c('RSV','Human.Metapneumovirus','Adenovirus','HRV.EV')

getParDist <- function(par,wd='oct12/mcmcoutput/',pathnms=goodpnms){
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

plotParDist <- function(pardist,labnms,save=F,filnm=F,parnm=F){
	if(save==T){pdf(filnm,height=4,width=4)}
	par(mfrow=c(1,1),mar=c(4,4,2,2))
	colnames(pardist) <- labnms	
	boxplot(pardist,las=2,main=parnm,cex.axis=1)
	#axis(side=1,labels=F)
	#text(x=1:5,y=par('usr')[3]-0.08,
	#	xpd=NA,srt=35,adj=0.965,cex=1.2,
	#	labels=labnms)
	if(save==T){dev.off()}
}

#write a function that extracts mles from the single-pathogen models and returns
#a vector for use in the coinfection models

getcoinfpar <- function(pathnm1,pathnm2){
	load('oct12/mcmcoutput/mles.Rdata')
	th1 = mles[pathnm1,]
	th2 = mles[pathnm2,]

	thetacoinf = c(th1[1],th2[1],th1[2],th2[2],th1[3],th2[3],mean(c(th1[4],th2[4])),
		th1[5],th2[5],0,0,mean(c(th1[6],th2[6])),th1[7],th2[7],mean(c(th1[8],th2[8])),
		mean(c(th1[9],th2[9])),mean(c(th1[10],th2[10])))
	return(thetacoinf)

}

plotFitCoinf <- function(pathnm1,pathnm2,tEnd=8000){
	theta = getcoinfpar(pathnm1,pathnm2)
	xpred = getSimDat(theta,tEnd=tEnd)
	xobs = getddat(pathnm1,pathnm2) 
	par(mfrow=c(2,2))
	plot(as.Date(16801:18199),xobs[,1],ylab='Coinfections',xlab='Time [years]',
		ylim=c(0,2*max(xobs[,1])))
	points(as.Date(16801:18199),xpred[,1],col='red')
	plot(as.Date(16801:18199),xobs[,3],ylab=paste0(pathnm1,' only'),xlab='Time [years]')
	points(as.Date(16801:18199),xpred[,3],col='red')
	plot(as.Date(16801:18199),xobs[,4],ylab=paste0(pathnm2,' only'),xlab='Time [years]')
	points(as.Date(16801:18199),xpred[,4],col='red')
	plot(as.Date(16801:18199),xobs[,2],ylab='Negative for both',xlab='Time [years]')
	points(as.Date(16801:18199),xpred[,2],col='red')

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
	sim=getSimPrev1(theta,tEnd=tEnd,nrows=length(tms),Y=Y)
	#tms = seq(start,end,by=1)
	#sim=tail(sim,n=length(tms))
	if(save==T){pdf(file=filnm,height=6,width=6)}
	par(mfrow=c(1,1))
	plot(log10(s1)~as.Date(tms),sim,type='l',lty='dashed',ylim=c(4,9),xlim=c(as.Date(tms[1]),as.Date(tail(tms,n=1)+xmax)),
		main=paste0(pathnm,' Estimated and Projected Prevalence'),
		ylab=expression(log[10]~'Population'),
		xlab='Time [years]',col='black')
	lines(log10(r1)~as.Date(tms),sim,lty='dotted',col='black')
	lines(log10(i1)~as.Date(tms),sim,col='black')
	legend('bottomleft',legend=c('S','I','R'),lty=c(2,1,3),col=c('black','black','black'))

	# save the initial condition to jitter as well:
	y=sim[dim(sim)[1],]
	y['S']=y['s1']
	y['I1']=y['i1']
	y['R1']=y['r1']
	y[c('I2','I12','R2','R1I2','R2I1','R1R2')] = 0
	y=y[2:10]

	inew = unlist(y['I1']) * c(1,0.5,0.2,0.05)
	pnew = unlist(y['S']/(y['S']+y['R1']))
	pnew = seq(pnew,pnew+(1-pnew)/2,length=3)
	newinits = data.frame(expand.grid(inew,pnew))
	colnames(newinits) <- c('i','p')
	newinits$s <- newinits$p * (327e6 - newinits$i)
	newinits$r <- 327e6 - newinits$i - newinits$s
	newinits <- newinits[,c('s','i','r')]

	switch <- rep(0,dim(newinits)[1])
	for(i in 1:dim(newinits)[1]){
		#theta=ensemble[,sample( seq(dim(ensemble)[2]-100, dim(ensemble)[2]), 1) ,sample(dim(ensemble)[3],1)]	
		#theta = ci99[,sample(1:dim(ci99)[2],1),sample(1:dim(ci99)[3],1)]	
		# comment this line out to draw parameters from the 99% confidence region boundary:	
		if(b==0){theta=getmcmcmle(ensemble)}

		#yinit = as.numeric(y)
		#jitter the initial condition - comment out to just visualize uncertainty in parameters
		#yinit[c(1,2,5)] = jitter(yinit[c(1,2,5)],amount=a*yinit[2])
		yinit = y
		yinit[c(1,2,5)] <- newinits[i,]
		#just run the simulation from Jan 1, 2016 to the end of the data fitting period (encoded in tms), 
		# which counts the number of weeks (but solver operates on timescale of days)
		sim1=getSimPrev1(theta,tInit = tail(sim$t,n=1),tEnd=tail(sim$t,n=1)+simmax,nrows=1000000,Y=as.numeric(yinit))
		tms2 = c(end:(end+dim(sim1)[1]-1))	
		#sim=tail(sim,n=length(tms))
		# color transients according to their attractor
	
		if(i==1){att=sim1$i1}		
		dist  = sum((tail(att,n=1000)-tail(sim1$i1,n=1000))^2)/327e6^2

		if(dist>thresh){cols='blue';switch[i]=1}else{cols = 'red'}
		if(i==1){cols='black'} 
		lines(log10(s1)~as.Date(tms2),sim1,lty='dashed',col=cols)				
		lines(log10(r1)~as.Date(tms2),sim1,lty='dotted',col=cols)
		lines(log10(i1)~as.Date(tms2),sim1,col=cols)	
	}
	if(save==T){dev.off()}
	return(list(switch,newinits))
}
