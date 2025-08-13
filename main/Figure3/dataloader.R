
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


source('model_functions.R')
source('plotfunctinos.R')


pnms <- names(dat)[7:28]
goodpnms <- pnms[1:20] #B..parapertussis and SARS.CoV.2 have lots of NAs




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


# extract the MLE from the mcmc ensemble

getmcmcmle <- function(ens){
        loglikind = dim(ens)[1]
        lmax = max(ens[loglikind,,])
        ind1 = which(apply(ens[loglikind,,], 1, function(x) sum(which(x==lmax)) )>0)[1]
        ind2 = which(ens[loglikind,ind1,] == lmax)[1]
        mcmcmle = ens[,ind1,ind2]
        return(mcmcmle)
}



goodpnms = setdiff(goodpnms, c('Influenza.A.H1','Mycoplasma.pneumoniae','Chlamydia.pneumoniae','B..pertussis')) #couldn't fit IAVH1, too few data

coronanms = c('Coronavirus.all',paste0('Coronavirus.',c('OC43','NL63','HKU1','229E')))
iavnms = c('Influenza.A.all','Influenza.A.H1.2009','Influenza.A.H3','Influenza.A..no.subtype','Influenza.B')
pivnms = paste0('Parainfluenza.Virus.',c(1:4))
othernms = c('RSV','Human.Metapneumovirus','Adenovirus','HRV.EV')




