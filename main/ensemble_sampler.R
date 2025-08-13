#MCMC scheme. There is a new-ish package called mcmcensemble, but at time of analysis
# we use the following code that uses the stretch move method from Goodman and Weare (2010)

# user supplies the following:
#	loglikelihood - function of parameters, theta
#	thetaCurrent - a realization of theta
#	dimTime - number of timesteps
#	dimEnsemble - number of ensemble members/particles

getEnsemble <- function(loglikelihood,thetaCurrent,dimTime,dimEnsemble){
        dimTheta = length(thetaCurrent)

        ensemble <- array( dim=c(dimTheta+1,dimTime,dimEnsemble) )

        #initialize the ensemble - Gaussian around the initial point
        for(i in 1:dimEnsemble){
                likvl <- -Inf

                while(likvl == -Inf ){
                prmvl <- rnorm(dimTheta,thetaCurrent,1e-2)
                likvl <- loglikelihood(prmvl)
        }

        ensemble[1:dimTheta,1,i] <- prmvl
        ensemble[dimTheta+1,1,i] <- likvl
        }
        #need to not have loglikelihoods = -Inf to start


        #generate all the independent stretches initially
        y <- runif(2*dimTime*dimEnsemble)
        a <- 10 #tuning parameter; stretches at most by a factor of a or compresses down by a factor of 1/a
        Z <- (y*(a-1)+1)^2/a #normalized the g(z) from Goodman and Weare (2010), generated unif(0,1)'s and inverted to get stretches
        Zind <- 1

        for(i in 1:(dimTime-1)){
         print(i)
         for(k in 1:dimEnsemble){
          Xk = ensemble[1:dimTheta,i,k]
          loglikelihood_Current = ensemble[1+dimTheta,i,k]
          j = sample(setdiff(1:dimEnsemble,k),1)
          Xj = ensemble[1:dimTheta,i,j]
          Y <<- Xj + (Xk-Xj)*Z[Zind]
	  loglikelihood_Proposed = loglikelihood(Y)
          if(is.na(loglikelihood_Proposed)){loglikelihood_Proposed=-Inf}
          proposal_Probability = Z[Zind]*exp( loglikelihood_Proposed - loglikelihood_Current ) #metropolis-hastings acceptance probability from Goodman and Weare (2010)
          if(runif(1)<proposal_Probability){
           ensemble[,i+1,k] <- c(Y,loglikelihood_Proposed)
          }else{
           ensemble[,i+1,k] <- c(Xk,loglikelihood_Current)
          }
          Zind <- Zind + 1
         }
        }
        return(ensemble)
}

