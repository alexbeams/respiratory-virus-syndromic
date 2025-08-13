# Codes for the manuscript: Estimating respiratory virus parameters from syndromic data: Using mathematical models to correct biases

The files in this repository reproduce benchmarking results (Figure 2 of the main text) and the two-virus simulation results from the Supplementary material. It also includes code to generate Figure 3 of the main text and the rest of the supplementary plots, but it is necessary to acquire data from BioMérieux to generate results (we cannot post the data publicly). Reasonable requests for access must be made to BioMérieux (ML-SyndromicTrendsDataRequest@biomerieux.com). We have included the codes to run the analysis for the version of the dataset they provided to us.


All codes are developed for R (version 4.4.2). We ran our analysis on Linux (Ubuntu 20.04.6 LTS (GNU/Linux 5.4.0-216-generic x86_64)) and MacOS (Sequoia 15.6). Slight modifications are needed for Windows (see just below).

Required R packages:
  - deSolve (simulates ODEs)
  - Rcpp (compiling and solving differential equations in C/C++)
  - parallel (allows for parallel computation)
  - xtable (display parameter tables in LaTeX)
  - scales (transparency in plots)

Installation of Rcpp may take several minutes. On Mac, it is convenient to run

brew install Rcpp

The other packages are easily installed in R using the install.packages command (e.g. install.packages('deSolve') )

# Note for Windows users:
The code model_functions.R contains the model definition for the differential equations. We solve the ODEs in C/C++, and this invovles creating a compiled object. The C file is sirmod/sir.c. On Linux and MacOS, the compiled object will be sir.so ("shared object"). On Windows, it is sir.dll ("dynamically loaded library"). The code is currently written to work for Linux/MacOS, but you will need to manually change the file extension from .so to .dll in model_functions.R (the R file is annotated so hopefully it is obvious where the edit needs to be made).

## General overview of codes in main:
You should begin by opening and running vignette.R. If you set your working directory to main, running source('vignette.R'), or running line-by-line, should do the following:
  1. produce a simulated dataset
  2. fit the model to that dataset using the ensemble MCMC method
  3. plot profile likelihoods of the ensemble MCMC run.

Codes in main produce results in Figures 2 and 3 of the main text. There are also three subdirectories taht contain more specific codes (Figure 2, Figure 3, and sirmod). The sirmod directory contains C/C++ files for rapidly simulating the differential equations (the main one being sir.c). The subdirectories called Figure2 and Figure3 contain some more specific codes useful for creating those figures.

  - vignette.R is relatively quick to run (<5 minutes using default settings). It simulates a dataset, and fits the SIRS model to that dataset using an ensemble MCMC sampling method (100 timesteps by default, which often produces clear profile loglikelihoods when initializing near the true values used in the simulations). Mainly this is to visualize the profile loglikelihoods, which reveals the amount of information available about parameters. Simply running source('vignette.R') in the termal should be sufficient, as it invokes all of the other codes upon which it depends.
    
  - model_functions.R: running this will load a number of useful functions (getSim, getSimDat, getSimMean, getSimPrev, getLogLik) that are the workhorse of the methods
    
  - plotfunctions.R: contains a number of useful functions for displaying results, notably plotEpiPar and plotPsiPar, which display profile loglikelihoods of ensemble MCMC output
    
  - ensemble_sampler.R: this defines the affine-invariant sampler that uses the ``stretch move'' described in Goodman and Weare (2010). Parameters of SIRS models have high posterior correlations, so an affine-invariant sampler is very useful.
    
  - getFigure2.R will reproduce Figure 2 of the main text. With default settings, it should run in less than 10 minutes. It simulates 180 datasets for 6 different parameter combinations (so, 30 simulations for each individual parameter set by default), and fits models using an ensemble MCMC sampling method. By default, this only runs for 5 timesteps so you can verify it works on your machine (but we use 2000 timesteps to generate the results for Figure 2 in the main text). The functions plotEpiPar and plotPsiPar defined in plotfunctions.R are useful for assessing convergence of the ensemble MCMC method (try increasing timesteps and ensemble size, ideally on a server that you can walk away from for an afternoon/overnight!). Obviously you should use plotEpiPar and plotPsiPar to visualize the MCMC output that gets produced in Figure2/simdatfits/ to ensure that Maximum Likelihood Estimates have been found before creating Figure 2.

other subdirectories in main:
  - Figure 2/: this folder contains some codes to generate simulations and model fits for the benchmarking results displayed in Figure 2 (benchmarking_sim.R and benchmarking_est.R). It should also contains folders called simdats and simdatfits that will contain your simulations and ensemble MCMC runs, respectively (but you will need to create these sub-directories on your machine). Invoke getFigure2.R produces the simulations and MCMC ensemble runs.
  - Figure 3/: this folder contains code to generate Figure 3, but will not run correctly unless you obtain the data from BioMérieux. It mostly calls the same functions as vignette and getFigure2.R (contained in model_functions.R, plotfunctions.R, and ensemble_sampler.R). Note: if you do obtain data, we recommend manually selecting a reasonable guess for parameters instead of initializing the ensemble MCMC method from random parameters. It is also useful to run the ensemble MCMC sampler for about 2000 timesteps, re-initialize near the current optimum, and run for another 2000 timesteps, repeating as necessary until convergence is achieved (or just start from the values we report in the main text).

## General overview of codes in supplement:
  - 2virusModels/: this directory contains the code to generate Figures S25-27, which describe results of fitting single-pathogen SIRS models to simulated data produced by 2-virus models. If you set your working directory in 2virusModels, running get2VirusSims.R will generate the simulations with two cocirculating viruses, and then running fitModels.R will fit single-virus SIRS models to aggregated detections of those simulations. It should produce figures analagous to Figures S25-27 of the supplement (they will differ slightly because of stochasticity).
  - ModelFiles directory: this conains C/C++ code to simulate a 9-state, two-virus SIRS model
  - sirmod directory: this is an exact copy of the sirmod subdirectory in main; needed for fitting single-virus models to two-virus simulations
  - get2VirusSims.R - this simulates the two-virus model and aggregates detections; stores simulated counts (simdats) and aggregated detections (dfs) in the 2virussims subdirectory
  - fitModels.R - this fits single-virus SIRS models to the simulated counts from the two-virus models
  - extra copies of model_functions.R, plotfunctions.R, ensemble_sampler.R for convenience
  - Note: the supplmental figures that display profile loglikelihoods are visualized using the plotEpiPar and plotPsiPar functions from plotfunctions.R for the ensemble MCMC output obtained by fitting models to the BioMéreiux data. 
