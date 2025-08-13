# Codes for the manuscript: Estimating respiratory virus parameters from syndromic data: Using mathematical models to correct biases

The files in this repository reproduce benchmarking results (Figure 2 of the main text) and the two-virus simulation results from the Supplementary material. It also includes code to generate Figure 3 of the main text, but it is necessary to acquire data from BioMérieux to generate results. 

All codes are developed for R (version 4.4.2). We ran our analysis on Linux and MacOS. Slight modifications are needed for Windows.

The data for the main analysis are not publicly available. Reasonable request for access must be made to BioMérieux. We have included the codes to run the analysis for the version of the dataset they provided to us.

Required R packages:
  - deSolve (simulates ODEs)
  - Rcpp (compiling and solving differential equations in C/C++)
  - parallel (allows for parallel computation)

# Note for Windows users:
The code model_functions.R contains the model definition for the differential equations. We solve the ODEs in C/C++, and this invovles creating a compiled object. The C file is sirmod/sir.c. On Linux and MacOS, the compiled object will be sir.so ("shared object"). On Windows, it is sir.dll ("dynamically loaded library"). The code is currently written to work for Linux/MacOS, but you will need to manually change the file extension from .so to .dll.

## General overview of codes:
We have also distributed code in two folders: main, and supplement. These contain more specific codes for generating figures in those sections of the manuscript. There is also a folder called sirmod/ that contains C/C++ files for rapidly simulating the differential equations (the main one being sir.c).

  - vignette.R is relatively quick to run. It simulates a dataset, and fits the SIRS model to that dataset using an ensemble MCMC sampling method (100 timesteps by default, which often produces clear profile loglikelihoods when initializing near the true values used in the simulations). Mainly this is to visualize the profile loglikelihoods, which reveals the amount of information available about parameters. Simply running source('vignette.R') in the termal should be sufficient, as it invokes all of the other codes upon which it depends.
  - model_functions.R: running this will load a number of useful functions (getSim, getSimDat, getSimMean, getSimPrev, getLogLik) that are the workhorse of the methods
  - plotfunctions.R: contains a number of useful functions for displaying results, notably plotEpiPar and plotPsiPar, which display profile loglikelihoods of ensemble MCMC output
  - ensemble_sampler.R: this defines the affine-invariant sampler that uses the ``stretch move'' described in Goodman and Weare (2010). Parameters of SIRS models have high posterior correlations, so an affine-invariant sampler is very useful.
  - getFigure2.R will reproduce Figure 2 of the main text. It simulates 180 datasets for 6 different parameter combinations (so, 30 simulations for each individual parameter set by default), and fits models using an ensemble MCMC sampling method. By default, this only runs for 5 timesteps so you can verify it works on your machine (but we use 2000 timesteps to generate the results for Figure 2 in the main text). The functions plotEpiPar and plotPsiPar defined in plotfunctions.R are useful for assessing convergence of the ensemble MCMC method (try increasing timesteps and ensemble size, ideally on a server that you can walk away from for an afternoon/overnight!). Obviously you should use plotEpiPar and plotPsiPar to visualize the MCMC output that gets produced in Figure2/simdatfits/ to ensure that Maximum Likelihood Estimates have been found before creating Figure 2.

main/:
  - Figure 2/: this folder contains some codes to generate simulations and model fits for the benchmarking results displayed in Figure 2 (benchmarking_sim.R and benchmarking_est.R). It also contains folders called simdats and simdatfits that will contain your simulations and ensemble MCMC runs, respectively. These get produced when you invoke getFigure2.R.
  - Figure 3/: this folder contains code to generate Figure 3, but will not run correctly unless you obtain the data from BioMérieux. It mostly calls the same functions as vignette and getFigure2.R (contained in model_functions.R, plotfunctions.R, and ensemble_sampler.R). Note: if you do obtain data, we recommend manually selecting a reasonable guess for parameters instead of initializing the ensemble MCMC method from random parameters. It is also useful to run the ensemble MCMC sampler for about 2000 timesteps, re-initialize near the current optimum, and run for another 2000 timesteps, repeating as necessary until convergence is achieved (or just start from the values we report in the main text).

supplement/:
  - 2virusModels/: this folder contains the code to generate Figures S25-27, which describe results of fitting single-pathogen SIRS models to simulated data produced by 2-virus models
  - Note: the supplmental figures that display profile loglikelihoods are visualized using plotEpiPar and plotPsiPar for the ensemble MCMC output obtained by fitting models to the BioMéreiux data. 
