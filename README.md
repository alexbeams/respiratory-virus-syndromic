# Codes for the manuscript: Estimating respiratory virus parameters from syndromic data: Using mathematical models to correct biases

The files in this repository reproduce benchmarking results (Figure 2 of the main text) and the two-virus simulation results from the Supplementary material. It also includes code to generate Figure 3 of the main text, but it is necessary to acquire data from BioMérieux to generate results. 

All codes are developed for R (version 4.4.2).

The data for the main analysis are not publicly available. Reasonable request for access must be made to BioMérieux. We have included the codes to run the analysis for the version of the dataset they provided to us.

Required R packages:
  -deSolve (simulates ODEs)
  -Rcpp (compiling and solving differential equations in C/C++)
  -parallel (allows for parallel computation)

main:
  - vignette.R is relatively quick to run. It simulates a dataset, and fits the SIRS model to that dataset using an ensemble MCMC sampling method (100 timesteps by default, which often produces clear profile loglikelihoods when initializing near the true values used in the simulations).
  - getFigure2.R will reproduce Figure 2 of the main text. It simulates 180 datasets for 6 different parameter combinations (so, 30 simulations for each individual parameter set by default), and fits models using an ensemble MCMC sampling method. By default, this only runs for 5 timesteps so you can verify it works on your machine (but we use 2000 timesteps to generate the results for Figure 2 in the main text). The functions plotEpiPar and plotPsiPar defined in plotfunctions.R are useful for assessing convergence of the ensemble MCMC method (try increasing timesteps and ensemble size!)
  - Figure 2/: this folder contains the code to generate simulations and model fits for the benchmarking results displayed in Figure 2
  - Figure 3/: this folder contains code to generate Figure 3, but will not run correctly unless you obtain the data from BioMérieux. It calls the same functions in main (contained in model_functions.R, plotfunctions.R, and ensemble_sampler.R).

supplement:
  - 2virusModels: this folder contains the code to generate Figures S25-27, which describe results of fitting single-pathogen SIRS models to simulated data produced by 2-virus models
  - Note: the supplmental figures that display profile loglikelihoods are visualized using plotEpiPar and plotPsiPar for the ensemble MCMC output obtained by fitting models to the BioMéreiux data. 
