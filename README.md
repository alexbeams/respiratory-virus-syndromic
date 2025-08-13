# Codes for the manuscript: Estimating respiratory virus parameters from syndromic data: Using mathematical models to correct biases

The files in this repository reproduce benchmarking results (Figure 3 of the main text) and the two-virus simulation results from the Supplementary material. All codes are in R (version 4.4.2).

The data for the main analysis are not publicly available. Reasonable request for access must be made to BioMérieux. We have included the codes to run the analysis for the version of the dataset they provided to us.

Required R packages:
-deSolve
-parallel

main:
  - vignette.R is relatively quick to run. It simulates a dataset, and fits the SIRS model to that dataset using an ensemble MCMC sampling method (100 timesteps by default, which often produces clear profile loglikelihoods when initializing near the true values used in the simulations).
  - getFigure2.R will reproduce Figure 2 of the main text. It simulates 180 datasets (by default), and fits models using an ensemble MCMC sampling method. By default, this only runs for 5 timesteps (but we used 2000 timesteps to generate the results for Figure 2 in the main text). The functions plotEpiPar and plotPsiPar defined in plotfunctions.R are useful for assessing convergence of the ensemble MCMC method (the solution is usually to increase timesteps and ensemble size)
  - Figure 2: this folder contains the code to generate Figure 2 (benchmarking results using simulations)
  - Figure 3: this folder contains code to generate Figure 3, but will not run correctly unless you obtain the data from BioMérieux

supplement:
  - 2virusModels: this folder contains the code to generate Figures S25-27, which describe results of fitting single-pathogen SIRS models to simulated data produced by 2-virus models
