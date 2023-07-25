# Supplemental Archive for "Phylogenies increase power to detect highly transmissible viral genome variants"

This repository contains code for recreating simulations and analyses from the manuscript "Phylogenies increase power to detect highly transmissible viral genome variants" by Michael R May and Bruce Rannala.

## Dependencies

The following `R` packages are required to run the code in this archive:
- `parallel`
- `ape`
- `phangorn`
- `TESS`
- `expm`
- `Rcpp`
- `RcppEigen`
- `matrixStats`
- `ape`
- `TreeTools`
- `RColorBrewer`
- `ggplot2`
- `viridis`
- `cubature`
- `dispRity`
- `Matrix`
- `deSolve`
- `pracma`
- `abind`

## Shared scripts

The directory `YuleSelectionArchive/scripts` are generic scripts that are used by subsequent simulations and analyses.

### Tree model

- `scripts/YuleLikelihoodBinary.R`: (biallelic phylogenetic model, for comparison against count model)
- `scripts/likelihood.R`: likelihood function for the phylogenetic model, and associated estimation functions

### Count model

- `scripts/countClass.R`: class for evaluating the likelihood function for the biallelic count model
- `scripts/countMatrixClass.R`: (deprecated) class for calculating transition probabilities with the biallelic count model
- `scripts/countMatrixClass2.R`: class for calculating transition probabilities with the biallelic count model
- `scripts/sampleMatrixClass.R`: class for calculating sampling densities with the biallelic count model
- `scripts/countClassBackward.R`: class for evaluating the likelihood function for the biallelic count model, uses backward algorithm to integrate root age
- `scripts/countMatrixClassBackward.R`: class for calculating transition probabilities with the biallelic count model, used with backward algorithm
- `scripts/sampleMatrixClassBackward.R`: class for calculating sampling densities with the biallelic count model, used with backward algorithm
- `scripts/computeDerivative.cpp`: C++ functions for evaluating ODE of biallelic count model

### Other function

- `scripts/AdaptiveSimpsonsIntegrator.R`: for numerical evaluation of posterior (marginal likelihood, posterior summaries)
- `scripts/utils.R`: various utility functions used by other scripts

## Simulation 1

In the first simulation, we compare the performance of the count and phylogenetic models for detecting transmission-enhancing variants as a function of time since its origin.
All of this code is intended to run from directory `YuleSelectionArchive/simulation_1`.
To simulate the data, run `Rscript scripts/simulate_data.R` from terminal (or run the Rscript however you want).
The simulated data will be generated in the `sims` folder, with the structure `sims/factor_X/sim_i/day_j` where `X` is the effect size, `i` is the simulation replicate, and `j` is the data sampled on day `j`.
To analyze a dataset, use `Rscript scripts/analysis.R sims/factor_1.25/sim_1/day_3` to, e.g., analyze the first replicate sampled on day 3 with effect size 1.25.

### scripts

- `scripts/simulate_data.R`: simulates viral outbreaks from the moment the variant arises, and draws samples at fixed time points until the end of the simulation
- `scripts/simulate_outbreak.R`: simulates a single outcome, used by `scripts/simulate_data.R`
- `scripts/compute_quantiles.R`: computes fraction of simulations that will be rejected for being too large for each effect size
- `scripts/summarize_rejects.R`: summarizes the fraction of rejected simulations of various types
- `scripts/analysis.R`: performs an analysis of a given dataset

## Simulation 2

In the second simulation, we analyze the power to reject the neutral model, and identify the true transmission-enhancing site, using the phylogenetic model.
All of this code is intended to run from directory `YuleSelectionArchive/simulation_2`.
To simulate the data, run `Rscript scripts/simulate_data.R` from terminal (or run the Rscript however you want).
The simulated data will be generated in the `sims` folder, with the structure `sims/tips_X_f_Y/rep_i` where `X` is the number of samples, `Y` is the effect size (factor, f), and `i` is the simulation replicate.
To analyze a dataset, use `Rscript scripts/analysis.R sims/tips_100_f_1.5/rep_1` to, e.g., analyze the first replicate with 100 tips and effect size 1.5.

### scripts

- `scripts/simulate_data.R`: simulates viral outbreaks of a given size
- `scripts/simulate_outbreak.R`: simulates a single outcome, used by `scripts/simulate_data.R`
- `scripts/analysis.R`: performs an analysis of a given dataset

## Validation

All of this code is intended to run from directory `YuleSelectionArchive/validation`.

### ODE validation

Here, we validate our implementation of the count model ODEs.

- `scripts/validate_odes.R`: validation of ODEs of the biallelic count model by comparing against Monte Carlo simulations
- `scripts/validate_backward_odes.R`: validation of ODEs of the biallelic count model (using the backward algorithm) by comparing against Monte Carlo simulations (assumes you have simulated data for simulation 2 in the appropriate location).

### Count and phylogenetic model comparison

Here, we validate the phylogenetic and count likelihoods by comparing the integrated phylogenetic likelihood against the count model.
Integrating the phylogenetic likelihood over all tree topologies and node ages should yield a likelihood identical to the count model.
Therefore, this validation demonstrates that: 1) our theoretical derivation of the likelihoods is correct, and 2) our implementations of both models are correct.

- `scripts/compare_likelihood_2_samples.R`: for two samples, compare the phylogenetic likelihood integrated over node ages and topologies to the likelihood under the count model
- `scripts/compare_likelihood_3_samples.R`: for three samples, compare the phylogenetic likelihood integrated over node ages and topologies to the likelihood under the count model
- `scripts/compare_likelihood_4_samples.R`: for four samples, compare the phylogenetic likelihood integrated over node ages and topologies to the likelihood under the count model






















