# Companion Code

This repository contains the companion code for the preprint **“Efficient Bayesian Inference for Multiple Network Data.”**

The repository provides an example illustrating part of the simulation study presented in the article. The code is implemented in R and relies on the **BayesCER** package, which is included in this repository.

## Installation

Before running the code, please install the **BayesCER** package from the repository.

## Running the Example

Once the package has been installed, run the following script from the root directory of the repository:

```r
source("reproduction/simulation_study.R")
```

The script runs the simulation study described in the article. The complete simulation study involves a substantially larger number of simulations and may require considerable computational time.
Once the code has been executed, run the following script to reproduce the plots showed in the article.

```r
source("data-raw/plots.R")
```
A subset of the simulation results is precomputed and included with the package for illustration in the vignette. This avoids rerunning the simulations when building the vignette.

```r
data("sim_data_result", package = "BayesCER")
```
