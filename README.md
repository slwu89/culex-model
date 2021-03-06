# culex-model

A discretized stochastic and deterministic model of temperature and photoperiod driven culex dynamics with movement of adults based on the mathematical model presented in:

1. Paper: Uncovering mechanisms behind mosquito seasonality by integrating mathematical models and daily empirical population data: Culex pipiens in the UK [here](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-019-3321-2), Code repo: [here](https://github.com/davewi13/Mosquito-seasonality-paper)
2. Paper: Modelling the effect of temperature on the seasonal population dynamics of temperate mosquitoes [here](https://www.sciencedirect.com/science/article/pii/S0022519316300285), Code repo: [here](https://github.com/davewi13/Temperate-Mosquito-DDE)

Running this model requires [RcppArmadillo](https://dirk.eddelbuettel.com/code/rcpp.armadillo.html) to be installed on the users machine, which provides fast linear algebra. The model is run in discrete time, but the size of the time step (`dt` in code) may be arbitrarily small. Rates are converted into probabilities by an exponential CDF.

## Files

Within the directory "src_R" are the following files:
* example.R: set up and run a simulation with 2 patches, and some small level of movement of adults between each patch. In order to transform the backward-looking developmental delays into forward-looking delays, we first need to integrate ODEs describing the change in delay as a function of time (temperature).
* culex.hpp: contains the main definition of the stochastic and deterministic models.
* culex.cpp: the wrapper Rcpp interface for users to control the models from R.
* tau_ODE.R: the ODEs describing the delay lengths.
