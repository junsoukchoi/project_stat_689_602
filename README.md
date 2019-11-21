ZIPBN: Zero-Inflated Poisson Bayesian Networks for Zero-Inflated Count
Data
================
Junsouk Choi

  - [Installation](#installation)

The R package ‘ZIPBN’ implements a Markov chain Monte Carlo (MCMC) to
fit zero-inflated poisson Bayesian networks for zero-inflated count data
such as scRNA-seq data. It is based on my current reseach on scalable
Bayesian networks (also known as directed acyclic graph) for scRNA-seq
data which comprise of counts with excessive zeros. A Bayesian network
factorizes the joint distribution of random variables into a set of
local distributions. To deal with an excess of zero counts in the data,
we model each local distribution to be a zero-inflated poisson model. We
further formulate a hierarchical model in Bayesian perspective by
specifing spike-and-slab priors for edge selection to determine graph
structure. The MCMC algorithm basically uses Metropolis-Within-Gibbs
sampling to sample parameters from the posterior distributions.

## Installation

To install the latest version from Github, use

``` r
library(devtools)
devtools::install_github("junsoukchoi/ZIPBN")
```
