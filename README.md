ZIPBN: Zero-Inflated Poisson Bayesian Networks for Zero-Inflated Count
Data
================
Junsouk Choi

  - [Installation](#installation)
  - [Usage](#usage)

The R package `ZIPBN` implements a Markov chain Monte Carlo (MCMC) to
fit zero-inflated poisson Bayesian networks for zero-inflated count data
such as scRNA-seq data. It is based on my current reseach on scalable
Bayesian networks (also known as directed acyclic graph) for scRNA-seq
data which comprise of counts with excessive zeros. A Bayesian network
factorizes the joint distribution of random variables into a set of
local distributions. To deal with an excess of zero counts in data, we
model each local distribution to be a zero-inflated poisson model. We
further formulate a hierarchical model in Bayesian perspective by
specifing spike-and-slab priors for edge selection to determine graph
structure. The MCMC algorithm basically uses Metropolis-Within-Gibbs
sampling to sample parameters from the posterior distributions.

## Installation

To install the latest version from Github, use

``` r
library(devtools)
devtools::install_github("junsoukchoi/ZIPBN", build_vignettes = TRUE)
```

## Usage

The following example describes how to use package `ZIPBN` to conduct
analysis on zero-inflated cound data. Function`mcmc_ZIPBN` takes
zero-inflated count data with some arguments needed for MCMC sampler,
and returns MCMC samples from posterior distributions of the parameters
defined by ZIPBN models. It also gives Metropolis sampling acceptance
rates, which indicates whether our MCMC sampler goes well. For more
information on `ZIPBN` package, please access the package documentation
or
[vignettes](https://github.com/junsoukchoi/ZIPBN/blob/master/vignettes/ZIPBN.pdf).
Please feel free to contact the author.

``` r
library(ZIPBN)


## Example data
set.seed(7)

# generate a simple graph: X1 -> X2 -> X3
p = 3
A = matrix(0, p, p)
A[3, 2] = A[2, 1] = 1

# parameters of the ZIPBN model, given graph A
alpha = matrix(0, p, p)
alpha[A == 1] = 0.3
beta = matrix(0, p, p)
beta[A == 1] = 0.2
delta = rep(1, p)
gamma = rep(1.5, p)

# generate data from the ZIPBN model
n = 200
x = matrix(0, n, p)
for (j in 1:p) {
    # calculate pi_j
    pi = exp(x %*% alpha[j, ] + delta[j])
    pi = pi/(1 + pi)
    # calculate mu_j
    mu = exp(x %*% beta[j, ] + gamma[j])
    # generate data for X_j
    x[, j] = rpois(n, mu) * (1 - rbinom(n, 1, pi))
}


## fit ZIPBN models create starting value list
m = colMeans(x)
v = apply(x, 2, var)
starting = list(alpha = matrix(0, p, p), beta = matrix(0, p, p), delta = log((v - 
    m)/(m * m)), gamma = log((v - m + m * m)/m), A = matrix(0, p, p), tau = c(10, 
    10, 1, 1), rho = 0.1)

# create tuning value list
tuning = list(phi_alpha = c(1e+08, 20), phi_beta = c(1e+08, 100), phi_delta = 5, 
    phi_gamma = 50, phi_A = c(1e+10, 10, 10, 1, 10))

# create priors list
priors = list(nu = 10000^2, tau_alpha = c(0.01, 0.01), tau_beta = c(0.01, 0.01), 
    tau_delta = c(0.01, 0.01), tau_gamma = c(0.01, 0.01), rho = c(0.5, 0.5))

# run mcmc_ZIPBN function
n_sample = 2000
n_burnin = 1000
out = mcmc_ZIPBN(x, starting, tuning, priors, n_sample, n_burnin)


## posterior inference via ZIPBN models report Metropolis sampling acceptance
## percents
out$acceptance

# recover garph structure
cutoff = 0.5
A_est = 1 * (apply(out$samples$A, c(1, 2), mean) > cutoff)

# calculate the posterior mean of each parameter, given the recovered graph
subset = apply(out$samples$A == array(A_est, dim = c(p, p, n_sample - n_burnin)), 
    3, all)
alpha_est = apply(out$samples$alpha[, , subset], c(1, 2), mean)
beta_est = apply(out$samples$beta[, , subset], c(1, 2), mean)
delta_est = rowMeans(out$samples$delta[, subset])
gamma_est = rowMeans(out$samples$gamma[, subset])

# report the posterior mean of each parameter with the recoverd graph
A_est
round(alpha_est, digits = 2)
round(beta_est, digits = 2)
round(delta_est, digits = 2)
round(gamma_est, digits = 2)
```
