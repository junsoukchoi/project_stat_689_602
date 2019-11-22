#' Function for evaluating log-likelihood of each observation for the j-th component of ZIPBN
#'
#' @param x_j data for the j-th variable (node)
#' @param logitPi logit(pi)
#' @param logLambda log(lambda)
#'
#' @return log-likelihood of each observation for the j-th component of ZIPBN 
#' @export
#'
#' @examples
llik_ZIPBN_j = function(x_j, logitPi, logLambda)
{
   # calculate pi and lambda
   pi     = exp(logitPi) / (1 + exp(logitPi))
   lambda = exp(logLambda)
   
   # evaluate and return log-likelihood of each observation
   llik   = log(pi + (1 - pi) * exp(-lambda)) * (x_j == 0) + 
      (log(1 - pi) + dpois(x_j, lambda, log = TRUE)) * (x_j > 0)
   return(llik)
}

# MCMC algorithm to sample parameters from posterior distributions for ZIPBN models
mcmc_ZIPBN = function(x, starting, tuning, priors, n.samples)
{
   # store the sample size and the number of variables (nodes)
   n = nrow(x)
   p = ncol(x)
   
   # initialize parameters with starting values supplied
   alpha = starting$alpha
   beta  = starting$beta
   delta = starting$delta
   gamma = starting$gamma
   tau   = starting$tau
   rho   = starting$rho
   
   # set precisions with supplied values for the Metropolis sampler Normal proposal distributions 
   phi_alpha = tuning$phi_alpha
   phi_beta  = tuning$phi_beta
   phi_delta = tuning$phi_delta
   phi_gamma = tuning$phi_gamma
   
   # set hyperparameters with supplied values
   nu = priors$nu
   b  = priors$b
   c  = priors$c
   
   # initialize MCMC samples
   A_MCMC     = array(NA, dim = c(p, p, n.samples))
   alpha_MCMC = array(NA, dim = c(p, p, n.samples))
   beta_MCMC  = array(NA, dim = c(p, p, n.samples))
   delta_MCMC = matrix(NA, p, n.samples)
   gamma_MCMC = matrix(NA, p, n.samples)
   tau_MCMC   = matrix(NA, 4, n.samples)
   rho_MCMC   = rep(NA, n.samples)
   
   # calculate logit(pi) and log(lambda)
   logitPi   = tcrossprod(x, alpha) + + matrix(delta, n, p, byrow = TRUE)
   logLambda = tcrossprod(x, beta) + matrix(gamma, n, p, byrow = TRUE)
   
   # do MCMC iterations
   for (t in 1 : n.samples)
   {
      # update alpha (Metropolis-Hastings step)
      # update beta (Metropolis-Hastings step)
      # update delta (Metropolis-Hastings step)
      # update gamma (Metropolis-Hastings step)
      # update A (Metropolis-Hastings step)
      # update tau (Gibbs sampling)
      tau[1] = rgamma(1, shape = b[1] + p * (p - 1) / 2, 
                      rate = c[1] + (nu * sum((A == 0) * alpha * alpha) + sum((A == 1) * alpha * alpha)) / 2)   # tau_alpha
      tau[2] = rgamma(1, shape = b[2] + p * (p - 1) / 2, 
                      rate = c[2] + (nu * sum((A == 0) * beta * beta) + sum((A == 1) * beta * beta)) / 2)   # tau_beta
      tau[3] = rgamma(1, shape = b[3] + p / 2, rate = c[3] + sum(delta * delta) / 2)   # tau_delta
      tau[4] = rgamma(1, shape = b[4] + p / 2, rate = c[4] + sum(gamma * gamma) / 2)   # tau_gamma
      
      # update rho (Gibbs sampling)
      rho = rbeta(1, shape1 = b[5] + sum(A == 1), shape2 = c[5] + sum(A == 0) - p)
      
      # store MCMC samples of iteration t
      alpha_MCMC[ , , t] = alpha
      beta_MCMC[ , , t]  = beta
      delta_MCMC[ , t]   = delta
      gamma_MCMC[ , t]   = gamma
      A_MCMC[ , , t]     = A
      tau_MCMC[ , t]     = tau
      rho_MCMC[t]        = rho
   }
   
   # return a list of MCMC samples
   return(list(alpha = alpha_MCMC,
               beta  = beta_MCMC,
               delta = delta_MCMC,
               gamma = gamma_MCMC,
               A     = A_MCMC,
               tau   = tau_MCMC,
               rho   = rho_MCMC))
}