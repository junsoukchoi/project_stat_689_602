# evaluate log-likelihood of each observation for the j-th component of ZIPBN with logit(pi) and log(lambda)
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
   
   # do MCMC iterations
   for (t in 1 : n.samples)
   {
      # update alpha (Metropolis-Hastings step)
      # update beta (Metropolis-Hastings step)
      # update delta (Metropolis-Hastings step)
      # update gamma (Metropolis-Hastings step)
      # update A (Metropolis-Hastings step)
      # update tau (Gibbs sampling)
      # update rho (Gibbs sampling)
      
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