#' Implementation of MCMC sampler for ZIPBN models
#'
#' @param x a matrix containing data
#' @param starting a list of parameters' starting values for MCMC
#' @param tuning a list of precision values for Metropolis sampler Normal proposal distribution
#' @param priors a list of hyperparameter values for priors
#' @param n_samples the number of MCMC iterations
#'
#' @return MCMC samples from posterior distributions of ZIPBN models
#' @export
#'
#' @examples
mcmc_ZIPBN = function(x, starting, tuning, priors, n_samples = 5000)
{
   # store the sample size and the number of variables (nodes)
   n = nrow(x)
   p = ncol(x)
   
   # initialize parameters with starting values supplied
   alpha = starting$alpha
   beta  = starting$beta
   delta = starting$delta
   gamma = starting$gamma
   A     = starting$A
   tau   = starting$tau
   rho   = starting$rho
   
   # set precisions with supplied values for the Metropolis sampler Normal proposal distributions 
   phi_alpha = tuning$phi_alpha
   phi_beta  = tuning$phi_beta
   phi_delta = tuning$phi_delta
   phi_gamma = tuning$phi_gamma
   phi_A     = tuning$phi_A
   
   # set hyperparameters with supplied values
   nu = priors$nu
   b  = priors$b
   c  = priors$c
   
   # initialize MCMC samples
   A_MCMC     = array(NA, dim = c(p, p, n_samples))
   alpha_MCMC = array(NA, dim = c(p, p, n_samples))
   beta_MCMC  = array(NA, dim = c(p, p, n_samples))
   delta_MCMC = matrix(NA, p, n_samples)
   gamma_MCMC = matrix(NA, p, n_samples)
   tau_MCMC   = matrix(NA, 4, n_samples)
   rho_MCMC   = rep(NA, n_samples)
   
   # initialize acceptance indicators
   accept_alpha = accept_beta = accept_A = array(0, dim = c(p, p, n_samples)) 
   accept_delta = accept_gamma = matrix(0, p, n_samples)
   
   # calculate logit(pi) and log(lambda)
   logitPi   = tcrossprod(x, alpha) + + matrix(delta, n, p, byrow = TRUE)
   logLambda = tcrossprod(x, beta) + matrix(gamma, n, p, byrow = TRUE)
   
   # do MCMC iterations
   for (t in 1 : n_samples)
   {
      # update alpha (Metropolis-Hastings step)
      updates = update_alpha(alpha, A, tau, x, logitPi, logLambda, phi_alpha, nu)
      alpha   = updates$alpha
      logitPi = updates$logitPi
      accept_alpha[ , , t] = updates$accept
      
      # update beta (Metropolis-Hastings step)
      updates   = update_beta(beta, A, tau, x, logitPi, logLambda, phi_beta, nu)
      beta      = updates$beta
      logLambda = updates$logLambda
      accept_beta[ , , t] = updates$accept
      
      # update delta (Metropolis-Hastings step)
      updates = update_delta(delta, tau, x, logitPi, logLambda, phi_delta, nu)
      delta   = updates$delta
      logitPi = updates$logitPi
      accept_delta[ , t] = updates$accept
      
      # update gamma (Metropolis-Hastings step)
      updates   = update_gamma(gamma, tau, x, logitPi, logLambda, phi_gamma, nu)
      gamma     = updates$gamma
      logLambda = updates$logLambda
      accept_gamma[ , t] = updates$accept
      
      # update A (Metropolis-Hastings step)
      updates   = update_A(A, alpha, beta, tau, rho, x, logitPi, logLambda, phi_A, nu)
      A         = updates$A
      alpha     = updates$alpha
      beta      = updates$beta
      logitPi   = updates$logitPi
      logLambda = updates$logLambda
      accept_A[ , , t] = updates$accept
      
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
      
      # print progress of the sampler
      if (t %% 100 == 0) 
         cat("iter=", t, "\n")
      if (t %% 1000 == 0)
      {
         cat("acceptance rates of alpha: \n")
         print(apply(accept_alpha[ , , (t - 999) : t], c(1, 2), mean))
         cat("acceptance rates of beta: \n")
         print(apply(accept_beta[ , , (t - 999) : t], c(1, 2), mean))
         cat("acceptance rates of delta: \n")
         print(apply(accept_delta[ , (t - 999) : t], 1, mean))
         cat("acceptance rates of gamma: \n")
         print(apply(accept_gamma[ , (t - 999) : t], 1, mean))
         cat("acceptance rates of A: \n")
         print(apply(accept_A[ , , (t - 999) : t], c(1, 2), mean))
      }
   }
   
   # return a list of MCMC samples and acceptance indicators
   return(list(alpha = alpha_MCMC,
               beta  = beta_MCMC,
               delta = delta_MCMC,
               gamma = gamma_MCMC,
               A     = A_MCMC,
               tau   = tau_MCMC,
               rho   = rho_MCMC,
               accept_alpha = accept_alpha,
               accept_beta  = accept_beta,
               accept_delta = accept_delta,
               accept_gamma = accept_gamma,
               accept_A = accept_A))
}

#' Evaluate log-likelihood of each observation for the j-th component of ZIPBN model
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

# update each element of alpha through Metropolis step
update_alpha = function(alpha, A, tau, x, logitPi, logLambda, phi_alpha, nu)
{
   p      = ncol(x)
   accept = matrix(0, p, p)
   
   for (j in 1 : p)
   {
      x_j = x[ , j]
      logitPi_j   = logitPi[ , j]
      logLambda_j = logLambda[ , j]
      
      for (k in 1 : p)
      {
         if (j != k)
         {
            alpha_old = alpha[j, k]
            
            if (A[j, k] == 0)
            {
               alpha_new   = rnorm(1, mean = alpha_old, sd = sqrt(1 / phi_alpha[1]))
               logitPi_new = logitPi_j + x[ , k] * (alpha_new - alpha_old)
               llik_old    = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
               llik_new    = llik_ZIPBN_j(x_j, logitPi_new, logLambda_j)
               ratio_MH    = exp(sum(llik_new - llik_old) - 0.5 * nu * tau[1] * (alpha_new * alpha_new - alpha_old * alpha_old))
            } else
            {
               alpha_new   = rnorm(1, mean = alpha_old, sd = sqrt(1 / phi_alpha[2]))
               logitPi_new = logitPi_j + x[ , k] * (alpha_new - alpha_old)
               llik_old    = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
               llik_new    = llik_ZIPBN_j(x_j, logitPi_new, logLambda_j)
               ratio_MH    = exp(sum(llik_new - llik_old) - 0.5 * tau[1] * (alpha_new * alpha_new - alpha_old * alpha_old))
            }
            
            if (runif(1) < min(1, ratio_MH))
            {
               alpha[j, k]  = alpha_new
               logitPi_j    = logitPi_new
               accept[j, k] = 1
            }
         }
      }
      
      logitPi[ , j] = logitPi_j
   }
   
   return(list(alpha   = alpha, 
               logitPi = logitPi,
               accept  = accept))
}

# update each element of beta through Metropolis step
update_beta = function(beta, A, tau, x, logitPi, logLambda, phi_beta, nu)
{
   p      = ncol(x)
   accept = matrix(0, p, p)
   
   for (j in 1 : p)
   {
      x_j = x[ , j]
      logitPi_j   = logitPi[ , j]
      logLambda_j = logLambda[ , j]
      
      for (k in 1 : p)
      {
         if (j != k)
         {
            beta_old = beta[j, k]
            
            if(A[j, k] == 0)
            {
               beta_new      = rnorm(1, mean = beta_old, sd = sqrt(1 / phi_beta[1]))
               logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old)
               llik_old      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
               llik_new      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_new)
               ratio_MH      = exp(sum(llik_new - llik_old) - 0.5 * nu * tau[2] * (beta_new * beta_new - beta_old * beta_old))
            } else
            {
               beta_new      = rnorm(1, mean = beta_old, sd = sqrt(1 / phi_beta[2]))
               logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old)
               llik_old      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
               llik_new      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_new)
               ratio_MH      = exp(sum(llik_new - llik_old) - 0.5 * tau[2] * (beta_new * beta_new - beta_old * beta_old)) 
            }
            
            if (runif(1) < min(1, ratio_MH))
            {
               beta[j, k]   = beta_new
               logLambda_j  = logLambda_new
               accept[j, k] = 1
            }
         }
      }
      
      logLambda[ , j] = logLambda_j
   }
   
   return(list(beta      = beta,
               logLambda = logLambda,
               accept    = accept))
}

# update each element of delta through Metropolis step
update_delta = function(delta, tau, x, logitPi, logLambda, phi_delta, nu)
{
   p      = ncol(x)
   accept = rep(0, p)
   
   for (j in 1 : p)
   {
      x_j         = x[ , j]
      logitPi_j   = logitPi[ , j]
      logLambda_j = logLambda[ , j]
      
      delta_old   = delta[j]
      
      delta_new   = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_delta))
      logitPi_new = logitPi_j + (delta_new - delta_old)
      llik_old    = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
      llik_new    = llik_ZIPBN_j(x_j, logitPi_new, logLambda_j)
      ratio_MH    = exp(sum(llik_new - llik_old) - 0.5 * tau[3] * (delta_new * delta_new - delta_old * delta_old))
      
      if (runif(1) < min(1, ratio_MH))
      {
         delta[j]      = delta_new
         logitPi[ , j] = logitPi_new
         accept[j]     = 1
      }
   }
   
   return(list(delta   = delta,
               logitPi = logitPi,
               accept  = accept))
}

# update each element of gamma through Metropolis step
update_gamma = function(gamma, tau, x, logitPi, logLambda, phi_gamma, nu)
{
   p      = ncol(x)
   accept = rep(0, p)
   
   for (j in 1 : p)
   {
      x_j         = x[ , j]
      logitPi_j   = logitPi[ , j]
      logLambda_j = logLambda[ , j]
      
      gamma_old = gamma[j]
      
      gamma_new     = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_gamma))
      logLambda_new = logLambda_j + (gamma_new - gamma_old)
      llik_old      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
      llik_new      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_new)
      ratio_MH      = exp(sum(llik_new - llik_old) - 0.5 * tau[4] * (gamma_new * gamma_new - gamma_old * gamma_old))
      
      if (runif(1) < min(1, ratio_MH))
      {
         gamma[j]        = gamma_new
         logLambda[ , j] = logLambda_new
         accept[j]       = 1
      }
   }
   
   return(list(gamma     = gamma,
               logLambda = logLambda,
               accept    = accept))
}

# update each element of A jointly with corresponding alpha and beta through Metropolis step
update_A = function(A, alpha, beta, tau, rho, x, logitPi, logLambda, phi_A, nu)
{
   p      = ncol(x)
   accept = matrix(0, p, p)
   
   for (j in 1 : p)
   {
      x_j = x[ , j]
      logitPi_j   = logitPi[ , j]
      logLambda_j = logLambda[ , j]
      
      for (k in 1 : p)
      {
         if (j != k)
         {
            alpha_old = alpha[j, k]
            beta_old  = beta[j, k]
            
            if (A[j, k] == 0)
            {
               A[j, k] = A_new = 1
               
               if (is_dag(graph_from_adjacency_matrix(A)))
               {
                  alpha_new = rnorm(1, mean = alpha_old, sd = sqrt(1 / phi_A[1]))
                  beta_new  = rnorm(1, mean = beta_old, sd = sqrt(1 / phi_A[2]))
                  
                  logitPi_new   = logitPi_j + x[ , k] * (alpha_new - alpha_old)
                  logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old)
                  llik_old  = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
                  llik_new  = llik_ZIPBN_j(x_j, logitPi_new, logLambda_new)
                  
                  ratio_MH = dnorm(alpha_old, mean = 0, sd = sqrt(1 / phi_A[3]), log = TRUE) + 
                     dnorm(beta_old, mean = 0, sd = sqrt(1 / phi_A[3]), log = TRUE) - 
                     dnorm(alpha_new, mean = alpha_old, sd = sqrt(1 / phi_A[1]), log = TRUE) - 
                     dnorm(beta_new, mean = beta_old, sd = sqrt(1 / phi_A[2]), log = TRUE)
                  ratio_MH = ratio_MH + sum(llik_new - llik_old) - log(nu) - 
                     0.5 * tau[1] * (alpha_new * alpha_new - nu * alpha_old * alpha_old) -
                     0.5 * tau[2] * (beta_new * beta_new - nu * beta_old * beta_old) + 
                     log(rho) - log(1 - rho)
                  ratio_MH = exp(ratio_MH)
               } else
               {
                  ratio_MH = 0
               }
               
               A[j, k] = 0
            } else
            {
               A_new     = 0
               alpha_new = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[3]))
               beta_new  = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[3]))
               
               logitPi_new   = logitPi_j + x[ , k] * (alpha_new - alpha_old)
               logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old)
               llik_old  = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
               llik_new  = llik_ZIPBN_j(x_j, logitPi_new, logLambda_new)
               
               ratio_MH = dnorm(alpha_old, mean = alpha_new, sd = sqrt(1 / phi_A[1]), log = TRUE) + 
                  dnorm(beta_old, mean = beta_new, sd = sqrt(1 / phi_A[2]), log = TRUE) - 
                  dnorm(alpha_new, mean = 0, sd = sqrt(1 / phi_A[3]), log = TRUE) - 
                  dnorm(beta_new, mean = 0, sd = sqrt(1 / phi_A[3]), log = TRUE)
               ratio_MH = ratio_MH + sum(llik_new - llik_old) + log(nu) - 
                  0.5 * tau[1] * (nu * alpha_new * alpha_new - alpha_old * alpha_old) -
                  0.5 * tau[2] * (nu * beta_new * beta_new - beta_old * beta_old) +
                  log(1 - rho) - log(rho)
               ratio_MH = exp(ratio_MH)
            }
            
            if (runif(1) < min(1, ratio_MH))
            {
               A[j, k]      = A_new
               alpha[j, k]  = alpha_new
               beta[j, k]   = beta_new
               logitPi_j    = logitPi_new
               logLambda_j  = logLambda_new
               accept[j, k] = 1
            }
         }
      }
      
      logitPi[ , j]   = logitPi_j
      logLambda[ , j] = logLambda_j
   }
   
   return(list(A         = A,
               alpha     = alpha,
               beta      = beta,
               logitPi   = logitPi,
               logLambda = logLambda,
               accept    = accept))
}