#' Implementation of MCMC sampler for ZIPBN models
#'
#' @param x a matrix containing data
#' @param starting a list of parameters' starting values for MCMC
#' @param tuning a list of precision values for Metropolis sampler Normal proposal distribution
#' @param priors a list of hyperparameter values for priors
#' @param n_sample the number of MCMC iterations
#'
#' @return MCMC samples from posterior distributions of ZIPBN models
#' @export
#'
#' @examples
mcmc_ZIPBN = function(x, starting, tuning, priors, n_sample = 5000, n_burnin = 3000, verbose = TRUE, n_report = 500)
{
   # check compatibility of x
   if (missing(x))
      stop("error: x should be specified")
   if (class(x) != "matrix")
      stop("error: x should be a matrix")
   if (any(x != as.integer(x) | x < 0))
      stop("error: Each element of x should be non-negative integer")
   
   # store the sample size and the number of nodes
   n = nrow(x)
   p = ncol(x)
   
   # check compatibility of starting value list and 
   # initialize parameters with supplied starting values 
   if (missing(starting))
      stop("error: starting value list for the parameters should be specified")
   if (!"alpha" %in% names(starting))
      stop("error: alpha should be specified in starting value list")
   else
   {
      alpha = starting$alpha
      if (nrow(alpha) != p | ncol(alpha) != p)
         stop(paste("error: alpha should be a", p, "x", p, "matrix", sep = " "))
      if (any(diag(alpha) != 0))
         stop("error: diagonal of alpha should be zeros")
   }
   if (!"beta" %in% names(starting))
      stop("error: beta should be specified in starting value list")
   else
   {
      beta = starting$beta
      if (nrow(beta) != p | ncol(beta) != p)
         stop(paste("error: beta should be a", p, "x", p, "matrix", sep = " "))
      if (any(diag(beta) != 0))
         stop("error: diagonal of beta should be zeros")
   }
   if (!"delta" %in% names(starting))
      stop("error: delta should be specified in starting value list")
   else
   {
      delta = starting$delta
      if (length(delta) != p)
         stop(paste("error: delta should be a vector of length", p, sep = ""))
   }
   if (!"gamma" %in% names(starting))
      stop("error: gamma should be specified in starting value list")
   else
   {
      gamma = starting$gamma
      if (length(gamma) != p)
         stop(paste("error: gamma should be a vector of length", p, sep = " "))
   }
   if (!"A" %in% names(starting))
      stop("error: A should be specified in starting value list")
   else
   {
      A = starting$A
      if (nrow(A) != p | ncol(A) != p | any(A != 0 & A != 1))
         stop(paste("error: A should be a", p, "x", p, "matrix, of which each element is 0 or 1", sep = " "))
      if (any(diag(A) != 0))
         stop("error: diagonal of A should be zeros")
   }
   if (!"tau" %in% names(starting))
      stop("error: tau should be specified in starting value list")
   else
   {
      tau = starting$tau
      if (length(tau) != 4 | any(tau <= 0))
         stop("error: tau should be a positive vector of length 4")
   }
   if (!"rho" %in% names(starting))
      stop("error: rho should be specified in starting value list")
   else
   {
      rho = starting$rho
      if (rho < 0 | rho > 1)
         stop("error: rho should be between 0 and 1")
   }
   
   
   # check compatibility of tuning value list and 
   # set precisions of Normal proposal distributions for Metropolis samplers with supplied tuning values 
   if (missing(tuning))
      stop("error: tuning value list for Metropolis sampler should be specified")
   if (!"phi_alpha" %in% names(tuning))
      stop("error: phi_alpha should be specified in tuning value list")
   else
   {
      phi_alpha = tuning$phi_alpha
      if (length(phi_alpha) != 2 | any(phi_alpha <= 0) | phi_alpha[1] < phi_alpha[2])
         stop("error: phi_alpha should be a positive vector of length 2, of which the first element is greater than the second element")
   }
   if (!"phi_beta" %in% names(tuning))
      stop("error: phi_beta should be specified in tuning value list")
   else
   {
      phi_beta = tuning$phi_beta
      if (length(phi_beta) != 2 | any(phi_beta <= 0) | phi_beta[1] <= phi_beta[2])
         stop("error: phi_beta should be a positive vector of length 2, of which the first element is greater than the second element")
   }
   if (!"phi_delta" %in% names(tuning))
      stop("error: phi_delta should be specified in tuning value list")
   else
   {
      phi_delta = tuning$phi_delta
      if (phi_delta <= 0)
         stop("error: phi_delta should be non-negative")
   }
   if (!"phi_gamma" %in% names(tuning))
      stop("error: phi_gamma should be specified in tuning value list")
   else
   {
      phi_gamma = tuning$phi_gamma
      if (phi_gamma <= 0)
         stop("error: phi_gamma should be non-negative")
   }
   if (!"phi_A" %in% names(tuning))
      stop("error: phi_A should be specified in tuning value list")
   else
   {
      phi_A = tuning$phi_A
      if (length(phi_A) != 5 | any(phi_A <= 0) | phi_A[1] <= phi_A[2] | phi_A[1] <= phi_A[3])
         stop("error: phi_A should be a positive vector of length 5, of which the first element is greater than the second and third elements")
   }

   # check compatibility of priors list and set hyperparameters with supplied prior list
   if (missing(priors)) 
      stop("error: prior list for the parameters should be specified")
   if (!"nu" %in% names(priors))
      stop("error: nu should be specified in prior list")
   if (priors$nu <= 0)
      stop("error: nu should be a positive value")
   if (!"tau_alpha" %in% names(priors))
      stop("error: tau_alpha should be specified in prior list")
   if (length(priors$tau_alpha) != 2 | any(priors$tau_alpha <= 0))
      stop("error: tau_alpha should be a positive vector of length 2")
   if (!"tau_beta" %in% names(priors))
      stop("error: tau_beta should be specified in prior list")
   if (length(priors$tau_beta) != 2 | any(priors$tau_beta <= 0))
      stop("error: tau_beta should be a positive vector of length 2")
   if (!"tau_delta" %in% names(priors))
      stop("error: tau_delta should be specified in prior list")
   if (length(priors$tau_delta) != 2 | any(priors$tau_delta <= 0))
      stop("error: tau_delta should be a positive vector of length 2")
   if (!"tau_gamma" %in% names(priors))
      stop("error: tau_gamma should be specified in prior list")
   if (length(priors$tau_gamma) != 2 | any(priors$tau_gamma <= 0))
      stop("error: tau_gamma should be a positive vector of length 2")
   if (!"rho" %in% names(priors))
      stop("error: rho should be specified in prior list")
   if (length(priors$rho) != 2 | any(priors$rho <= 0))
      stop("error: rho should be a positive vector of length 2")
   nu = priors$nu
   b  = c(priors$tau_alpha[1], priors$tau_beta[1], priors$tau_delta[1], priors$tau_gamma[1], priors$rho[1])
   c  = c(priors$tau_alpha[2], priors$tau_beta[2], priors$tau_delta[2], priors$tau_gamma[2], priors$rho[2])

   # check compatibility of n_sample, n_burnin, verbose, and n_report
   if (n_sample != as.integer(n_sample) | n_sample <= 0)
      stop("error: n_sample should be a natural number")
   if (n_burnin != as.integer(n_burnin) | n_burnin <= 0 | n_burnin >= n_sample)
      stop("error: n_burnin should be a natural number less than n_sample")
   if (class(verbose) != "logical")
      stop("error: verbose should be a logical value")
   if (n_report != as.integer(n_report) | n_report <= 0 | n_report >= n_sample)
      stop("error: n_report should be a natural number less than n_sample")
   
   # initialize MCMC samples
   A_MCMC     = array(NA, dim = c(p, p, n_sample))
   alpha_MCMC = array(NA, dim = c(p, p, n_sample))
   beta_MCMC  = array(NA, dim = c(p, p, n_sample))
   delta_MCMC = matrix(NA, p, n_sample)
   gamma_MCMC = matrix(NA, p, n_sample)
   tau_MCMC   = matrix(NA, 4, n_sample)
   rho_MCMC   = rep(NA, n_sample)
   
   # initialize acceptance indicators
   accept_alpha = accept_beta = array(0, dim = c(p, p, n_sample)) 
   accept_delta = accept_gamma = matrix(0, p, n_sample)
   
   # calculate logit(pi) and log(lambda) with the starting values
   logitPi   = tcrossprod(x, alpha) + + matrix(delta, n, p, byrow = TRUE)
   logLambda = tcrossprod(x, beta) + matrix(gamma, n, p, byrow = TRUE)
   
   
   ###--------------------------- Metropolis-within-Gibbs sampler for ZIPBN models -------------------------###
   # sample from the ZIPBN posterior through MCMC iterations
   for (t in 1 : n_sample)
   {
      # sample alpha given the other parameters (Metropolis step)
      samp_alpha = MH_alpha(alpha, A, tau, x, logitPi, logLambda, phi_alpha, nu)
      alpha      = samp_alpha$alpha
      logitPi    = samp_alpha$logitPi
      accept_alpha[ , , t] = samp_alpha$accept
      
      # sample beta given the other parameters (Metropolis step)
      samp_beta = MH_beta(beta, A, tau, x, logitPi, logLambda, phi_beta, nu)
      beta      = samp_beta$beta
      logLambda = samp_beta$logLambda
      accept_beta[ , , t] = samp_beta$accept
      
      # sample delta given the other parameters (Metropolis step)
      samp_delta = MH_delta(delta, tau, x, logitPi, logLambda, phi_delta, nu)
      delta      = samp_delta$delta
      logitPi    = samp_delta$logitPi
      accept_delta[ , t] = samp_delta$accept
      
      # sample gamma given the other parameters (Metropolis step)
      samp_gamma = MH_gamma(gamma, tau, x, logitPi, logLambda, phi_gamma, nu)
      gamma      = samp_gamma$gamma
      logLambda  = samp_gamma$logLambda
      accept_gamma[ , t] = samp_gamma$accept
      
      # sample A jointly with alpha, beta, delta, and gamma, given tau and rho (Metropolis step)
      # toss a coin to determine the strategy of sampling A 
      # strategy 1: Metropolis sampler proposing addition or deletion of edges 
      # strategy 2: Metropolis sampler proposing reversal of edges   
      if (runif(1) < 0.5)
         samp_A = MH_A_each(A, alpha, beta, delta, gamma, tau, rho, x, logitPi, logLambda, phi_A, nu)
      else
         samp_A = MH_A_rev(A, alpha, beta, delta, gamma, tau, x, logitPi, logLambda, phi_A, nu)
      A         = samp_A$A
      alpha     = samp_A$alpha
      beta      = samp_A$beta
      delta     = samp_A$delta
      gamma     = samp_A$gamma
      logitPi   = samp_A$logitPi
      logLambda = samp_A$logLambda
      
      # sample taus from their full conditionals
      tau[1] = rgamma(1, shape = b[1] + p * (p - 1) / 2, 
                      rate = c[1] + (nu * sum((A == 0) * alpha * alpha) + sum((A == 1) * alpha * alpha)) / 2)   # tau_alpha
      tau[2] = rgamma(1, shape = b[2] + p * (p - 1) / 2, 
                      rate = c[2] + (nu * sum((A == 0) * beta * beta) + sum((A == 1) * beta * beta)) / 2)   # tau_beta
      tau[3] = rgamma(1, shape = b[3] + p / 2, rate = c[3] + sum(delta * delta) / 2)   # tau_delta
      tau[4] = rgamma(1, shape = b[4] + p / 2, rate = c[4] + sum(gamma * gamma) / 2)   # tau_gamma
      
      # sample rho from its full conditional
      rho = rbeta(1, shape1 = b[5] + sum(A == 1), shape2 = c[5] + sum(A == 0) - p)
      
      # store MCMC samples of iteration t
      alpha_MCMC[ , , t] = alpha
      beta_MCMC[ , , t]  = beta
      delta_MCMC[ , t]   = delta
      gamma_MCMC[ , t]   = gamma
      A_MCMC[ , , t]     = A
      tau_MCMC[ , t]     = tau
      rho_MCMC[t]        = rho
      
      # print progress of the sampler and Metropolis sampler acceptance
      if (verbose)
      {
         if (t %% n_report == 0)
         {
            cat("\n")
            cat("iter=", t, "\n")
            cat("acceptance rates of alpha: \n")
            print(round(100 * apply(accept_alpha[ , , 1 : t], c(1, 2), mean), digits = 2))
            cat("acceptance rates of beta: \n")
            print(round(100 * apply(accept_beta[ , , 1 : t], c(1, 2), mean), digits = 2))
            cat("acceptance rates of delta: \n")
            print(round(100 * apply(accept_delta[ , 1 : t], 1, mean), digits = 2))
            cat("acceptance rates of gamma: \n")
            print(round(100 * apply(accept_delta[ , 1 : t], 1, mean), digits = 2))
            cat("\n")
         }
      }
   }
   
   # return a list of 1. MCMC samples (after burn-in period) for each parameter 
   #                  2. Metropolis sampler acceptance rates for alpha, beta, delta, and gamma
   results = list()
   results$alpha = alpha_MCMC[ , , (n_burnin + 1) : n_sample]
   results$beta  = beta_MCMC[ , , (n_burnin + 1) : n_sample]
   results$delta = delta_MCMC[ , (n_burnin + 1) : n_sample]
   results$gamma = gamma_MCMC[ , (n_burnin + 1) : n_sample]
   results$A     = A_MCMC[ , , (n_burnin + 1) : n_sample]
   results$tau   = tau_MCMC[ , (n_burnin + 1) : n_sample]
   results$rho   = rho_MCMC[(n_burnin + 1) : n_sample]
   results$acceptance = list(alpha = 100 * apply(accept_alpha, c(1, 2), mean),
                             beta  = 100 * apply(accept_beta, c(1, 2), mean),
                             gamma = 100 * apply(accept_delta, 1, mean),
                             delta = 100 * apply(accept_delta, 1, mean))
   class(results) = "ZIPBN"
   return(results)
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
MH_alpha = function(alpha, A, tau, x, logitPi, logLambda, phi_alpha, nu)
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

# sample each element of beta through Metropolis step
MH_beta = function(beta, A, tau, x, logitPi, logLambda, phi_beta, nu)
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

# sample each element of delta through Metropolis step
MH_delta = function(delta, tau, x, logitPi, logLambda, phi_delta, nu)
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

# sample each element of gamma through Metropolis step
MH_gamma = function(gamma, tau, x, logitPi, logLambda, phi_gamma, nu)
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

# sample each element of A through Metropolis step (propose addition or deletion of an edge)
# corresponding alpha, beta, delta, and gamma are jointly proposed with A 
MH_A_each = function(A, alpha, beta, delta, gamma, tau, rho, x, logitPi, logLambda, phi_A, nu)
{
   p = ncol(x)
   
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
            delta_old = delta[j]
            gamma_old = gamma[j]
            
            if (A[j, k] == 0)
            {
               A[j, k] = A_new = 1
               graph   = graph_from_adjacency_matrix(A)
               A[j, k] = 0
               
               if (is_dag(graph))
               {
                  alpha_new = rnorm(1, mean = alpha_old, sd = sqrt(1 / phi_A[2]))
                  beta_new  = rnorm(1, mean = beta_old, sd = sqrt(1 / phi_A[3]))
                  delta_new = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_A[4]))
                  gamma_new = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_A[5]))
                  
                  logitPi_new   = logitPi_j + x[ , k] * (alpha_new - alpha_old) + (delta_new - delta_old)
                  logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old) + (gamma_new - gamma_old)
                  llik_old  = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
                  llik_new  = llik_ZIPBN_j(x_j, logitPi_new, logLambda_new)
                  
                  ratio_MH = dnorm(alpha_old, mean = 0, sd = sqrt(1 / phi_A[1]), log = TRUE) + 
                     dnorm(beta_old, mean = 0, sd = sqrt(1 / phi_A[1]), log = TRUE) - 
                     dnorm(alpha_new, mean = alpha_old, sd = sqrt(1 / phi_A[2]), log = TRUE) - 
                     dnorm(beta_new, mean = beta_old, sd = sqrt(1 / phi_A[3]), log = TRUE)
                  ratio_MH = ratio_MH + sum(llik_new - llik_old) - log(nu) - 
                     0.5 * tau[1] * (alpha_new * alpha_new - nu * alpha_old * alpha_old) -
                     0.5 * tau[2] * (beta_new * beta_new - nu * beta_old * beta_old) - 
                     0.5 * tau[3] * (delta_new * delta_new - delta_old * delta_old) - 
                     0.5 * tau[4] * (gamma_new * gamma_new - gamma_old * gamma_old) +
                     log(rho) - log(1 - rho)
                  ratio_MH = exp(ratio_MH)
                  
                  if (runif(1) < min(1, ratio_MH))
                  {
                     cat("An edge", k, "->", j, "is added \n")
                     A[j, k]      = A_new
                     alpha[j, k]  = alpha_new
                     beta[j, k]   = beta_new
                     delta[j]     = delta_new
                     gamma[j]     = gamma_new
                     logitPi_j    = logitPi_new
                     logLambda_j  = logLambda_new
                  }
               } 
            } else
            {
               A_new = 0
               alpha_new = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
               beta_new  = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
               delta_new = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_A[4]))
               gamma_new = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_A[5]))
               
               logitPi_new   = logitPi_j + x[ , k] * (alpha_new - alpha_old) + (delta_new - delta_old)
               logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old) + (gamma_new - gamma_old)
               llik_old  = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
               llik_new  = llik_ZIPBN_j(x_j, logitPi_new, logLambda_new)
               
               ratio_MH = dnorm(alpha_old, mean = alpha_new, sd = sqrt(1 / phi_A[2]), log = TRUE) + 
                  dnorm(beta_old, mean = beta_new, sd = sqrt(1 / phi_A[3]), log = TRUE) - 
                  dnorm(alpha_new, mean = 0, sd = sqrt(1 / phi_A[1]), log = TRUE) - 
                  dnorm(beta_new, mean = 0, sd = sqrt(1 / phi_A[1]), log = TRUE)
               ratio_MH = ratio_MH + sum(llik_new - llik_old) + log(nu) - 
                  0.5 * tau[1] * (nu * alpha_new * alpha_new - alpha_old * alpha_old) -
                  0.5 * tau[2] * (nu * beta_new * beta_new - beta_old * beta_old) - 
                  0.5 * tau[3] * (delta_new * delta_new - delta_old * delta_old) - 
                  0.5 * tau[4] * (gamma_new * gamma_new - gamma_old * gamma_old) +
                  log(1 - rho) - log(rho)
               ratio_MH = exp(ratio_MH)
               
               if (runif(1) < min(1, ratio_MH))
               {
                  cat("An edge", k, "->", j, "is deleted \n")
                  A[j, k]      = A_new
                  alpha[j, k]  = alpha_new
                  beta[j, k]   = beta_new
                  delta[j]     = delta_new
                  gamma[j]     = gamma_new
                  logitPi_j    = logitPi_new
                  logLambda_j  = logLambda_new
               }
            }
         }
      }
      
      logitPi[ , j]   = logitPi_j
      logLambda[ , j] = logLambda_j
   }
   
   return(list(A         = A,
               alpha     = alpha,
               beta      = beta,
               delta     = delta,
               gamma     = gamma,
               logitPi   = logitPi,
               logLambda = logLambda))
}

# sample A based on proposal of reversing an edge through Metropolis step 
# corresponding alpha, beta, delta, and gamma are jointly proposed with A 
MH_A_rev = function(A, alpha, beta, delta, gamma, tau, x, logitPi, logLambda, phi_A, nu)
{
   ids_A1 = which(A == 1, arr.ind = TRUE)
   n_swap = nrow(ids_A1)
   
   if (n_swap > 0)
   {
      
      for (s in 1 : n_swap)
      {
         j1 = k0 = ids_A1[s, 1]
         j0 = k1 = ids_A1[s, 2]
         
         A[j1, k1] = 0
         A[j0, k0] = 1
         graph    = graph_from_adjacency_matrix(A)
         A[j1, k1] = 1
         A[j0, k0] = 0
         
         if (is_dag(graph))
         {
            x_j1 = x_k0 = x[ , j1]
            x_j0 = x_k1 = x[ , j0]
            logitPi_j1   = logitPi[ , j1]
            logitPi_j0   = logitPi[ , j0]
            logLambda_j1 = logLambda[ , j1]
            logLambda_j0 = logLambda[ , j0]
            
            alpha_old1 = alpha[j1, k1]
            alpha_old0 = alpha[j0, k0]
            beta_old1  = beta[j1, k1]
            beta_old0  = beta[j0, k0]
            delta_old1 = delta[j1]
            delta_old0 = delta[j0]
            gamma_old1 = gamma[j1]
            gamma_old0 = gamma[j0]
            
            alpha_new1 = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
            alpha_new0 = rnorm(1, mean = alpha_old0, sd = sqrt(1 / phi_A[2]))
            beta_new1  = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
            beta_new0  = rnorm(1, mean = beta_old0, sd = sqrt(1 / phi_A[3]))
            delta_new1 = rnorm(1, mean = delta_old1, sd = sqrt(1 / phi_A[4]))
            delta_new0 = rnorm(1, mean = delta_old0, sd = sqrt(1 / phi_A[4]))
            gamma_new1 = rnorm(1, mean = gamma_old1, sd = sqrt(1 / phi_A[5]))
            gamma_new0 = rnorm(1, mean = gamma_old0, sd = sqrt(1 / phi_A[5]))
            
            logitPi_new1   = logitPi_j1 + x_k1 * (alpha_new1 - alpha_old1) + (delta_new1 - delta_old1)
            logitPi_new0   = logitPi_j0 + x_k0 * (alpha_new0 - alpha_old0) + (delta_new0 - delta_old0)
            logLambda_new1 = logLambda_j1 + x_k1 * (beta_new1 - beta_old1) + (gamma_new1 - gamma_old1)
            logLambda_new0 = logLambda_j0 + x_k0 * (beta_new0 - beta_old0) + (gamma_new0 - gamma_old0)
            
            llik_old1 = llik_ZIPBN_j(x_j1, logitPi_j1, logLambda_j1)
            llik_old0 = llik_ZIPBN_j(x_j0, logitPi_j0, logLambda_j0)
            llik_new1 = llik_ZIPBN_j(x_j1, logitPi_new1, logLambda_new1)
            llik_new0 = llik_ZIPBN_j(x_j0, logitPi_new0, logLambda_new0)
            
            ratio_MH = dnorm(alpha_old0, mean = 0, sd = sqrt(1 / phi_A[1]), log = TRUE) + 
               dnorm(beta_old0, mean = 0, sd = sqrt(1 / phi_A[1]), log = TRUE) + 
               dnorm(alpha_old1, mean = alpha_new1, sd = sqrt(1 / phi_A[2]), log = TRUE) + 
               dnorm(beta_old1, mean = beta_new1, sd = sqrt(1 / phi_A[3]), log = TRUE) - 
               dnorm(alpha_new0, mean = alpha_old0, sd = sqrt(1 / phi_A[2]), log = TRUE) - 
               dnorm(beta_new0, mean = beta_old0, sd = sqrt(1 / phi_A[3]), log = TRUE) - 
               dnorm(alpha_new1, mean = 0, sd = sqrt(1 / phi_A[1]), log = TRUE) - 
               dnorm(beta_new1, mean = 0, sd = sqrt(1 / phi_A[1]), log = TRUE)
            ratio_MH = ratio_MH + sum(llik_new0 - llik_old0) + sum(llik_new1 - llik_old1) - 
               0.5 * (tau[1] * (alpha_new0 * alpha_new0 - alpha_old1 * alpha_old1) + 
                         tau[2] * (beta_new0 * beta_new0 - beta_old1 * beta_old1)) - 
               0.5 * nu * (tau[1] * (alpha_new1 * alpha_new1 - alpha_old0 * alpha_old0) + 
                              tau[2] * (beta_new1 * beta_new1 - beta_old0 * beta_old0)) - 
               0.5 * tau[3] * (delta_new0 * delta_new0 - delta_old0 * delta_old0 + 
                                  delta_new1 * delta_new1 - delta_old1 * delta_old1) - 
               0.5 * tau[4] * (gamma_new0 * gamma_new0 - gamma_old0 * gamma_old0 +
                                  gamma_new1 * gamma_new1 - gamma_old1 * gamma_old1)
            ratio_MH = exp(ratio_MH)
            
            if (runif(1) < min(1, ratio_MH))
            {
               cat("An edge", k1, "->", j1, "is reversed to", k0, "->", j0, " \n")
               A[j1, k1] = 0
               A[j0, k0] = 1
               alpha[j1, k1] = alpha_new1
               alpha[j0, k0] = alpha_new0
               beta[j1, k1]  = beta_new1
               beta[j0, k0]  = beta_new0
               delta[j1]     = delta_new1
               delta[j0]     = delta_new0
               gamma[j1]     = gamma_new1
               gamma[j0]     = gamma_new0
               logitPi[ , j1]   = logitPi_new1
               logitPi[ , j0]   = logitPi_new0
               logLambda[ , j1] = logLambda_new1
               logLambda[ , j0] = logLambda_new0
            }
         } 
      } 
   }

   return(list(A         = A,
               alpha     = alpha,
               beta      = beta,
               delta     = delta,
               gamma     = gamma,
               logitPi   = logitPi,
               logLambda = logLambda))
}