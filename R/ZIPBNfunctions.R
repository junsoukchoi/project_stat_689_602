#' Implementation of MCMC sampler for ZIPBN models
#'
#' @param x a matrix containing data.
#' @param starting a list with each tag corresponding to a parameter name. 
#' Valid tags are 'alpha', 'beta', 'delta', 'gamma', 'A', 'tau', and 'rho'. 
#' The value portion of each tag is the parameters' starting values for MCMC.
#' @param tuning a list with each tag corresponding to a parameter name.
#' Valid tags are 'phi_alpha', 'phi_beta', 'phi_delta', 'phi_gamma', and 'phi_A'.
#' The value portion of each tag defines the precision of Normal proposal distribution for the Metropolis sampler. 
#' @param priors a list with each tag corresponding to a parameter name. 
#' Valid tags are 'nu', 'tau_alpha', 'tau_beta', 'tau_delta', 'tau_gamma', and 'rho'.
#' The value portion of each tag defines hyperparameters of the priors specified for ZIPBN models.
#' @param n_sample the number of MCMC iterations.
#' @param n_burnin the number of burn-in samples.
#' @param verbose if TRUE, progress of the sampler is printed to the screen. 
#' Otherwise, nothing is printed to the screen.
#' @param n_report the interval to report Metropolis sampler acceptance rates and MCMC progress.
#'
#' @return An object of class ZIPBN, which is a list with the tags 'samples' and 'acceptance'. 
#' The value portion of the tag 'samples' gives MCMC samples from posterior distributions for the defined parameters of ZIPBN models.
#' The value portion of the tag 'acceptance' shows the Metropolis sampling acceptance percents for alpha, beta, delta, and gamma. 
#' 
#' @export
#'
#' @examples
#'library(ZIPBN)
#'
#'
#'## Example data
#'set.seed(7)
#'
#'# generate a simple graph: X1 -> X2 -> X3
#'p = 3
#'A = matrix(0, p, p)
#'A[3, 2] = A[2, 1] = 1
#'
#'# parameters of the ZIPBN model, given graph A
#'alpha = matrix(0, p, p)
#'alpha[A == 1] = 0.3
#'beta  = matrix(0, p, p)
#'beta[A == 1] = 0.2
#'delta = rep(1, p)
#'gamma = rep(1.5, p)
#'
#'# generate data from the ZIPBN model
#'n = 200
#'x = matrix(0, n, p)
#'for (j in 1 : p)
#'{
#'   # calculate pi_j
#'   pi = exp(x %*% alpha[j, ] + delta[j])
#'   pi = pi / (1 + pi)
#'   # calculate mu_j
#'   mu = exp(x %*% beta[j, ] + gamma[j])
#'   # generate data for X_j
#'   x[ , j] = rpois(n, mu) * (1 - rbinom(n, 1, pi))
#'}
#'
#'
#'## fit ZIPBN models
#'# create starting value list
#'m = colMeans(x)
#'v = apply(x, 2, var)
#'starting = list(alpha = matrix(0, p, p),
#'                beta  = matrix(0, p, p),
#'                delta = log((v - m) / (m * m)),
#'                gamma = log((v - m + m * m) / m),
#'                A     = matrix(0, p, p),
#'                tau   = c(10, 10, 1, 1),
#'                rho   = 0.1)
#'
#'# create tuning value list
#'tuning = list(phi_alpha = c(1e+8, 20),
#'              phi_beta  = c(1e+8, 100),
#'              phi_delta = 5,
#'              phi_gamma = 50,
#'              phi_A     = c(1e+10, 10, 10, 1, 10))
#'
#'# create priors list
#'priors = list(nu        = 10000^2,
#'              tau_alpha = c(0.01, 0.01),
#'              tau_beta  = c(0.01, 0.01),
#'              tau_delta = c(0.01, 0.01),
#'              tau_gamma = c(0.01, 0.01),
#'              rho       = c(0.5, 0.5))
#'
#'# run mcmc_ZIPBN function
#'n_sample = 2000
#'n_burnin = 1000
#'out = mcmc_ZIPBN(x, starting, tuning, priors, n_sample, n_burnin)
#'
#'
#'## posterior inference via ZIPBN models
#'# report Metropolis sampling acceptance percents
#'out$acceptance
#'
#'# recover garph structure
#'cutoff = 0.5
#'A_est  = 1 * (apply(out$samples$A, c(1, 2), mean) > cutoff)
#'
#'# calculate the posterior mean of each parameter, given the recovered graph 
#'subset = apply(out$samples$A == array(A_est, dim = c(p, p, n_sample - n_burnin)), 3, all)
#'alpha_est = apply(out$samples$alpha[ , , subset], c(1, 2), mean)
#'beta_est  = apply(out$samples$beta[ , , subset], c(1, 2), mean)
#'delta_est = rowMeans(out$samples$delta[ , subset])
#'gamma_est = rowMeans(out$samples$gamma[ , subset])
#'
#'# report the posterior mean of each parameter with the recoverd graph
#'A_est
#'round(alpha_est, digits = 2)
#'round(beta_est, digits = 2)
#'round(delta_est, digits = 2)
#'round(gamma_est, digits = 2)
mcmc_ZIPBN = function(x, starting, tuning, priors, n_sample = 5000, n_burnin = 2500, verbose = TRUE, n_report = 500)
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
   
   
   ## Metropolis-within-Gibbs sampler for ZIPBN models
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
         samp_A = MH_A_each(A, alpha, beta, delta, gamma, tau, rho, x, logitPi, logLambda, phi_A, nu, verbose)
      else
         samp_A = MH_A_rev(A, alpha, beta, delta, gamma, tau, x, logitPi, logLambda, phi_A, nu, verbose)
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
      
      # save MCMC samples of iteration t
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

   
   # return a list of 
   # 1. MCMC samples (after burn-in period) for each parameter 
   # 2. Metropolis sampler acceptance rates for alpha, beta, delta, and gamma
   results = list()
   results$samples = list(alpha = alpha_MCMC[ , , (n_burnin + 1) : n_sample],
                          beta  = beta_MCMC[ , , (n_burnin + 1) : n_sample],
                          delta = delta_MCMC[ , (n_burnin + 1) : n_sample],
                          gamma = gamma_MCMC[ , (n_burnin + 1) : n_sample],
                          A     = A_MCMC[ , , (n_burnin + 1) : n_sample],
                          tau   = tau_MCMC[ , (n_burnin + 1) : n_sample],
                          rho   = rho_MCMC[(n_burnin + 1) : n_sample])
   results$acceptance = list(alpha = 100 * apply(accept_alpha, c(1, 2), mean),
                             beta  = 100 * apply(accept_beta, c(1, 2), mean),
                             gamma = 100 * apply(accept_delta, 1, mean),
                             delta = 100 * apply(accept_delta, 1, mean))
   class(results) = "ZIPBN"
   return(results)
}


# evaluate log-likelihood of each observation for the j-th component of ZIPBN model
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


# sample each element of alpha through Metropolis step
MH_alpha = function(alpha, A, tau, x, logitPi, logLambda, phi_alpha, nu)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = matrix(0, p, p)
   
   for (j in 1 : p)
   {
      # get data, logit(pi) and log(lambda) for the node j
      x_j = x[ , j]
      logitPi_j   = logitPi[ , j]
      logLambda_j = logLambda[ , j]
      
      for (k in 1 : p)
      {
         # if j = k, do not sample 
         if (j == k) next
         
         # current value of \alpha_{jk}
         alpha_old = alpha[j, k]
         
         # divide into two cases with and without the edge k -> j
         if (A[j, k] == 0)
         {
            # propose new value of \alpha_{jk} using Normal proposal distribution
            alpha_new   = rnorm(1, mean = alpha_old, sd = sqrt(1 / phi_alpha[1]))
            
            # calculate MH ratio
            logitPi_new = logitPi_j + x[ , k] * (alpha_new - alpha_old)
            llik_old    = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
            llik_new    = llik_ZIPBN_j(x_j, logitPi_new, logLambda_j)
            ratio_MH    = exp(sum(llik_new - llik_old) - 0.5 * nu * tau[1] * (alpha_new * alpha_new - alpha_old * alpha_old))
         } 
         else
         {
            # propose new value of \alpha_{jk} using Normal proposal distribution
            alpha_new   = rnorm(1, mean = alpha_old, sd = sqrt(1 / phi_alpha[2]))
            
            # calculate MH ratio
            logitPi_new = logitPi_j + x[ , k] * (alpha_new - alpha_old)
            llik_old    = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
            llik_new    = llik_ZIPBN_j(x_j, logitPi_new, logLambda_j)
            ratio_MH    = exp(sum(llik_new - llik_old) - 0.5 * tau[1] * (alpha_new * alpha_new - alpha_old * alpha_old))
         }
         
         # accept the proposed value with probability of min(1, MH ratio)
         if (runif(1) < min(1, ratio_MH))
         {
            alpha[j, k]  = alpha_new
            logitPi_j    = logitPi_new   # update logit(pi) for the node j 
            accept[j, k] = 1             # 1 if proposal accepted
         }
      }
      
      # save the updated logit(pi) for the node j
      logitPi[ , j] = logitPi_j
   }
   
   # return the updated alpha and logit(pi), and the acceptance indicators 
   return(list(alpha   = alpha, 
               logitPi = logitPi,
               accept  = accept))
}


# sample each element of beta through Metropolis step
MH_beta = function(beta, A, tau, x, logitPi, logLambda, phi_beta, nu)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = matrix(0, p, p)
   
   for (j in 1 : p)
   {
      # get data, logit(pi) and log(lambda) for the node j
      x_j = x[ , j]
      logitPi_j   = logitPi[ , j]
      logLambda_j = logLambda[ , j]
      
      for (k in 1 : p)
      {
         # if j = k, do not sample 
         if (j == k) next
         
         # current value of \beta_{jk}
         beta_old = beta[j, k]
         
         # divide into two cases with and without the edge k -> j
         if(A[j, k] == 0)
         {
            # propose new value of \beta_{jk} using Normal proposal distribution
            beta_new      = rnorm(1, mean = beta_old, sd = sqrt(1 / phi_beta[1]))
            
            # calculate MH ratio
            logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old)
            llik_old      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
            llik_new      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_new)
            ratio_MH      = exp(sum(llik_new - llik_old) - 0.5 * nu * tau[2] * (beta_new * beta_new - beta_old * beta_old))
         } 
         else
         {
            # propose new value of \beta_{jk} using Normal proposal distribution
            beta_new      = rnorm(1, mean = beta_old, sd = sqrt(1 / phi_beta[2]))
            
            # calculate MH ratio
            logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old)
            llik_old      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
            llik_new      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_new)
            ratio_MH      = exp(sum(llik_new - llik_old) - 0.5 * tau[2] * (beta_new * beta_new - beta_old * beta_old)) 
         }
         
         # accept the proposed value with probability of min(1, MH ratio)
         if (runif(1) < min(1, ratio_MH))
         {
            beta[j, k]   = beta_new
            logLambda_j  = logLambda_new   # update log(lambda) for the node j
            accept[j, k] = 1               # 1 if proposal accepted
         }
      }
      
      # save the updated log(lambda) for the node j
      logLambda[ , j] = logLambda_j
   }
   
   # return the updated beta and log(lambda), and the acceptance indicators 
   return(list(beta      = beta,
               logLambda = logLambda,
               accept    = accept))
}


# sample each element of delta through Metropolis step
MH_delta = function(delta, tau, x, logitPi, logLambda, phi_delta, nu)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = rep(0, p)
   
   for (j in 1 : p)
   {
      # get data, logit(pi) and log(lambda) for the node j
      x_j         = x[ , j]
      logitPi_j   = logitPi[ , j]
      logLambda_j = logLambda[ , j]
      
      # current value of \delta_j
      delta_old   = delta[j]
      
      # propose new value of \delta_j using Normal proposal distribution
      delta_new   = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_delta))
      
      # calculate MH ratio
      logitPi_new = logitPi_j + (delta_new - delta_old)
      llik_old    = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
      llik_new    = llik_ZIPBN_j(x_j, logitPi_new, logLambda_j)
      ratio_MH    = exp(sum(llik_new - llik_old) - 0.5 * tau[3] * (delta_new * delta_new - delta_old * delta_old))
      
      # accept the proposed value with probability of min(1, MH ratio)
      if (runif(1) < min(1, ratio_MH))
      {
         delta[j]      = delta_new
         logitPi[ , j] = logitPi_new   # update and save logit(pi) for the node j 
         accept[j]     = 1             # 1 if proposal accepted
      }
   }
   
   # return the updated delta and logit(pi), and the acceptance indicators
   return(list(delta   = delta,
               logitPi = logitPi,
               accept  = accept))
}


# sample each element of gamma through Metropolis step
MH_gamma = function(gamma, tau, x, logitPi, logLambda, phi_gamma, nu)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = rep(0, p)
   
   for (j in 1 : p)
   {
      # get data, logit(pi) and log(lambda) for the node j
      x_j         = x[ , j]
      logitPi_j   = logitPi[ , j]
      logLambda_j = logLambda[ , j]
      
      # current value of \gamma_j
      gamma_old = gamma[j]
      
      # propose new value of \gamma_j using Normal proposal distribution
      gamma_new     = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_gamma))
      
      # calculate MH ratio
      logLambda_new = logLambda_j + (gamma_new - gamma_old)
      llik_old      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
      llik_new      = llik_ZIPBN_j(x_j, logitPi_j, logLambda_new)
      ratio_MH      = exp(sum(llik_new - llik_old) - 0.5 * tau[4] * (gamma_new * gamma_new - gamma_old * gamma_old))
      
      # accept the proposed value with probability of min(1, MH ratio)
      if (runif(1) < min(1, ratio_MH))
      {
         gamma[j]        = gamma_new
         logLambda[ , j] = logLambda_new   # update and save log(lambda) for the node j
         accept[j]       = 1               # 1 if proposal accepted
      }
   }
   
   # return the updated gamma and log(lambda), and the acceptance indicators
   return(list(gamma     = gamma,
               logLambda = logLambda,
               accept    = accept))
}


# sample each element of A through Metropolis step (propose addition or deletion of an edge)
# corresponding alpha, beta, delta, and gamma are jointly proposed with A 
MH_A_each = function(A, alpha, beta, delta, gamma, tau, rho, x, logitPi, logLambda, phi_A, nu, verbose)
{
   # get the number of nodes
   p = ncol(x)
   
   for (j in 1 : p)
   {
      # get data, logit(pi) and log(lambda) for the node j
      x_j = x[ , j]
      logitPi_j   = logitPi[ , j]
      logLambda_j = logLambda[ , j]
      
      for (k in 1 : p)
      {
         # if j = k, do not sample 
         if (j == k) next
         
         # current value of \alpha_{jk}, \beta_{jk}, \delta_j, and \gamma_j
         alpha_old = alpha[j, k]
         beta_old  = beta[j, k]
         delta_old = delta[j]
         gamma_old = gamma[j]
         
         # divide into two cases: 1. there is currently no edge k -> j, 2. there exist an edge k -> j now 
         if (A[j, k] == 0)
         {
            # if there is no edge k -> j, propose addition of the edge unless it makes a cycle
            A[j, k] = A_new = 1
            graph   = graph_from_adjacency_matrix(A)
            A[j, k] = 0
            if (!is_dag(graph)) next
            
            # propose new value of \alpha_{jk}, \beta_{jk}, \delta_j, and \gamma_j using Normal proposal distribution
            alpha_new = rnorm(1, mean = alpha_old, sd = sqrt(1 / phi_A[2]))
            beta_new  = rnorm(1, mean = beta_old, sd = sqrt(1 / phi_A[3]))
            delta_new = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_A[4]))
            gamma_new = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_A[5]))
            
            # calculate MH ratio
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
            
            # accept the proposed values with probabiliof min(1, MH ratio)ty 
            if (runif(1) < min(1, ratio_MH))
            {
               A[j, k]      = A_new
               alpha[j, k]  = alpha_new
               beta[j, k]   = beta_new
               delta[j]     = delta_new
               gamma[j]     = gamma_new
               logitPi_j    = logitPi_new     # update logit(pi) for the node j
               logLambda_j  = logLambda_new   # update log(lambda) for the node j
               
               # print addition of the edge if verbose = TRUE
               if (verbose)
                  cat("An edge", k, "->", j, "is added \n")
            }
         } 
         else
         {
            # if there is an edge k -> j, propose deletion of the edge
            A_new = 0
            
            # propose new value of \alpha_{jk}, \beta_{jk}, \delta_j, and \gamma_j using Normal proposal distribution
            # Normal proposal distributions for \alpha_{jk} and \beta_{jk} have mean 0
            alpha_new = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
            beta_new  = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
            delta_new = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_A[4]))
            gamma_new = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_A[5]))
            
            # calculate MH ratio
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
            
            # accept the proposed value with probability of min(1, MH ratio)
            if (runif(1) < min(1, ratio_MH))
            {
               A[j, k]      = A_new
               alpha[j, k]  = alpha_new
               beta[j, k]   = beta_new
               delta[j]     = delta_new
               gamma[j]     = gamma_new
               logitPi_j    = logitPi_new     # update logit(pi) for the node j
               logLambda_j  = logLambda_new   # update log(lambda) for the node j
               
               # print deletion of the edge if verbose = TRUE
               if (verbose)
                  cat("An edge", k, "->", j, "is deleted \n")
            }
         }
      }
      
      # save the updated logit(pi) and log(lambda) for the node j
      logitPi[ , j]   = logitPi_j
      logLambda[ , j] = logLambda_j
   }
   
   # return the updated A, alpha, beta, delta, gamma,logit(pi), and log(lambda) 
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
MH_A_rev = function(A, alpha, beta, delta, gamma, tau, x, logitPi, logLambda, phi_A, nu, verbose)
{
   # get indices of existing edges 
   ids_A1 = which(A == 1, arr.ind = TRUE)
   n_swap = nrow(ids_A1)
   
   # if there exists no edge, don't do anything
   if (n_swap > 0)
   {
      for (s in 1 : n_swap)
      {
         # index of an edge which will be reversed
         j1 = k0 = ids_A1[s, 1]
         k1 = j0 = ids_A1[s, 2]
         
         # propose reversal of the edge unless it makes a cycle
         A[j1, k1] = 0
         A[j0, k0] = 1
         graph    = graph_from_adjacency_matrix(A)
         A[j1, k1] = 1
         A[j0, k0] = 0
         if (!is_dag(graph)) next
         
         # get data, logit(pi) and log(lambda) for reversal of the edge
         x_j1 = x_k0 = x[ , j1]
         x_j0 = x_k1 = x[ , j0]
         logitPi_j1   = logitPi[ , j1]
         logitPi_j0   = logitPi[ , j0]
         logLambda_j1 = logLambda[ , j1]
         logLambda_j0 = logLambda[ , j0]
         
         # current values of elements of alpha, beta, delta, and gamma corresponding to reversing
         alpha_old1 = alpha[j1, k1]
         alpha_old0 = alpha[j0, k0]
         beta_old1  = beta[j1, k1]
         beta_old0  = beta[j0, k0]
         delta_old1 = delta[j1]
         delta_old0 = delta[j0]
         gamma_old1 = gamma[j1]
         gamma_old0 = gamma[j0]
         
         # propose new values of elements of alpha, beta, delta, and gamma corresponding to reversing
         alpha_new1 = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
         alpha_new0 = rnorm(1, mean = alpha_old0, sd = sqrt(1 / phi_A[2]))
         beta_new1  = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
         beta_new0  = rnorm(1, mean = beta_old0, sd = sqrt(1 / phi_A[3]))
         delta_new1 = rnorm(1, mean = delta_old1, sd = sqrt(1 / phi_A[4]))
         delta_new0 = rnorm(1, mean = delta_old0, sd = sqrt(1 / phi_A[4]))
         gamma_new1 = rnorm(1, mean = gamma_old1, sd = sqrt(1 / phi_A[5]))
         gamma_new0 = rnorm(1, mean = gamma_old0, sd = sqrt(1 / phi_A[5]))
         
         # calculate MH ratio
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
         
         # accept the proposed value with probability of min(1, MH ratio)
         if (runif(1) < min(1, ratio_MH))
         {
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
            logitPi[ , j0]   = logitPi_new0     # update logit(pi)'s, following reversal of the edge
            logLambda[ , j1] = logLambda_new1
            logLambda[ , j0] = logLambda_new0   # update log(lambda)'s, following reversal of the edge
            
            # print reversal of the edge if verbose = TRUE
            if (verbose)
               cat("An edge", k1, "->", j1, "is reversed to", k0, "->", j0, " \n")
         } 
      } 
   }
   
   # return the updated A, alpha, beta, delta, gamma,logit(pi), and log(lambda)
   return(list(A         = A,
               alpha     = alpha,
               beta      = beta,
               delta     = delta,
               gamma     = gamma,
               logitPi   = logitPi,
               logLambda = logLambda))
}