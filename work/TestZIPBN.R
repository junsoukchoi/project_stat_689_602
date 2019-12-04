# clear
rm(list = ls())
# set random seed
set.seed(7)

# initialize true parameters with zero matrices
n = 1000
p = 3
A_true     = matrix(0, p, p)
alpha_true = matrix(0, p, p)
beta_true  = matrix(0, p, p)
delta_true = rep(0, p)
gamma_true = rep(0, p)

# generate a random graph
n_edges = 2
while (sum(A_true == 1) < n_edges)
{
   id_edge = matrix(sample(1 : p, 2), ncol = 2)
   A_true[id_edge] = 1
   g_true = igraph::graph_from_adjacency_matrix(t(A_true))
   # if selected edge make a directed cycle, discard the edge 
   if (!(igraph::is_dag(g_true)))
      A_true[id_edge] = 0
}
A_true
sum(A_true)   # number of edges

# determine true parameters given graph A
alpha_true[A_true == 1] = 0.2
beta_true[A_true == 1]  = -0.2
delta_true[ ] = -1
gamma_true[ ] = 2

# generate data with true parameters
x = matrix(0, n, p)
g_true = igraph::graph_from_adjacency_matrix(t(A_true))
order_nodes = igraph::as_ids(igraph::topo_sort(g_true))
for (j in order_nodes)
{
   pi = exp(x %*% alpha_true[j, ] + delta_true[j])
   pi = pi / (1 + pi)
   pi[is.nan(pi)] = 1
   lambda  = exp(x %*% beta_true[j, ] + gamma_true[j])
   x[ , j] = rpois(n, lambda) * (1 - rbinom(n, 1, pi))
}
head(x)
table(x)
mean(x == 0)   # proprtion of zeros

# evaluate log-likelihood for ZIPBN with true parameter values
llik_true = 0
logitPi_true   = tcrossprod(x, alpha_true) + matrix(delta_true, n, p, byrow = TRUE)
logLambda_true = tcrossprod(x, beta_true) + matrix(gamma_true, n, p, byrow = TRUE)
for (j in 1 : p)
{
   llik_true = llik_true + sum(llik_ZIPBN_j(x[ , j], logitPi_true[ , j], logLambda_true[ , j]))
}
llik_true

# use the mcmc_ZIPBN function to sample parameters from posterior distributions for ZIPBN
starting = tuning = priors = list()

## set starting values with true parameter values
#starting$alpha = alpha_true
#starting$beta  = beta_true
#starting$delta = delta_true
#starting$gamma = gamma_true
#starting$A     = A_true
#starting$tau   = c(10, 10, 1, 1)
#starting$rho   = 0.1

# set starting values with zero matrices and vectors
m = colMeans(x)
v = apply(x, 2, var)
delta = (v - m) / (v - m + m * m)
delta = log(delta / (1 - delta))
gamma = (v - m + m * m) / m
gamma = log(gamma)
starting$alpha = matrix(0, p, p)
starting$beta  = matrix(0, p, p)
starting$delta = delta
starting$gamma = gamma
starting$A     = matrix(0, p, p)
starting$tau   = c(10, 10, 1, 1)
starting$rho   = 0.1

# set precision values for Metropolis sampler Normal proposal distribution
tuning$phi_alpha = c(100000000, 300)
tuning$phi_beta  = c(100000000, 300)
tuning$phi_delta = 20
tuning$phi_gamma = 400
tuning$phi_A     = c(100000000, 50, 50, 5, 50)

# set hyperparameter values
priors$nu = 10000^2
priors$b  = c(0.01, 0.01, 0.01, 0.01, 0.5)
priors$c  = c(0.01, 0.01, 0.01, 0.01, 0.5)

# run the mcmc_ZIPBN function
out = mcmc_ZIPBN(x, starting, tuning, priors)
apply(out$alpha, c(1, 2), mean)
apply(out$beta, c(1, 2), mean)
apply(out$delta, 1, mean)
apply(out$gamma, 1, mean)
apply(out$A, c(1, 2), mean)
apply(out$tau, 1, mean)
mean(out$rho)
