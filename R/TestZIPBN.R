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
for (i in 1 : 5)
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