# clear
rm(list = ls())
# set random seed
set.seed(05012020)
# load necessary libraries
library(igraph)
library(MXM)

#----------------------------------------------------------------------------------------------------#
# simulate data
n_rep = 30

# set true parameters
n = 500   # sample size
p = 50    # number of node
E_true     = matrix(0, p, p)
alpha_true = matrix(0, p, p)
beta_true  = matrix(0, p, p)
delta_true = rep(0, p)
gamma_true = rep(0, p)

# generate a DAG with p edges
n_edge = p
while (sum(E_true == 1) < n_edge)
{
   id_edge = matrix(sample(1 : p, 2), ncol = 2)
   E_true[id_edge] = 1
   g_true = graph_from_adjacency_matrix(t(E_true))
   if (!(is_dag(g_true)))
      E_true[id_edge] = 0
}

alpha_true[E_true == 1] = runif(n_edge, 0.3, 1)
beta_true[E_true == 1]  = runif(n_edge, -1, -0.3)
delta_true[ ] = runif(p, 0.5, 1)
gamma_true[ ] = runif(p, 1, 1.5)

# generate n_rep datasets
dat = list()
for (rep in 1 : n_rep)
{
   data = matrix(0, n, p)
   order_nodes = as_ids(topo_sort(g_true))
   for (j in order_nodes)
   {
      pi = exp(data %*% alpha_true[j, ] + delta_true[j])
      pi = pi / (1 + pi)
      pi[is.nan(pi)] = 1
      lambda  = exp(data %*% beta_true[j, ] + gamma_true[j])
      data[ , j] = rpois(n, lambda) * (1 - rbinom(n, 1, pi))
   }
   
   dat[[rep]] = data
}

#----------------------------------------------------------------------------------------------------#
# implementation of [Park & Park, 2019]
# assume hyper-Poisson 
# tuning parameter
Nmin = 2

# define the CMR function
CMRfunction = function(x, a, b)
{
   x * x * ((a + 1) / a) * (b / (b + 1))
}

# MRS algorithm with r = 2
mrs = list()
for (iter in 1 : n_rep)
{
   data = dat[[iter]]
   
   # step 1: estimate the skeleton of the graph 
   skeleton = mmhc.skel(data, test = "testIndPois", nc = 10)$G

   # step 2: estimate an ordering of the graph using the second moments ratio scores
   dag      = matrix(0, p, p)
   ordering = rep(NA, p)
   
   # 1st element of the ordering
   score = rep(0, p) 
   b     = apply(data, 2, var) / colMeans(data)
   for (j in 1 : p)
   {
      moment1  = mean(data[ , j])
      moment2  = mean(data[ , j] * data[ , j])
      score[j] = moment2 / (CMRfunction(moment1, 1, b[j]) + moment1)
   }
   ordering[1] = which.min(score)
   
   # 2nd to (p-1)-th elements of the ordering
   for (m in 2 : (p - 1))
   {
      score = rep(0, p)
      score[ordering[1 : (m - 1)]] = Inf
      
      for (j in (1 : p)[-ordering[1 : (m - 1)]])
      {
         x_j     = data[ , j]
         parents = intersect(which(skeleton[j, ] != 0), ordering[1 : (m - 1)])
         
         if (length(parents) == 0)
         {
            moment1  = mean(x_j)
            moment2  = mean(x_j * x_j)
            score[j] = moment2 / (CMRfunction(moment1, 1, b[j]) + moment1)
         } else
         {
            id      = apply(cbind(rep("x", n), data[ , parents]) , 1 , paste, collapse = "-")
            uniq_id = unique(id)
            num = denom = 0
            
            for (s in 1 : length(uniq_id))
            {
               n_Xs = sum(id == uniq_id[s])
               if (n_Xs < Nmin)
                  next
              
               denom   = denom + n_Xs 
               moment1 = mean(x_j[id == uniq_id[s]])
               moment2 = mean(x_j[id == uniq_id[s]] * x_j[id == uniq_id[s]])
               if (moment1 == 0)
                  num = num + n_Xs
               else
                  num = num + n_Xs * moment2 / (CMRfunction(moment1, 1, b[j]) + moment1)
            }

            score[j] = num / denom
         }
      }
      
      ordering[m] = which.min(score)
   }
   
   # last element of the ordering
   ordering[p] = (1 : p)[-ordering[1 : (p - 1)]]
   
   # estimate the edge sets
   for (m in 2 : p)
   {
      edges = intersect(which(skeleton[ordering[m], ] != 0), ordering[1 : (m - 1)])
      dag[ordering[m], edges] = 1
   }
   
   # save the estimated dag
   mrs[[iter]] = dag
   
   # print progress
   cat("iter=", iter, "\n")
}

save(mrs, file = "mrs_z75.RData")
