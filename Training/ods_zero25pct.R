# clear
rm(list = ls())
# set random seed
set.seed(05012020)
# load necessary libraries
library(foreach)
library(doParallel)
library(doRNG)
library(igraph)
library(glmnet)

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

alpha_true[E_true == 1] = runif(n_edge, -1.0, -0.3)
beta_true[E_true == 1]  = runif(n_edge, -0.4, -0.2)
delta_true[ ] = runif(p, -2, -1.5)
gamma_true[ ] = runif(p, 1.5, 2)

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
# implementation of [Park & Raskutti, 2015]
# tuning parameter
c0 = 0.004

cluster = makeCluster(10)
registerDoParallel(cluster)
ods = foreach(iter = 1 : n_rep, .packages = "glmnet") %dorng%
{
   data = dat[[iter]]
   
   # step 1: estimate the undirected edges corresponding to the moralized graph with neighborhood sets
   moral_graph = matrix(0, p, p)
   for (j in 1 : p)
   {
      cv  = cv.glmnet(data[ , -j], data[ , j], family = "poisson")
      out = glmnet(data[ , -j], data[ , j], family = "poisson", lambda = cv$lambda.1se)
      moral_graph[j, -j] = as.vector(out$beta != 0) + 0
   }
   moral_graph = ((moral_graph != 0) & (t(moral_graph) != 0)) + 0
   
   # step 2: estimate causal ordering using overdispersion score
   dag = matrix(0, p, p)
   ordering = rep(NA, p)
   
   # 1st element of the causal ordering
   score = apply(data, 2, var) - colMeans(data)
   ordering[1] = which.min(score)
   
   # 2nd to (p-1)-th elements of the causal ordering
   for (j in 2 : (p - 1))
   {
      score = rep(0, p)
      score[ordering[1 : (j - 1)]] = Inf
      
      for (k in (1 : p)[-ordering[1 : (j - 1)]])
      {
         x_k = data[ , k]
         parents = intersect(which(moral_graph[k, ] != 0), ordering[1 : (j - 1)])
         
         if (length(parents) == 0)
         {
            score[k] = var(x_k) - mean(x_k)
         } else
         {
            id      = apply(cbind(rep("x", n), data[ , parents]) , 1 , paste, collapse = "-")
            uniq_id = unique(id)
            num = denom = 0
            
            for (l in 1 : length(uniq_id))
            {
               n_Xs = sum(id == uniq_id[l])
               if (n_Xs < c0 * n)
                  next
               num   = num + n_Xs * (var(x_k[id == uniq_id[l]]) - mean(x_k[id == uniq_id[l]]))
               denom = denom + n_Xs 
            }

            score[k] = num / denom
         }
      }
      
      pi_j = which.min(score)
      ordering[j] = pi_j
      
      # step 3: estimate directed edges toward the j-th element of the causal ordering
      candidates = intersect(which(moral_graph[pi_j, ] != 0), ordering[1 : (j - 1)])
      
      if (length(candidates) <= 1)
      {
         dag[pi_j, candidates] = 1
      } else 
      {
         cv  = cv.glmnet(data[ , candidates, drop = FALSE], data[ , pi_j], family = "poisson")
         out = glmnet(data[ , candidates], data[ , pi_j], family = "poisson", lambda = cv$lambda.1se)
         dag[pi_j, candidates] = as.vector(out$beta != 0) + 0
      }
   }
   
   # p-th element of the causal ordering
   pi_p = (1 : p)[-ordering[1 : (p - 1)]]
   ordering[p] = pi_p
   
   # directed edges toward the p-th element of the causal ordering
   candidates = intersect(which(moral_graph[pi_p, ] != 0), ordering[1 : (p - 1)])
   if (length(candidates) <= 1)
   {
      dag[pi_j, candidates] = 1
   } else 
   {
      cv  = cv.glmnet(data[ , candidates], data[ , pi_p], family = "poisson")
      out = glmnet(data[ , candidates], data[ , pi_p], family = "poisson", lambda = cv$lambda.1se)
      dag[pi_p, candidates] = as.vector(out$beta != 0) + 0
   }
   
   # print progress
   cat("iter=", iter, "\n")
   
   # return the estimated dag
   dag
}

stopCluster(cluster)

save(ods, file = "ods_z25.RData")
