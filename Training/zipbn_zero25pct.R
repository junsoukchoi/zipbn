# clear
rm(list = ls())
# set random seed
set.seed(05012020)
# load necessary libraries
library(foreach)
library(doParallel)
library(doRNG)
library(igraph)
library(pscl)
# load necessary functions
source("ZIPBNfunctions.R")

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
# run the parallel-tempered MCMC 
cluster = makeCluster(30)
registerDoParallel(cluster)

zipbn = foreach(rep = 1 : n_rep, .packages = c("igraph", "pscl")) %dorng%
{
   data = dat[[rep]]
      
   # necessary for the MCMC implementation
   MCMC      = 3000
   BURNIN    = 1500
   n_chains  = 10
   prob_swap = 0.1
   n = nrow(data); p = ncol(data)

   # sequence of temperatures for the parallel tempering
   temps = exp(seq(0, 1, length.out = n_chains))
   
   # hyperparameters for priors
   priors = list()
   priors$tau = c(0.01, 0.01)
   priors$rho = c(0.5, 0.5)

   # sd of the Metropolis sampler Normal proposal distribution
   tuning = list()
   for (m in 1 : n_chains)
   {
      tuning[[m]] = list()
      tuning[[m]]$sigma_alpha = 0.7  * 2^((m - 1) / (n_chains - 1))
      tuning[[m]]$sigma_beta  = 0.05 * 2^((m - 1) / (n_chains - 1))
      tuning[[m]]$sigma_delta = 0.7  * 2^((m - 1) / (n_chains - 1))
      tuning[[m]]$sigma_gamma = 0.1  * 2^((m - 1) / (n_chains - 1))
      tuning[[m]]$sigma_E     = c(0.05, 0.05, 0.7, 0.1)
      tuning[[m]]$sigma_E_rev = c(0.5, 0.5, 1.4, 0.2) 
   }
   
   # starting values for MCMC
   starting = list()
   starting$alpha = matrix(0, p, p)
   starting$beta  = matrix(0, p, p)
   starting$delta = rep(0, p)
   starting$gamma = rep(0, p)
   starting$E     = matrix(0, p, p)
   starting$tau   = c(10, 10, 1, 1)
   starting$rho   = 1 / p
   for (j in 1 : p)
   {
      mle = zeroinfl(data[ , j] ~ 1 | 1, dist = "poisson")
      starting$delta[j] = mle$coefficients$zero
      starting$gamma[j] = mle$coefficients$count
   }
   
   # calculate logitPi and logLambda for each group with starting values
   starting$logitPi   = tcrossprod(data, starting$alpha) + matrix(starting$delta, n, p, byrow = TRUE)
   starting$logLambda = tcrossprod(data, starting$beta) + matrix(starting$gamma, n, p, byrow = TRUE)
   
   # initialize parameters
   param = list()
   for (m in 1 : n_chains)
   {
      param[[m]] = starting
   }
   
   # initialize MCMC samples
   mcmc_samples = list()
   mcmc_samples$alpha = array(NA, dim = c(p, p, MCMC))
   mcmc_samples$beta  = array(NA, dim = c(p, p, MCMC))
   mcmc_samples$delta = matrix(NA, p, MCMC)
   mcmc_samples$gamma = matrix(NA, p, MCMC)
   mcmc_samples$E     = array(NA, dim = c(p, p, MCMC))
   mcmc_samples$tau   = matrix(NA, 4, MCMC)
   mcmc_samples$rho   = rep(NA, MCMC)

   # run MCMC iterations
   for (t in 1 : MCMC)
   {
      if (runif(1) > prob_swap)
      {
         # perform one-step update for all chains 
         for (m in 1 : n_chains)
         {
            # update each chain
            out = update_chain(data, param[[m]], tuning[[m]], priors, temps[m])
            param[[m]] = out$param
         }
      } else
      {
         # propose a swapping move
         # randomly choose chains to swap
         rand = sort(sample(1 : n_chains, 2))
         m1   = rand[1]
         m2   = rand[2]
         
         # calculate MH ratio for swapping
         ratio_swap = log_dZIPBN(data, param[[m1]], priors, temps[m2]) + 
            log_dZIPBN(data, param[[m2]], priors, temps[m1]) -
            log_dZIPBN(data, param[[m1]], priors, temps[m1]) - 
            log_dZIPBN(data, param[[m2]], priors, temps[m2])
         ratio_swap = exp(ratio_swap)
   
         # accept swapping
         if (is.nan(ratio_swap)) ratio_swap = 0
         if (runif(1) < min(1, ratio_swap))
         {
            cat("swap the states of chains", m1, "and", m2, "\n")
            param_m1    = param[[m1]]
            param[[m1]] = param[[m2]]
            param[[m2]] = param_m1
         }   
      }

      # save MCMC samples for the cold chain
      mcmc_samples$alpha[ , , t] = param[[1]]$alpha
      mcmc_samples$beta[ , , t]  = param[[1]]$beta
      mcmc_samples$delta[ , t]   = param[[1]]$delta
      mcmc_samples$gamma[ , t]   = param[[1]]$gamma
      mcmc_samples$E[ , , t]     = param[[1]]$E
      mcmc_samples$tau[ , t]     = param[[1]]$tau
      mcmc_samples$rho[t]        = param[[1]]$rho
      
      # print progress
      if (t %% 100 == 0)
         cat("iter =", t, "\n")
   }
   
   # return MCMC samples for the cold chain
   mcmc_samples
}

# save simulation results
save(zipbn, E_true, alpha_true, beta_true, delta_true, gamma_true, file = "zipbn_z25.RData")

