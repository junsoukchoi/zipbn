# NOTE: The following codes use "A" instead of "E" to denote the graph structure of ZIPBN
# update a Markov chain once according to MH algorithm
update_chain = function(data, param, tuning, priors, temp = 1, prob_rev = 0.5, verbose = FALSE)
{
   # parameters to be updated
   alpha = param$alpha
   beta  = param$beta
   delta = param$delta
   gamma = param$gamma
   A     = param$E               
   tau   = param$tau
   rho   = param$rho
   logitPi   = param$logitPi
   logLambda = param$logLambda
   
   # tuning parameters
   phi_alpha = 1 / (tuning$sigma_alpha)^2
   phi_beta  = 1 / (tuning$sigma_beta)^2
   phi_delta = 1 / (tuning$sigma_delta)^2
   phi_gamma = 1 / (tuning$sigma_gamma)^2
   phi_A     = 1 / (tuning$sigma_E)^2
   phi_A_rev = 1 / (tuning$sigma_E_rev)^2
   
   # update alpha 
   update1 = update_alpha(data, alpha, A, tau, logitPi, logLambda, phi_alpha, temp)
   alpha   = update1$alpha
   logitPi = update1$logitPi
   
   # update beta
   update2   = update_beta(data, beta, A, tau, logitPi, logLambda, phi_beta, temp)
   beta      = update2$beta
   logLambda = update2$logLambda
   
   # update delta
   update3 = update_delta(data, delta, tau, logitPi, logLambda, phi_delta, temp)
   delta   = update3$delta
   logitPi = update3$logitPi
   
   # update gamma 
   update4   = update_gamma(data, gamma, tau, logitPi, logLambda, phi_gamma, temp)
   gamma     = update4$gamma
   logLambda = update4$logLambda
   
   # update A
   if (runif(1) > prob_rev)
      update5 = update_A_each(data, alpha, beta, delta, gamma, A, tau, rho, logitPi, logLambda, phi_A, temp, verbose)
   else
      update5 = update_A_rev(data, alpha, beta, delta, gamma, A, tau, logitPi, logLambda, phi_A_rev, temp, verbose)
   A         = update5$A
   alpha     = update5$alpha
   beta      = update5$beta
   delta     = update5$delta
   gamma     = update5$gamma
   logitPi   = update5$logitPi
   logLambda = update5$logLambda
   
   # update tau
   tau = update_tau(alpha, beta, delta, gamma, A, priors$tau[1], priors$tau[2], temp)
   
   # uphdate rho
   rho = update_rho(A, priors$rho[1], priors$rho[2], temp)
   
   # return updated parameters and acceptance record
   result = list()
   result$param = list()
   result$param$alpha = alpha
   result$param$beta  = beta
   result$param$delta = delta
   result$param$gamma = gamma
   result$param$E     = A
   result$param$tau   = tau
   result$param$rho   = rho
   result$param$logitPi   = logitPi
   result$param$logLambda = logLambda
   
   result$acceptance = list()
   result$acceptance$alpha = update1$accept
   result$acceptance$beta  = update2$accept
   result$acceptance$delta = update3$accept
   result$acceptance$gamma = update4$accept
   
   return(result)
}

# compute log-joint density of ZIPBN models (needed for the parallel tempering)
log_dZIPBN = function(data, param, priors, temp = 1)
{
   p = ncol(data)
   
   # parameters 
   alpha = param$alpha
   beta  = param$beta
   delta = param$delta
   gamma = param$gamma
   A     = param$E
   tau   = param$tau
   rho   = param$rho
   logitPi   = param$logitPi
   logLambda = param$logLambda

   # calculate the likelihood
   llik = 0
   for (j in 1 : p)
      llik = llik + sum(llik_ZIPBN_j(data[ , j], logitPi[ , j], logLambda[ , j]))
   
   # calculate the prior density
   lprior = sum(dnorm(alpha[A == 1], mean = 0, sd = sqrt(1 / tau[1]), log = TRUE)) + 
      sum(dnorm(beta[A == 1], mean = 0, sd = sqrt(1 / tau[2]), log = TRUE)) + 
      sum(dnorm(delta, mean = 0, sd = sqrt(1 / tau[3]), log = TRUE)) + 
      sum(dnorm(gamma, mean = 0, sd = sqrt(1 / tau[4]), log = TRUE)) + 
      sum(A) * log(rho) + (p * (p - 1) - sum(A)) * log(1 - rho) + 
      dgamma(tau[1], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) + 
      dgamma(tau[2], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) + 
      dgamma(tau[3], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) + 
      dgamma(tau[4], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) + 
      dbeta(rho, shape1 = priors$rho[1], shape2 = priors$rho[2], log = TRUE)
   
   return((llik + lprior) / temp)
}

# evaluate log-likelihood of each observation for the j-th component of ZIPBN model
llik_ZIPBN_j = function(x, logitPi, logLambda)
{
   # calculate pi and lambda
   pi     = exp(logitPi) / (1 + exp(logitPi))
   lambda = exp(logLambda)
   pi[is.nan(pi)] = 1
   lambda[lambda == Inf] = .Machine$double.xmax
   
   # evaluate and return log-likelihood of each observation
   llik = rep(0, length(x))
   llik[x == 0] = log(pi[x == 0] + (1 - pi[x == 0]) * dpois(0, lambda[x == 0]))
   llik[x > 0]  = log(1 - pi[x > 0]) + dpois(x[x > 0], lambda[x > 0], log = TRUE)
   return(llik)
}

# update each element of alpha through Metropolis-Hastings step
update_alpha = function(x, alpha, A, tau, logitPi, logLambda, phi_alpha, temp = 1)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = matrix(NA, p, p)
   
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
         
         if (A[j, k] == 1)
         {
            accept[j, k] = 0
            
            # current value of \alpha_{jk}
            alpha_old = alpha[j, k]
            
            # propose new value of \alpha_{jk} using Normal proposal distribution
            alpha_new   = rnorm(1, mean = alpha_old, sd = sqrt(1 / phi_alpha))
            logitPi_new = logitPi_j + x[ , k] * (alpha_new - alpha_old)
            
            # calculate MH ratio
            llik_old = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
            llik_new = llik_ZIPBN_j(x_j, logitPi_new, logLambda_j)
            ratio_MH = exp((sum(llik_new - llik_old) - 0.5 * tau[1] * (alpha_new * alpha_new - alpha_old * alpha_old)) / temp)
         
            # accept the proposed value with probability of min(1, MH ratio)
            if (is.nan(ratio_MH)) ratio_MH = 0
            if (runif(1) < min(1, ratio_MH))
            {
               alpha[j, k]  = alpha_new
               logitPi_j    = logitPi_new   # update logit(pi) for the node j 
               accept[j, k] = 1             # 1 if proposal accepted
            }
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

# update each element of beta through Metropolis-Hastings step
update_beta = function(x, beta, A, tau, logitPi, logLambda, phi_beta, temp = 1)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = matrix(NA, p, p)
   
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
         
         if(A[j, k] == 1)
         {
            accept[j, k] = 0
            
            # current value of \beta_{jk}
            beta_old = beta[j, k]
      
            # propose new value of \beta_{jk} using Normal proposal distribution
            beta_new      = rnorm(1, mean = beta_old, sd = sqrt(1 / phi_beta))
            logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old)
            
            # calculate MH ratio
            llik_old = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
            llik_new = llik_ZIPBN_j(x_j, logitPi_j, logLambda_new)
            ratio_MH = exp((sum(llik_new - llik_old) - 0.5 * tau[2] * (beta_new * beta_new - beta_old * beta_old)) / temp)
            
            # accept the proposed value with probability of min(1, MH ratio)
            if (is.nan(ratio_MH)) ratio_MH = 0
            if (runif(1) < min(1, ratio_MH))
            {
               beta[j, k]   = beta_new
               logLambda_j  = logLambda_new   # update log(lambda) for the node j
               accept[j, k] = 1               # 1 if proposal accepted
            }
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

# update each element of delta through Metropolis-Hastings step
update_delta = function(x, delta, tau, logitPi, logLambda, phi_delta, temp = 1)
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
      logitPi_new = logitPi_j + (delta_new - delta_old)
      
      # calculate MH ratio
      llik_old = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
      llik_new = llik_ZIPBN_j(x_j, logitPi_new, logLambda_j)
      ratio_MH = exp((sum(llik_new - llik_old) - 0.5 * tau[3] * (delta_new * delta_new - delta_old * delta_old)) / temp)
      
      # accept the proposed value with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
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

# update each element of gamma through Metropolis-Hastings step
update_gamma = function(x, gamma, tau, logitPi, logLambda, phi_gamma, temp = 1)
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
      logLambda_new = logLambda_j + (gamma_new - gamma_old)
      
      # calculate MH ratio
      llik_old = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
      llik_new = llik_ZIPBN_j(x_j, logitPi_j, logLambda_new)
      ratio_MH = exp((sum(llik_new - llik_old) - 0.5 * tau[4] * (gamma_new * gamma_new - gamma_old * gamma_old)) / temp)
      
      # accept the proposed value with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
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

# update each element of A through Metropolis-Hastings step (propose addition or deletion of an edge)
# corresponding alpha, beta, delta, and gamma are jointly proposed with A 
update_A_each = function(x, alpha, beta, delta, gamma, A, tau, rho, logitPi, logLambda, phi_A, temp = 1, verbose = TRUE)
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
            
            # propose new value of \alpha_{jk}, \beta_{jk}, \delta_j, and \gamma_j 
            alpha_new = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
            beta_new  = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[2]))
            delta_new = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_A[3]))
            gamma_new = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_A[4]))
            logitPi_new   = logitPi_j + x[ , k] * (alpha_new - alpha_old) + (delta_new - delta_old)
            logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old) + (gamma_new - gamma_old)
            
            # calculate MH ratio
            llik_old   = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
            llik_new   = llik_ZIPBN_j(x_j, logitPi_new, logLambda_new)
            ratio_post = sum(llik_new - llik_old) + 
               0.5 * log(tau[1]) - 0.5 * tau[1] * alpha_new * alpha_new + 
               0.5 * log(tau[2]) - 0.5 * tau[2] * beta_new * beta_new -
               0.5 * tau[3] * (delta_new * delta_new - delta_old * delta_old) - 
               0.5 * tau[4] * (gamma_new * gamma_new - gamma_old * gamma_old) +
               log(rho) - log(1 - rho)
            ratio_prop = - 0.5 * log(phi_A[1]) + 0.5 * phi_A[1] * alpha_new * alpha_new - 
               0.5 * log(phi_A[2]) + 0.5 * phi_A[2] * beta_new * beta_new
            ratio_MH = exp(ratio_post / temp + ratio_prop)
            
            # accept the proposed values with probabiliof min(1, MH ratio) 
            if (is.nan(ratio_MH)) ratio_MH = 0
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
                  cat("An edge", j, "<-", k, "is added \n")
            }
         } 
         else
         {
            # if there is an edge k -> j, propose deletion of the edge
            A_new = 0
            
            # propose new value of \alpha_{jk}, \beta_{jk}, \delta_j, and \gamma_j 
            alpha_new = beta_new = 0
            delta_new = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_A[3]))
            gamma_new = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_A[4]))
            logitPi_new   = logitPi_j + x[ , k] * (alpha_new - alpha_old) + (delta_new - delta_old)
            logLambda_new = logLambda_j + x[ , k] * (beta_new - beta_old) + (gamma_new - gamma_old)
            
            # calculate MH ratio
            llik_old   = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
            llik_new   = llik_ZIPBN_j(x_j, logitPi_new, logLambda_new)
            ratio_post = sum(llik_new - llik_old) - 
               0.5 * log(tau[1]) + 0.5 * tau[1] * alpha_old * alpha_old - 
               0.5 * log(tau[2]) + 0.5 * tau[2] * beta_old * beta_old - 
               0.5 * tau[3] * (delta_new * delta_new - delta_old * delta_old) - 
               0.5 * tau[4] * (gamma_new * gamma_new - gamma_old * gamma_old) +
               log(1 - rho) - log(rho)
            ratio_prop = 0.5 * log(phi_A[1]) - 0.5 * phi_A[1] * alpha_old * alpha_old + 
               0.5 * log(phi_A[2]) - 0.5 * phi_A[2] * beta_old * beta_old
            ratio_MH = exp(ratio_post / temp + ratio_prop)
            
            # accept the proposed value with probability of min(1, MH ratio)
            if (is.nan(ratio_MH)) ratio_MH = 0
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
                  cat("An edge", j, "<-", k, "is deleted \n")
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

# update A based on proposal of reversing an edge through Metropolis-Hastings step 
# corresponding alpha, beta, delta, and gamma are jointly proposed with A 
update_A_rev = function(x, alpha, beta, delta, gamma, A, tau, logitPi, logLambda, phi_A, temp = 1, verbose = TRUE)
{
   # get indices of existing edges 
   id_edges = which(A == 1, arr.ind = TRUE)
   n_rev    = nrow(id_edges)
   
   # if there exists no edge, don't do anything
   if (n_rev > 0)
   {
      for (s in 1 : n_rev)
      {
         # index of an edge which will be reversed
         j = id_edges[s, 1]
         k = id_edges[s, 2]
         
         # propose reversal of the edge unless it makes a cycle
         A[j, k] = 0
         A[k, j] = 1
         graph   = graph_from_adjacency_matrix(A)
         A[j, k] = 1
         A[k, j] = 0
         if (!is_dag(graph)) next
         
         # get data, logit(pi) and log(lambda) for reversal of the edge
         x_j = x[ , j]
         x_k = x[ , k]
         logitPi_j   = logitPi[ , j]
         logitPi_k   = logitPi[ , k]
         logLambda_j = logLambda[ , j]
         logLambda_k = logLambda[ , k]
         
         # current values of elements of alpha, beta, delta, and gamma corresponding to reversing
         alpha_jk_old = alpha[j, k]
         alpha_kj_old = alpha[k, j]
         beta_jk_old  = beta[j, k]
         beta_kj_old  = beta[k, j]
         delta_j_old  = delta[j]
         delta_k_old  = delta[k]
         gamma_j_old  = gamma[j]
         gamma_k_old  = gamma[k]
         
         # propose new values of elements of alpha, beta, delta, and gamma corresponding to reversing
         alpha_jk_new = beta_jk_new  = 0
         alpha_kj_new = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
         beta_kj_new  = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[2]))
         delta_j_new  = rnorm(1, mean = delta_j_old, sd = sqrt(1 / phi_A[3]))
         delta_k_new  = rnorm(1, mean = delta_k_old, sd = sqrt(1 / phi_A[3]))
         gamma_j_new  = rnorm(1, mean = gamma_j_old, sd = sqrt(1 / phi_A[4]))
         gamma_k_new  = rnorm(1, mean = gamma_k_old, sd = sqrt(1 / phi_A[4]))
         logitPi_j_new   = logitPi_j + x_k * (alpha_jk_new - alpha_jk_old) + (delta_j_new - delta_j_old)
         logitPi_k_new   = logitPi_k + x_j * (alpha_kj_new - alpha_kj_old) + (delta_k_new - delta_k_old)
         logLambda_j_new = logLambda_j + x_k * (beta_jk_new - beta_jk_old) + (gamma_j_new - gamma_j_old)
         logLambda_k_new = logLambda_k + x_j * (beta_kj_new - beta_kj_old) + (gamma_k_new - gamma_k_old)
         
         # calculate MH ratio
         llik_j_old = llik_ZIPBN_j(x_j, logitPi_j, logLambda_j)
         llik_k_old = llik_ZIPBN_j(x_k, logitPi_k, logLambda_k)
         llik_j_new = llik_ZIPBN_j(x_j, logitPi_j_new, logLambda_j_new)
         llik_k_new = llik_ZIPBN_j(x_k, logitPi_k_new, logLambda_k_new)
         ratio_post = sum(llik_j_new - llik_j_old) + sum(llik_k_new - llik_k_old) - 
            0.5 * tau[1] * (alpha_kj_new * alpha_kj_new - alpha_jk_old * alpha_jk_old) - 
            0.5 * tau[2] * (beta_kj_new * beta_kj_new - beta_jk_old * beta_jk_old) -
            0.5 * tau[3] * (delta_j_new * delta_j_new - delta_j_old * delta_j_old + 
                               delta_k_new * delta_k_new - delta_k_old * delta_k_old) - 
            0.5 * tau[4] * (gamma_j_new * gamma_j_new - gamma_j_old * gamma_j_old +
                               gamma_k_new * gamma_k_new - gamma_k_old * gamma_k_old)
         ratio_prop = - 0.5 * phi_A[1] * (alpha_jk_old * alpha_jk_old - alpha_kj_new * alpha_kj_new) - 
            0.5 * phi_A[2] * (beta_jk_old * beta_jk_old - beta_kj_new * beta_kj_new)
         ratio_MH   = exp(ratio_post / temp + ratio_prop)
         
         # accept the proposed value with probability of min(1, MH ratio)
         if (is.nan(ratio_MH)) ratio_MH = 0   
         if (runif(1) < min(1, ratio_MH))
         {
            A[j, k] = 0
            A[k, j] = 1
            alpha[j, k] = alpha_jk_new
            alpha[k, j] = alpha_kj_new
            beta[j, k]  = beta_jk_new
            beta[k, j]  = beta_kj_new
            delta[j]    = delta_j_new
            delta[k]    = delta_k_new
            gamma[j]    = gamma_j_new
            gamma[k]    = gamma_k_new
            logitPi[ , j]   = logitPi_j_new
            logitPi[ , k]   = logitPi_k_new     # update logit(pi)'s, following reversal of the edge
            logLambda[ , j] = logLambda_j_new
            logLambda[ , k] = logLambda_k_new   # update log(lambda)'s, following reversal of the edge
            
            # print reversal of the edge if verbose = TRUE
            if (verbose)
               cat("An edge", j, "<-", k, "is reversed to", k, "<-", j, " \n")
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

# update tau via Gibbs sampling
update_tau = function(alpha, beta, delta, gamma, A, a, b, temp = 1)
{
   p   = nrow(A)
   tau = rep(NA, 4)
   
   # sample tau_alpha
   tau[1] = rgamma(1, shape = (a + 0.5 * sum(A) + temp - 1) / temp, rate = (b + 0.5 * sum(alpha * alpha)) / temp)
   
   # sample tau_beta
   tau[2] = rgamma(1, shape = (a + 0.5 * sum(A) + temp - 1) / temp, rate = (b + 0.5 * sum(beta * beta)) / temp)
   
   # sample tau_delta
   tau[3] = rgamma(1, shape = (a + 0.5 * p + temp - 1) / temp, rate = (b + 0.5 * sum(delta * delta)) / temp)
   
   # sample tau_gamma
   tau[4] = rgamma(1, shape = (a + 0.5 * p + temp - 1) / temp, rate = (b + 0.5 * sum(gamma * gamma)) / temp)
   
   # return the sampled tau
   return(tau)
}

# update rho via Gibbs sampling
update_rho = function(A, a, b, temp = 1)
{
   p = nrow(A)
   
   # sample rho
   rho = rbeta(1, shape1 = (a + sum(A) + temp - 1) / temp, shape2 = (b + p * (p - 1) - sum(A) + temp - 1) / temp)
   
   # return the sampled rho
   return(rho)
}

