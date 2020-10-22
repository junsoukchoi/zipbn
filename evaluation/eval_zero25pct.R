# clear
rm(list = ls())

#----------------------------------------------------------------------------------------------------#
# load ZIPBN results
load("pre-trained/zipbn_zero25pct.RData")

# load ODS results
load("pre-trained/ods_zero25pct.RData")

# load MRS results
load("pre-trained/mrs_zero25pct.RData")

#----------------------------------------------------------------------------------------------------#
# summarize simulation results
n_sim  = 30
MCMC   = 3000
BURNIN = 1500

ZIPBN = ODS = MRS = list()
ZIPBN$TPR = ODS$TPR = MRSTPR = rep(NA, n_sim)
ZIPBN$FDR = ODS$FDR = MRSFDR = rep(NA, n_sim)
ZIPBN$MCC = ODS$MCC = MRSMCC = rep(NA, n_sim)
for (iter in 1 : n_sim)
{
   ### evaluate ZIPBN results ###
   # burn-in
   E_MCMC     = zipbn[[iter]]$E[ , , (BURNIN + 1) : MCMC]
   alpha_MCMC = zipbn[[iter]]$alpha[ , , (BURNIN + 1) : MCMC]
   beta_MCMC  = zipbn[[iter]]$beta[ , , (BURNIN + 1) : MCMC]
   delta_MCMC = zipbn[[iter]]$delta[ , (BURNIN + 1) : MCMC]
   gamma_MCMC = zipbn[[iter]]$gamma[ , (BURNIN + 1) : MCMC]
   tau_MCMC   = zipbn[[iter]]$tau[ , (BURNIN + 1) : MCMC]
   rho_MCMC   = zipbn[[iter]]$rho[(BURNIN + 1) : MCMC]
   
   # estimate parameters for ZIPBN
   E_est     = (apply(E_MCMC, c(1, 2), mean) > 0.5) + 0
   alpha_est = apply(alpha_MCMC, c(1, 2), mean) * E_est
   beta_est  = apply(beta_MCMC, c(1, 2), mean) * E_est
   delta_est = rowMeans(delta_MCMC)
   gamma_est = rowMeans(gamma_MCMC)
   tau_est   = rowMeans(tau_MCMC)
   rho_est   = mean(rho_MCMC)
   
   # calculate TPR, FDR, and MCC from ZIPBN estimates for the true graph
   tab = table(E_est, E_true)
   TP  = tab[2, 2]
   TN  = tab[1, 1]
   FP  = tab[2, 1]
   FN  = tab[1, 2]
   ZIPBN$TPR[iter] = TP / (TP + FN)
   ZIPBN$FDR[iter] = FP / (FP + TP)
   ZIPBN$MCC[iter] = (TP * TN - FP * FN) / sqrt(TP + FP) / sqrt(TP + FN) / sqrt(TN + FP) / sqrt(TN + FN)
   
   ### evaluate ODS results ###
   # calculate TPR, FDR, and MCC from ODS estimates for the true graph
   tab = table(ods[[iter]], E_true)
   TP  = tab[2, 2]
   TN  = tab[1, 1]
   FP  = tab[2, 1]
   FN  = tab[1, 2]
   ODS$TPR[iter] = TP / (TP + FN)
   ODS$FDR[iter] = FP / (FP + TP)
   ODS$MCC[iter] = (TP * TN - FP * FN) / sqrt(TP + FP) / sqrt(TP + FN) / sqrt(TN + FP) / sqrt(TN + FN)
   
   ### evaluate MRS results ###
   # calculate TPR, FDR, and MCC from MRS estimates for the true graph
   tab = table(mrs[[iter]], E_true)
   TP  = tab[2, 2]
   TN  = tab[1, 1]
   FP  = tab[2, 1]
   FN  = tab[1, 2]
   MRS$TPR[iter] = TP / (TP + FN)
   MRS$FDR[iter] = FP / (FP + TP)
   MRS$MCC[iter] = (TP * TN - FP * FN) / sqrt(TP + FP) / sqrt(TP + FN) / sqrt(TN + FP) / sqrt(TN + FN)
}

# visualize the simulation results
par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
boxplot(ZIPBN$TPR, ODS$TPR, MRS$TPR, main = "TPR")
boxplot(ZIPBN$FDR, ODS$FDR, MRS$FDR, main = "FDR")
boxplot(ZIPBN$MCC, ODS$MCC, MRS$MCC, main = "MCC")
mtext("Simulations with ~25% zeros", outer = TRUE, cex = 1.2)
