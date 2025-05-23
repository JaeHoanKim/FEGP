library(Rcpp)
#########################################################################
########## Function for MCMC samples using ESS and WC algoritm ##########
#########################################################################

## knitting v (size N vector) and g_N (size N+1 vector) in an alternative order
sample.ESS.Nfixed2D = function(Z, X, l.in, nu.in, mcmc, brn, thin, sigsq, N.init, tausq = 1){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   ## l: range of the fields::Matern function
   if(length(Z) != dim(X)[1])
      stop("Z and X should have the same length!")
   em = mcmc + brn # total number of sampling
   ## compatibility check when predicting
   N_list = list()
   prob_list = c()
   g_list = list()
   # starting with the initial N and g
   # setting up for iterative work
   g.in = matrix(0, N.init+1, N.init + 1)
   N = N.init
   # since N does not change, we can fix 'result'.
   result = eigvals_exact(ndim = N, nu = nu.in, lambda_g = l.in, tausq)
   for(i in 1:em){
      ### update knot_N for updated N
      ## sample new g given N, D_n!= 102
      nu.ess = matrix(samp_from_grid(N, mdim = result$mvec, egs = result$egvals, nu, lambda_g), 
                      nrow = N[1] + 1, ncol = N[2]+1, byrow = TRUE)
      g.out = ESS(g.in, nu.ess, z = Z, x = X, sigsq)
      N.out = N
      # to make an output as a vector form! returning as a matrix occurs an error in glist_to_plotdf function
      g_list[[i]] = as.vector(t(g.out))
      # only take one N out of N.out, since two elements would be the same
      N_list[i] = N.out[1]
      ### preparation for iteration
      g.in = g.out
      N = N.out
      if (i%%100 == 0){
         print(c("interation number:", i))
      }
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list))
}



## goal: draw samples of (N, g^{(N)}, \theta) given the data D_n (X_n, Y_n)
## input: X, Y, sigma 
## output: mcmc samples after burn-in
## without considering the error term sigma

#### ESS nested method after multiple ESS samples when considering N to 2N
sample.ESS.nested2D = function(Z, X, N.pr, mcmc, brn, thin, Ndoubling, l.in = NULL, nu.in = NULL, sigsq, 
                               N.init, tausq){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   if(length(Z) != dim(X)[1])
      stop("Z and X should have the same length!")
   ## initial value of g
   g.init = matrix(0, N.init + 1, N.init + 1)
   em = mcmc + brn # total number of sampling
   ## compatibility check when predicting
   N_list = c()
   prob_list = c()
   g_list = list()
   # starting with the initial N and g
   g.in = g.init
   N = N.init
   # store result_N and result_2N to prevent re-calculation of the inverse of BTTB
   result_N_list = list()
   for(i in 1:5){
      N1 = 2 * 2^(i-1)
      N2 = 2 * 2^(i-1)
      X_N = cbind(rep(c(0:N1)/N1, each = N2 + 1), 
                  rep(c(0:N2)/N2, N1 + 1)) # (N1+1)*(N2+1) by 2 matrix
      Sigma_N = thekernel(X_N, nu.in, lambda = l.in)
      result_N_list[[N1]] = inv_BTTB(Sigma_N, N1 + 1)
   }
   # initial eigenvalue result
   result = eigvals_exact(ndim = N, nu = nu.in, lambda_g = l.in, tausq)
   for(i in 1:em){
      ### update knot_N for updated N
      ## sample new g given N, D_n
      uu = ((N[1] == 2) - 3) * runif(1)
      if (uu < -2){
         ## half N
         result_half = eigvals_exact(ndim = N / 2, nu = nu.in, lambda_g = l.in, tausq)
         g.cand.half = g.in[2 * c(0:(N[1]/2)) + 1, 2 * c(0:(N[2]/2)) + 1]
         log.test.prob = min(0, loglik(Z, X, g.cand.half, sigsq) - loglik(Z, X, g.in, sigsq) +
            log(N.pr(N[1]/2)) - log(N.pr(N[1])))
         u = runif(1)
         if (u < exp(log.test.prob)){
            g.out = g.cand.half
            N.out = N / 2
            result = result_half
         } else{
            g.out = g.in
            N.out = N
         }
      } else if(uu < -1){
         ## double N
         ## 1. sample from posterior distribution of g, (2N+1) by (2N+1) matrix
         result2 = eigvals_exact(ndim = 2 * N, nu = nu.in, lambda_g = l.in, tausq)
         gdouble = matrix(samp_from_grid(2*N, mdim = result2$mvec, egs = result2$egvals, nu.in, lambda_g = l.in), 
                          nrow = 2 * N[1] + 1, byrow = TRUE)
         for(kk in 1:Ndoubling){
            gdouble_ess = matrix(samp_from_grid(2*N, mdim = result2$mvec, egs = result2$egvals, nu.in, lambda_g = l.in), 
                                 nrow = 2 * N[1] + 1, byrow = TRUE)
            gdouble = ESS_double(gdouble, g.in, gdouble_ess, nu.in, l.in, 
                                  result_N = result_N_list[[nrow(g.in) - 1]], 
                                  result_2N = result_N_list[[nrow(gdouble) - 1]])
         }
         g.cand.2 = gdouble
         ## 2. likelihood check
         log.test.prob = min(0, loglik(Z, X, g.cand.2, sigsq) - loglik(Z, X, g.in, sigsq) +
            log(N.pr(N[1]*2)) - log(N.pr(N[1])) )
         u = runif(1)
         if (u < exp(log.test.prob)){
            g.out = g.cand.2
            N.out = 2 * N
            result = result2
         }else{
            g.out = g.in
            N.out = N
         }
      } else if (uu < 0){
         ## same N
         nu.ess = matrix(samp_from_grid(N, mdim = result$mvec, egs = result$egvals, nu, lambda_g), 
                         nrow = N[1] + 1, ncol = N[2]+1, byrow = TRUE)
         g.out = ESS(g.in, nu.ess, z = Z, x = X, sigsq)
         N.out = N
         log.test.prob = 3
      }
      ### preparation for iteration
      g.in = g.out
      N = N.out
      N_list[i] = N[1]
      g_list[[i]] = as.vector(t(g.out))
      if (i%%20 == 0){
         print(c("interation number:", i))
         print(N)
      }
      prob_list[i] = log.test.prob
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list, prob_list = prob_list))
}





sample.PTESS2D = function(Z, X, Nk, Tk, N.pr, mcmc, brn, thin, l.in = NULL, nu.in = NULL, sigsq, 
                            N.init, tausq){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   if(length(Z) != dim(X)[1])
      stop("Z and X should have the same length!")
   if(length(Nk)!=length(Tk)){
      stop("Nk and Tk should have the same length!")
   }
   em = mcmc + brn # total number of sampling
   ## compatibility check when predicting
   N_list = c()
   g_list = list()
   g.out = vector("list", em)
   ## g matrices for multiple N - before mixing
   for(k in 1:length(Nk)){
      N = Nk[k]
      g.out[[k]][[1]] = matrix(0, N+1, N+1)
   }
   ## generate ESS samples for the PT
   for(k in 1:length(Nk)){
      N = rep(Nk[k], 2)
      result = eigvals_exact(ndim = N, nu = nu.in, lambda_g = l.in, tausq)
      for(i in 2:em){
         nu.ess = matrix(samp_from_grid(N, mdim = result$mvec, egs = result$egvals, nu, lambda_g), 
                         nrow = N[1] + 1, ncol = N[2] + 1, byrow = TRUE) * sqrt(Tk[k])
         g.out[[k]][[i]] = ESS(g.out[[k]][[i-1]], nu.ess, z = Z, x = X, sigsq = sigsq, Temper = Tk[k])
         if (i%%100 == 0){
            print(c(i, N[1]))
         }
      }
   }
   for(i in 1:em){
      # mixing for PT
      for (kk in 1:1){
         ## swap states - Sambridge (2014)! randomly choose two chains and decide
         chains = sample(1:(length(Tk)), 2, replace = F)
         k1 = chains[1]
         k2 = chains[2]
         log.prob.swap = min(0, (Tk[k2] - Tk[k1])/(Tk[k1] * Tk[k2]) *
                                (loglik(Z, X, g.out[[k2]][[i]], sigsq) + log(N.pr(Nk[k2])) -
                                    loglik(Z, X, g.out[[k1]][[i]], sigsq) - log(N.pr(Nk[k1]))))
         u = runif(1)
         if (u < exp(log.prob.swap)){
            # swapping the sample - should change the remaining thing as a whole. 
            # Since we keep the required info of the previous steps, we swap the list as a hole.
            x = g.out[[k1]]
            g.out[[k1]] = g.out[[k2]]
            g.out[[k2]] = x
            # swapping the N
            y = Nk[k1]
            Nk[k1] = Nk[k2]
            Nk[k2] = y
         }
      }
      N_list[i] = Nk[1]
      g_list[[i]] = as.vector(t(g.out[[1]][[i]]))
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list))
}


sample.RJESS2D = function(Z, X, Nk, N.pr, mcmc, brn, thin, l.in = NULL, nu.in = NULL, sigsq, 
                          N.init, tausq){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   if(length(Z) != dim(X)[1])
      stop("Z and X should have the same length!")
   em = mcmc + brn # total number of sampling
   ## compatibility check when predicting
   N_list = c()
   g_list = list()
   g.out = vector("list", em)
   ## g matrices for multiple N - before mixing
   for(k in 1:length(Nk)){
      N = Nk[k]
      g.out[[k]][[1]] = matrix(0, N+1, N+1)
   }
   ## generate ESS samples for the PT
   for(k in 1:length(Nk)){
      N = rep(Nk[k], 2)
      result = eigvals_exact(ndim = N, nu = nu.in, lambda_g = l.in, tausq)
      for(i in 2:em){
         nu.ess = matrix(samp_from_grid(N, mdim = result$mvec, egs = result$egvals, nu, lambda_g), 
                         nrow = N[1] + 1, ncol = N[2] + 1, byrow = TRUE)
         g.out[[k]][[i]] = ESS(g.out[[k]][[i-1]], nu.ess, z = Z, x = X, sigsq = sigsq, Temper = 1)
         if (i%%100 == 0){
            print(c(i, N[1]))
         }
      }
   }
   for(i in 1:em){
      # mixing for PT
      for (kk in 1:1){
         ## swap states - Sambridge (2014)! randomly choose two chains and decide
         chains = sample(2:(length(Nk)), 1, replace = F)
         k = chains
         log.prob.swap = min(0, loglik(Z, X, g.out[[k]][[i]], sigsq) + log(N.pr(Nk[k])) -
                                    loglik(Z, X, g.out[[1]][[i]], sigsq) - log(N.pr(Nk[1])))
         u = runif(1)
         if (u < exp(log.prob.swap)){
            # swapping the sample - should change the remaining thing as a whole. 
            # Since we keep the required info of the previous steps, we swap the list as a hole.
            x = g.out[[1]]
            g.out[[1]] = g.out[[k]]
            g.out[[k]] = x
            # swapping the N
            y = Nk[1]
            Nk[1] = Nk[k]
            Nk[k] = y
         }
      }
      N_list[i] = Nk[1]
      g_list[[i]] = as.vector(t(g.out[[1]][[i]]))
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list))
}

# nu.in : BesselK(d, nu)
# l.in : K_nu(d/lambda) = K_nu(d/l.in)
# kappa = 1/l.in should hold
sample.RJESS2D.seq = function(Z, X, Nk, N.pr, kappak, kappa.pr, tausqk, tausq.pr, 
                              brn.ESS = 100, beta, mcmc, brn, sigsq, seed = 1234){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   if(length(Z) != dim(X)[1])
      stop("Z and X should have the same length!")
   n = nrow(X)
   em = mcmc + brn # total number of sampling
   g_list = list()
   N1 = length(Nk)
   N2 = length(kappak)
   N3 = length(tausqk)
   log_prob_N_list = vector(length = N1 * N2 * N3)
   N_list = kappa_list = tausq_list = vector(length = em)
   for(k1 in 1:N1){
      N = Nk[k1]
      Phi = Phi_2D(X, N)
      PhiTPhi = t(Phi) %*% Phi
      gridmat = cbind(rep(c(0:N)/N, each = N + 1), 
                      rep(c(0:N)/N, N + 1))
      for(k2 in 1:N2){
         for(k3 in 1:N3){
            index = (k1 - 1) * N2 * N3 + (k2 - 1) * N3 + k3 # index from 1 to N1 * N2 * N3
            kappa = kappak[k2]
            tausq = tausqk[k3]
            Sigma_N = thekernel(gridmat, nu = beta - 1, lambda = 1/kappa, tausq = tausq)
            Q_N = solve(Sigma_N, tol = 1e-30)
            Q_N_star = Q_N + PhiTPhi/sigsq
            mu_star = solve(Q_N_star, t(Phi) %*% Z) / sigsq
            log_prob_N_list[index] = log(N.pr(N)) + log(kappa.pr(kappa)) + log(tausq.pr(tausq)) -
               1/2 * log(det(diag((N+1)^2) + Sigma_N %*% PhiTPhi / sigsq)) +
               # a) 1/2 * log(det(diag((N+1)^2) + Sigma_N %*% PhiTPhi / sigsq)) +
               # b) 1/2 * (determinant(Q_N, logarithm = TRUE) - determinant(Q_N_star, logarithm = TRUE)) +
               1/2 * t(mu_star) %*% Q_N_star %*% mu_star - t(Z) %*% Z / 2 / sigsq
         }
      }
   }
   log_prob_N_list[is.na(log_prob_N_list)] = -Inf
   # sampling from p(N|D) - with fixed seed
   set.seed(seed)
   param_index_list = sample(1:(N1 * N2 * N3), size = em, replace = TRUE, prob = exp(log_prob_N_list - max(log_prob_N_list)))
   ## generate ESS samples
   for(param_index in 1:(N1 * N2 * N3)){
      index = which(param_index_list == param_index)
      if(length(index >= 1)){
         set.seed(seed * param_index)
         N = Nk[(param_index - 1) %/% (N2 * N3) + 1]
         kappa = kappak[((param_index - 1) %% (N2 * N3)) %/% N3 + 1]
         tausq = tausqk[(param_index - 1) %% N3 + 1]
         g.out = matrix(0, N+1, N+1)
         ##### procedure for the GPI sampling #####
         result = eigvals_exact(ndim = c(N, N), nu = beta - 1, lambda_g = 1/kappa, tausq = tausq)
         for(a in 1:(brn.ESS + length(index))){
            nu.ess = matrix(samp_from_grid(ndim = c(N, N), mdim = result$mvec, egs = result$egvals, nu = beta-1, lambda_g = 1/kappa, seed = seed * a * param_index), 
                            nrow = N + 1, ncol = N + 1, byrow = TRUE)
            g.out = ESS(g.out, nu.ess, z = Z, x = X, sigsq, Temper = 1, seed = seed * (param_index + 100))
            if(a > brn.ESS){
               g_list[[(index[a - brn.ESS])]] = t(g.out)
               print(a - brn.ESS)
            }
         }
         for(j in 1:length(index)){
            N_list[index[j]] = N
            kappa_list[index[j]] = kappa
            tausq_list[index[j]] = tausq
         }
      }
   }
   return(list(g_list = g_list, N_list = N_list, kappa_list = kappa_list,
               tausq_list = tausq_list, log_prob_N_list = log_prob_N_list))
}

# Sampling from continuous kappa incorporated 
# No more discretization
sample.GPI2D = function(Z, X, Nk, N.pr, kappa.pr, kappa.sampler, tausq.pr, tausq.sampler, 
                              brn.ESS = 100, beta, mcmc, brn, sigsq, seed = 1234){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   if(length(Z) != dim(X)[1])
      stop("Z and X should have the same length!")
   n = nrow(X)
   em = mcmc + brn # total number of sampling
   g_list = list()
   N_list = kappa_list = tausq_list = vector(length = em)
   # Sample initial N, kappa, tau
   N = sample(Nk, 1, prob = N.pr(Nk) / sum(N.pr(Nk)))
   kappa = kappa.sampler()
   tausq = tausq.sampler()
   param_now = param_post_2D(X, Z, beta = beta, N = N, kappa = kappa, tausq = tausq, 
                                 N.pr = N.pr, kappa.pr = kappa.pr, tausq.pr = tausq.pr, sigsq = sigsq)
   g.out = matrix(0, N+1, N+1)
   # result = eigvals_exact(ndim = c(N, N), nu = beta - 1, lambda_g = 1/kappa, tausq = tausq)
   for (i in 1:(mcmc+brn)){
      N.cand = sample(Nk, 1, prob = N.pr(Nk) / sum(N.pr(Nk)))
      kappa.cand = kappa.sampler()
      tausq.cand = tausq.sampler()
      
      param_cand = param_post_2D(X, Z, beta = beta, N = N.cand, kappa = kappa.cand, tausq = tausq.cand, 
                                     N.pr = N.pr, kappa.pr = kappa.pr, tausq.pr = tausq.pr, sigsq = sigsq)
      u = runif(1)
      if (u < exp(param_cand$log_prob[1] - param_now$log_prob[1])){
         N = N.cand
         kappa = kappa.cand
         tausq = tausq.cand
         param_now = param_cand
      }
      g.out = param_now$mu_star + backsolve(param_now$chol_Q_N_star, rnorm(length(param_now$mu_star)))
      g_list[[i]] = t(g.out)
      N_list[i] = N
      kappa_list[i] = kappa
      tausq_list[i] = tausq
      print(i)
   }
   
   return(list(g_list = g_list, N_list = N_list, kappa_list = kappa_list,
               tausq_list = tausq_list))
}

## log_prob_param function defined which is required for sample.GPI2D
param_post_2D = function(X, Z, beta, N, kappa, tausq, N.pr, kappa.pr, tausq.pr, sigsq){
   Phi = Phi_2D(X, N)
   PhiTPhi = t(Phi) %*% Phi
   gridmat = cbind(rep(c(0:N)/N, each = N + 1), 
                   rep(c(0:N)/N, N + 1))
   Sigma_N = thekernel(gridmat, nu = beta - 1, lambda = 1/kappa, tausq = tausq)
   Q_N = solve(Sigma_N, tol = 1e-30)
   Q_N_star = Q_N + PhiTPhi/sigsq
   chol_Q_N = chol(Q_N)
   chol_Q_N_star = chol(Q_N_star)
   logdet_diff = sum(2 * log(diag(chol_Q_N))) -  sum(2 * log(diag(chol_Q_N_star)))
   mu_star = solve(Q_N_star, t(Phi) %*% Z) / sigsq
   log_prob = log(N.pr(N)) + log(kappa.pr(kappa)) + log(tausq.pr(tausq)) -
      1/2 * logdet_diff +
      1/2 * t(mu_star) %*% Q_N_star %*% mu_star - t(Z) %*% Z / 2 / sigsq
   return(list(log_prob = log_prob, mu_star = mu_star, chol_Q_N_star = chol_Q_N_star))
}

param_check_2D = function(X, Z, beta, N_supp, kappa_supp, tausq_supp,  N.pr, kappa.pr, tausq.pr, sigsq, M = 10000){
   N1 = length(N_supp)
   N2 = length(kappa_supp)
   N3 = length(tausq_supp)
   N_list = rep(0, M)
   kappa_list = rep(0, M)
   tausq_list = rep(0, M)
   log_prob_N_list = vector(length = N1 * N2 * N3)
   for(k1 in 1:N1){
      for(k2 in 1:N2){
         for(k3 in 1:N3){
            index = (k1 - 1) * N2 * N3 + (k2 - 1) * N3 + k3 # index from 1 to N1 * N2 * N3
            N = N_supp[k1]
            kappa = kappa_supp[k2]
            tausq = tausq_supp[k3]
            log_prob_N_list[index] = param_post_2D(X, Z, beta = beta, N = N, kappa = kappa, tausq = tausq, 
                                                N.pr = N.pr, kappa.pr = kappa.pr, tausq.pr = tausq.pr, sigsq = sigsq)$log_prob
            }
      }
   }
   param_index_list = sample(1:(N1 * N2 * N3), size = M, replace = TRUE, prob = exp(log_prob_N_list - max(log_prob_N_list)))
   for(param_index in 1:(N1 * N2 * N3)){
      idx = which(param_index_list == param_index)
      N = N_supp[(param_index - 1) %/% (N2 * N3) + 1]
      kappa = kappa_supp[((param_index - 1) %% (N2 * N3)) %/% N3 + 1]
      tausq = tausq_supp[(param_index - 1) %% N3 + 1]
      N_list[idx] = N
      kappa_list[idx] = kappa
      tausq_list[idx] = tausq
   }
   return(list(N_list = N_list, kappa_list = kappa_list, tausq_list = tausq_list))
}

############### function splitting for one time & iterative time comparison ###################


sample.RJESS2D.onetime = function(Z, X, Nk, N.pr, kappak, kappa.pr, tausqk, tausq.pr, 
                                  brn.ESS = 100, beta, mcmc, brn, sigsq, seed = 1234){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   if(length(Z) != dim(X)[1])
      stop("Z and X should have the same length!")
   n = nrow(X)
   em = mcmc + brn # total number of sampling
   g_list = list()
   N1 = length(Nk)
   N2 = length(kappak)
   N3 = length(tausqk)
   log_prob_N_list = vector(length = N1 * N2 * N3)
   N_list = kappa_list = tausq_list = vector(length = em)
   result_list = vector("list", em)
   for(k1 in 1:N1){
      for(k2 in 1:N2){
         for(k3 in 1:N3){
            index = (k1 - 1) * N2 * N3 + (k2 - 1) * N3 + k3 # index from 1 to N1 * N2 * N3
            N = Nk[k1]
            kappa = kappak[k2]
            tausq = tausqk[k3]
            Phi = Phi_2D(X, N)
            PhiTPhi = t(Phi) %*% Phi
            gridmat = cbind(rep(c(0:N)/N, each = N + 1), 
                            rep(c(0:N)/N, N + 1))
            Sigma_N = thekernel(gridmat, nu = beta - 1, lambda = 1/kappa, tausq = tausq)
            Q_N = solve(Sigma_N, tol = 1e-30)
            Q_N_star = Q_N + PhiTPhi/sigsq
            mu_star = solve(Q_N_star, t(Phi) %*% Z) / sigsq
            log_prob_N_list[index] = log(N.pr(N)) + log(kappa.pr(kappa)) + log(tausq.pr(tausq)) +
               1/2 * log(det(diag((N+1)^2) + Sigma_N %*% PhiTPhi / sigsq)) +
               1/2 * t(mu_star) %*% Q_N_star %*% mu_star - t(Z) %*% Z / 2 / sigsq
         }
      }
   }
   # sampling from p(N|D) - with fixed seed
   set.seed(seed)
   param_index_list = sample(1:(N1 * N2 * N3), size = em, replace = TRUE, prob = exp(log_prob_N_list - max(log_prob_N_list)))
   ## generate ESS samples for the PT
   for(param_index in 1:(N1 * N2 * N3)){
      index = which(param_index_list == param_index)
      if(length(index >= 1)){
         N = Nk[(param_index - 1) %/% (N2 * N3) + 1]
         kappa = kappak[((param_index - 1) %% (N2 * N3)) %/% N3 + 1]
         tausq = tausqk[(param_index - 1) %% N3 + 1]
         result = eigvals_exact(ndim = c(N, N), nu = beta - 1, lambda_g = 1/kappa, tausq = tausq)
         result_list[[param_index]] = result
      }
   }
   return(list(result_list = result_list, param_index_list = param_index_list))
}

# once one time calculation is done
sample.RJESS2D.iter = function(Z, X, Nk, N.pr, kappak, kappa.pr, tausqk, tausq.pr, result_list, param_index_list, brn.ESS, sigsq, seed = 1234){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   ## generate ESS samples for the PT
   em = length(param_index_list)
   g_list = list()
   N1 = length(Nk)
   N2 = length(kappak)
   N3 = length(tausqk)
   log_prob_N_list = vector(length = N1 * N2 * N3)
   N_list = kappa_list = tausq_list = vector(length = em)
   for(param_index in 1:(N1 * N2 * N3)){
      index = which(param_index_list == param_index)
      if(length(index >= 1)){
         N = Nk[(param_index - 1) %/% (N2 * N3) + 1]
         kappa = kappak[((param_index - 1) %% (N2 * N3)) %/% N3 + 1]
         tausq = tausqk[(param_index - 1) %% N3 + 1]
         result = result_list[[param_index]]
         g.out = matrix(0, N+1, N+1)
         for(a in 1:(brn.ESS + length(index))){
            nu.ess = matrix(samp_from_grid(ndim = c(N, N), mdim = result$mvec, egs = result$egvals, nu = beta-1, lambda_g = 1/kappa, seed = seed * a * param_index), 
                            nrow = N + 1, ncol = N + 1, byrow = TRUE)
            g.out = ESS(g.out, nu.ess, z = Z, x = X, sigsq, Temper = 1, seed = seed * (param_index + 100))
            if(a > brn.ESS){
               g_list[[(index[a - brn.ESS])]] = t(g.out)
            }
         }
         for(j in 1:length(index)){
            N_list[index[j]] = N
            kappa_list[index[j]] = kappa
            tausq_list[index[j]] = tausq
         }
      }
   }
   return(list(g_list = g_list, N_list = N_list, kappa_list = kappa_list))
}
