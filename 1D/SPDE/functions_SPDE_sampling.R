## Sample exact using the original GP

sample.exact.seq = function(X, Y, kappa.pr = function(x){return(1)}, Nk, N.pr,
                        beta, mcmc, brn, sigsq = 0.01, kappa.init, seed = 1000){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(length(X)!=length(Y))
      stop("X and Y should be of same length!")
   n = length(Y)
   em = mcmc + brn # total number of sampling
   g_list = list()
   log_jump_prob_N = vector(length = length(Nk))
   var_grid = list()
   mean_grid = list()
   prec_grid = list()
   for(k in 1:length(Nk)){
      N = Nk[k]
      Omega = Q1D(N, kappa, beta)
      Phi = Phi_1D(X, N)
      # computation of the mean and the variance vector
      prec_grid = Omega + t(Phi) %*% Phi / sigsq
      var_grid[[k]] = solve(Omega + t(Phi) %*% Phi / sigsq)
      mean_grid[[k]] = var_grid[[k]] %*% t(Phi) %*% Y / sigsq
      # computation of p(N | D)
      log_jump_prob_N[k] = log(N.pr(Nk[k])) - 1/2 * log(det(diag((N+1)) + solve(Omega) %*% t(Phi) %*% Phi / sigsq)) +
         1/2 * t(mean_grid[[k]]) %*% (Omega + t(Phi) %*% Phi / sigsq) %*% mean_grid[[k]] - t(Y) %*% Y/(2*sigsq)
   }
   # sample from p(N|D)
   set.seed(seed)
   N_list = sample(Nk, size = em, replace = TRUE, prob = exp(log_jump_prob_N - log_jump_prob_N[length(Nk)]))
   
   for(k in 1:length(Nk)){
      N = Nk[k]
      index = which(N_list == N)
      if (length(index) >= 1){
         set.seed(seed * k)
         g_samples <- mvtnorm::rmvnorm(n = length(index), mean = mean_grid[[k]], sigma = var_grid[[k]],
                                      checkSymmetry = FALSE)
         for(j in 1:length(index)){
            g_list[[(index[j])]] = g_samples[j, ]
         }
      }
   }
   return(list(g_list = g_list, N_list = N_list, log_jump_prob_N = log_jump_prob_N))
}


