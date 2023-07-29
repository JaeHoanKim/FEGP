## Sample exact using the original GP

sample.exact.seq = function(X, Y, kappa.pr = function(x){return(1)}, Nk, N.pr, kappak, kappa.pr, tausqk, tausq.pr,
                        beta, mcmc, brn, sigsq = 0.01, kappa.init, seed = 1000){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(length(X)!=length(Y))
      stop("X and Y should be of same length!")
   n = length(Y)
   em = mcmc + brn # total number of sampling
   g_list = list()
   log_prob_N_list = vector(length = length(Nk * kappak * tausqk))
   var_grid = list()
   mean_grid = list()
   prec_grid = list()
   N1 = length(Nk)
   N2 = length(kappak)
   N3 = length(tausqk)
   for(k1 in 1:N1){
      for(k2 in 1:N2){
         for(k3 in 1:N3){
            index = (k1 - 1) * N2 * N3 + (k2 - 1) * N3 + k3 # index from 1 to N1 * N2 * N3
            N = Nk[k1]
            kappa = kappak[k2]
            tausq = tausqk[k3]
            Omega = Q1D(N, kappa, beta, tausq)
            Phi = Phi_1D(X, N)
            # computation of the mean and the variance vector
            prec_grid = Omega + t(Phi) %*% Phi / sigsq
            var_grid[[index]] = solve(Omega + t(Phi) %*% Phi / sigsq)
            mean_grid[[index]] = var_grid[[index]] %*% t(Phi) %*% Y / sigsq
            # computation of p(N | D)
            log_prob_N_list[index] = log(N.pr(N)) + log(kappa.pr(kappa)) + log(tausq.pr(tausq)) - 
               1/2 * log(det(diag((N+1)) + solve(Omega) %*% t(Phi) %*% Phi / sigsq)) +
               1/2 * t(mean_grid[[index]]) %*% (Omega + t(Phi) %*% Phi / sigsq) %*% mean_grid[[index]] - t(Y) %*% Y/(2*sigsq)
         }
      }
   }
   # sample from p(N|D)
   set.seed(seed)
   N_list = sample(Nk, size = em, replace = TRUE, prob = exp(log_prob_N_list - log_prob_N_list[length(Nk)]))
   
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
   return(list(g_list = g_list, N_list = N_list, log_prob_N_list = log_prob_N_list))
}


