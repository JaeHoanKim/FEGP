## Sample exact using the original GP

sample.exact.seq = function(X, Y, Nk, N.pr, kappak, kappa.pr, tausqk, tausq.pr,
                        beta, mcmc, brn, sigsq = 0.01, kappa.init, seed = 1000){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(length(X)!=length(Y))
      stop("X and Y should be of same length!")
   n = length(Y)
   em = mcmc + brn # total number of sampling
   g_list = list()
   N1 = length(Nk)
   N2 = length(kappak)
   N3 = length(tausqk)
   log_prob_N_list = vector(length = N1 * N2 * N3)
   N_list = kappa_list = tausq_list = vector(length = em)
   var_grid = list()
   mean_grid = list()
   prec_grid = list()
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
               1/2 * log(det(prec_grid)) + 1/2 * log(det(Omega)) +
               1/2 * t(mean_grid[[index]]) %*% prec_grid %*% mean_grid[[index]] - t(Y) %*% Y/(2*sigsq)
         }
      }
   }
   # sample from p(N|D)
   set.seed(seed)
   param_index_list = sample(1:(N1 * N2 * N3), size = em, replace = TRUE, prob = exp(log_prob_N_list - log_prob_N_list[N1 * N2 * N3]))
   for(param_index in 1:(N1 * N2 * N3)){
      index = which(param_index_list == param_index)
      if(length(index >= 1)){
         set.seed(seed * param_index)
         g_samples <- mvtnorm::rmvnorm(n = length(index), mean = mean_grid[[param_index]], sigma = var_grid[[param_index]],
                                       checkSymmetry = FALSE)
         for(j in 1:length(index)){
            g_list[[(index[j])]] = g_samples[j, ]
            N_list[index[j]] = Nk[(param_index - 1) %/% (N2 * N3) + 1]
            kappa_list[index[j]] = kappak[((param_index - 1) %% (N2 * N3)) %/% N3 + 1]
            tausq_list[index[j]] = tausqk[(param_index - 1) %% N3 + 1]
         }
      }
   }
   return(list(g_list = g_list, N_list = N_list, kappa_list = kappa_list,
               tausq_list = tausq_list, log_prob_N_list = log_prob_N_list))
}
