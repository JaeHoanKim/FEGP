sample.ESS.Nfixed2D = function(X, Z, kappa.pr = function(x){return(1)}, 
                               beta = 2, mcmc, brn, thin, sigsq = 0.01, kappa.init = 2, N.init = 20, 
                               g.init = rep(0, N.init+1)){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa init should be 2 to match with GPI method
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(dim(X)[1]!=length(Z))
      stop("X and Z should be of same length!")
   n = length(Z)
   em = mcmc + brn # total number of sampling
   N_list = c()
   prob_list = c()
   g_list = list()
   kappa_list = list()
   # starting with the initial N and g
   # setting up for iterative work
   # g.in = g.init
   N = N.init
   kappa = kappa.init
   mv = 1/(4*pi*kappa^2)
   g.in = Sampling_N_new2D(N, kappa) / sqrt(mv)
   for(i in 1:em){
      g_ESS = Sampling_N_new2D(N, kappa) / sqrt(mv)
      # Here, byrow should be F to be aligned with f_N_2D_multi ordering
      g.out = ESS_post2D(Z, X, matrix(g.in, N+1, N+1, byrow = F), matrix(g_ESS, N+1, N+1, byrow = F), sigsq)
      N.out = N
      kappa.out = kappa
      ### track g.out and N.out 
      g_list[[i]] = g.out
      N_list[i] = N.out
      kappa_list[[i]] = kappa.out
      ### preparation for iteration
      g.in = g.out
      N = N.out
      kappa = kappa.out
      if (i%%100 == 0){
         print(c("interation number:", i))
         print(c("N: ", N))
      }
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list, prob_list = prob_list, kappa_list = kappa_list))
}

## goal: draw samples of (kappa, N, w) given the data D_n (X_n, Y_n) with updating N
## data: X (n by 2 matrix, given by x, y coordinate of all observations), Z: n by 1 vector


sample.ESS.2D = function(X, Z, kappa.pr = function(x){return(1)},N.pr = function(x){return(exp(-5) * 5^x / factorial(x))},
                         beta = 2, mcmc, brn, thin, sigsq = 0.01, kappa.init = 100, N.init = 20, 
                         g.init = rep(0, N.init+1), gridsize = 15){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(dim(X)[1]!=length(Z))
      stop("X and Z should be of same length!")
   n = length(Z)
   em = mcmc + brn # total number of sampling
   N_list = c()
   prob_list = c()
   g_list = list()
   kappa_list = list()
   # starting with the initial N and g
   # setting up for iterative work
   # g.in = g.init
   N = N.init
   kappa = kappa.init
   mv2 = 1/(4*pi*kappa^2)
   mv1 = 1/(4*kappa^3)
   g.in = as.vector(Sampling_N_new2D(N, kappa)) / sqrt(mv2)  # (N+1)^2 vector
   for(i in 1:em){
      ### which stage to run 
      ## if u < -2 : consider N/2; else if u < -1: consider staying; else, consider double N
      uu = ((N==2) - 3) * runif(1)
      u = runif(1) # random uniform r.v. for the jump
      g.in.mat = matrix(g.in, N+1, N+1, byrow = F)
      if(uu > -1){
         # N increase
         g.cand.mat = matrix(0, N+2, N+2)
         g.cand.mat[1:(N+1), 1:(N+1)] = g.in.mat
         g.cand.mat[N+2, 1:(N+1)] = as.vector(g.in.mat[N+1, 1:(N+1)] + Sampling_N_new1D(N, kappa) / sqrt(mv1))
         g.cand.mat[1:(N+1), N+2] = as.vector(g.in.mat[1:(N+1), N+1] + Sampling_N_new1D(N, kappa) / sqrt(mv1))
         g.cand.mat[N+2, N+2] = 0.5 * (g.cand.mat[N+1, N+2] + g.cand.mat[N+2, N+1]) +
            rnorm(1, sd = 0.5*kappa^(-1.5))
         g.cand = matrix(g.cand.mat, ncol = 1, byrow = F) # as (N+2)^2 by 1 vector
         log.test.prob = min(0, log(N.pr(N+1) / N.pr(N)) +
                                loglik2D(Z, X, g.cand.mat, sigsq) - loglik2D(Z, X, g.in.mat, sigsq) +
                                loglik_w2D(N+1, kappa, g.cand) - loglik_w2D(N, kappa, g.in) -
                                loglik_w1D(N, kappa, g.cand.mat[N+2, 1:(N+1)] - g.in.mat[N+1, 1:(N+1)]) -
                                loglik_w1D(N, kappa, g.cand.mat[1:(N+1), N+2] - g.in.mat[1:(N+1), N+1]) - 
                                log(dnorm(g.cand.mat[N+2, N+2] - 0.5*g.cand.mat[N+1, N+2] - 0.5*g.cand.mat[N+2, N+1],
                                          sd = 0.5*kappa^(-1.5))))
         if (exp(log.test.prob) > u){
            g.out = g.cand
            N.out = N+1
         } else{
            g.out = g.in
            N.out = N
         }
      } else if(uu > -2){
         v1.mat = matrix(Sampling_N_new2D(N, kappa), N+1, N+1, byrow = F) / sqrt(mv1)
         g.out = matrix(ESS_post2D(Z, X, g.in.mat, v1.mat, sigsq), ncol = 1, byrow = F)
         N.out = N
         log.test.prob = 2
      } else if (uu > -3){
         # N decrease
         g.cand.mat = g.in.mat[1:N, 1:N]
         g.cand = matrix(g.cand.mat, ncol = 1) # as N^2 by 1 vector
         log.test.prob = min(0, log(N.pr(N-1) / N.pr(N)) +
                                loglik2D(Z, X, g.cand.mat, sigsq) - loglik2D(Z, X, g.in.mat, sigsq) +
                                loglik_w2D(N-1, kappa, g.cand) - loglik_w2D(N, kappa, g.in) +
                                loglik_w1D(N-1, kappa, g.in.mat[N+1, 1:N] - g.in.mat[N, 1:N]) +
                                loglik_w1D(N-1, kappa, g.in.mat[1:N, N+1] - g.in.mat[1:N, N]) +
                                log(dnorm(g.in.mat[N+1, N+1] - 0.5*g.in.mat[N, N+1] - 0.5*g.in.mat[N+1, N],
                                          sd = 0.5*kappa^(-1.5))))
         u = runif(1)
         if (exp(log.test.prob) > u){
            g.out = g.cand
            N.out = N-1
         }
         else{
            g.out = g.in
            N.out = N
         }
      }
      ### track g.out and N.out 
      g_list[[i]] = g.out
      N_list[i] = N.out
      prob_list[i] = log.test.prob
      ### preparation for iteration
      g.in = g.out
      N = N.out
      # kappa = kappa.out
      if (i%%100 == 0){
         print(c("interation number:", i))
         print(c("N: ", N))
      }
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list, prob_list = prob_list))
}


sample.exact2D= function(X, Z, kappa.pr = function(x){return(1)}, 
                               beta = 2, mcmc, brn, thin, sigsq = 0.01, kappa.init = 2, 
                               gridsize){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa init should be 2 to match with GPI method
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(dim(X)[1]!=length(Z))
      stop("X and Z should be of same length!")
   n = length(Z)
   em = mcmc + brn # total number of sampling
   N_list = c()
   prob_list = c()
   g_list = list()
   kappa_list = list()
   # starting with the initial N and g
   # setting up for iterative work
   # g.in = g.init
   kappa = kappa.init
   gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                   rep(c(0:gridsize)/gridsize, gridsize + 1))
   N = gridsize
   # sparse matrix Omega and Phi
   Omega = Q2D(N, kappa)
   Phi = Phi_2D(X, N)
   # computation of the mean and the variance vector
   var_grid = solve(Omega + t(Phi) %*% Phi / sigsq)
   mean_grid = var_grid %*% t(Phi) %*% Z / sigsq
   # symmetrize due to prevent the numerical error
   var_grid = (var_grid + t(var_grid)) / 2
   g_samples = mvtnorm::rmvnorm(n = em, mean = mean_grid, sigma = var_grid,
                                checkSymmetry = FALSE)
   for(i in 1:em){
      g_list[[i]] = g_samples[i, ] 
         # as.vector(t(matrix(g_samples[i, ], nrow = sqrt(length(g_samples[i, ])), byrow = TRUE)))
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list))
}

sample.RJexact = function(X, Z, kappa.pr = function(x){return(1)}, Nk, N.pr,
                          beta = 2, mcmc, brn, thin, sigsq = 0.01, kappa.init = 2, N.init = 20, 
                          g.init = rep(0, N.init+1)){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa init should be 2 to match with GPI method
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(dim(X)[1]!=length(Z))
      stop("X and Z should be of same length!")
   n = length(Z)
   em = mcmc + brn # total number of sampling
   N_list = c()
   prob_list = c()
   g_list = vector("list", em)
   kappa_list = list()
   # g: double-indexed list, first: about chain, second: about i-th sample
   # g.in = list()
   g.out = vector("list", em)
   g.out.mat = vector("list", em)
   # starting with the initial N and g
   # setting up for iterative work
   # g.in = g.init
   kappa = kappa.init
   mv = 1/(4*pi*kappa^2)
   log_jump_prob_N = vector(length = length(Nk))
   # 1. Progress within the chain
   for(k in 1:length(Nk)){
      N = Nk[k]
      kappa = kappa.init
      gridsize = N
      gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                      rep(c(0:gridsize)/gridsize, gridsize + 1))
      # sparse matrix Omega and Phi
      Omega = Q2D(N, kappa)
      Phi = Phi_2D(X, N)
      # computation of the mean and the variance vector
      var_grid = solve(Omega + t(Phi) %*% Phi / sigsq)
      mean_grid = var_grid %*% t(Phi) %*% Z / sigsq
      # symmetrize due to prevent the numerical error
      var_grid = (var_grid + t(var_grid)) / 2
      g_samples = mvtnorm::rmvnorm(n = em, mean = mean_grid, sigma = var_grid,
                                   checkSymmetry = FALSE)
      log_jump_prob_N[k] = log(N.pr(Nk[k])) - 1/2 * log(det(diag((N+1)^2) + solve(Omega) %*% t(Phi) %*% Phi))
      for(i in 1:em){
         g.out[[k]][[i]] = # as.vector(t(matrix(g_samples[i, ], nrow = sqrt(length(g_samples[i, ])), byrow = FALSE)))
            g_samples[i, ]
         g.out.mat[[k]][[i]] = matrix(g.out[[k]][[i]], N+1, N+1, byrow = T)
      }
   }
   for(i in 1:em){
      # 2. Swapping between the chain
      # For the RJMCMC, randomly choose another chain other than the first one and compare
      for (kk in 1:1){
         ## swap states - Sambridge (2014)! randomly choose two chains and decide
         chains = sample(2:(length(Nk)), 1, replace = F)
         k = chains
         # jump probability - not from Brooks et al., but from the direct calculation of the detailed balance condition
         log.jump.prob = min(0, log_jump_prob_N[k]- log_jump_prob_N[1])
         u = runif(1)
         if (u < exp(log.jump.prob)){
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
      ### preparation for iteration
      g.in = g.out
      if (i %% 100 == 0){
         print(c("iteration number: ", i))
         print(c("N: ", sqrt(length(g.out[[1]][[i]])) - 1))
         print(c("N mixing status :", Nk))
      }
      g_list[[i]] = g.in[[1]][[i]]
      N_list[i] = Nk[1]
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list, kappa_list = kappa_list, log_jump_prob_N = log_jump_prob_N))
}

sample.exact2D.seq = function(X, Z, Nk, N.pr, kappak, kappa.pr, tausqk, tausq.pr,
                          beta = 2, mcmc, brn, sigsq = 0.01, seed = 1234){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa init should be 2 to match with GPI method
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(dim(X)[1]!=length(Z))
      stop("X and Z should be of same length!")
   n = length(Z)
   em = mcmc + brn # total number of sampling
   N_list = c()
   prob_list = c()
   g_list = vector("list", em)
   N1 = length(Nk)
   N2 = length(kappak)
   N3 = length(tausqk)
   log_prob_N_list = vector(length = N1 * N2 * N3)
   N_list = kappa_list = tausq_list = vector(length = em)
   mean_grid = list()
   prec_grid = list()
   chol_prec_grid = list()
   nu = beta - 1
   tausq0 = gamma(nu) / gamma(nu + 1) / (4*pi) / kappa^(2*nu)
   # 1. Sample from N | D
   for(k1 in 1:N1){
      for (k2 in 1:N2){
         for(k3 in 1:N3){
            index = (k1 - 1) * N2 * N3 + (k2 - 1) * N3 + k3
            N = Nk[k1]
            kappa = kappak[k2]
            tausq = tausqk[k3]
            Omega = Q2D(N, kappa, tausq, beta = 2)
            Phi = Phi_2D(X, N)
            # computation of the mean and the variance vector
            prec_grid[[index]] = Omega + t(Phi) %*% Phi / sigsq
            chol_prec_grid[[index]] = chol(prec_grid[[index]])
            mean_grid[[index]] = solve(prec_grid[[index]], t(Phi) %*% Z / sigsq)
            log_prob_N_list[index] = log(N.pr(N)) + log(kappa.pr(kappa)) + log(tausq.pr(tausq)) - 
               1/2 * log(det(prec_grid[[index]])) + 1/2 * log(det(Omega)) +
               1/2 * t(mean_grid[[index]]) %*% prec_grid[[index]] %*% mean_grid[[index]] - t(Z) %*% Z/(2*sigsq)
         }
      }
   }
   # sample from p(N|D)
   set.seed(seed)
   param_index_list = sample(1:(N1 * N2 * N3), size = em, replace = TRUE, prob = exp(log_prob_N_list - max(log_prob_N_list)))
   for(param_index in 1:(N1 * N2 * N3)){
      index = which(param_index_list == param_index)
      if(length(index >= 1)){
         set.seed(seed * param_index)
         stdnorms = matrix(rnorm(length(index) * length(mean_grid[[param_index]])), nrow = length(mean_grid[[param_index]]))
         g_samples = Matrix::solve(t(chol_prec_grid[[param_index]]), stdnorms) + mean_grid[[param_index]]
         g_samples = t(g_samples)
         for(j in 1:length(index)){
            g_list[[(index[j])]] = g_samples[j, ]
            N_list[index[j]] = Nk[(param_index - 1) %/% (N2 * N3) + 1]
            kappa_list[index[j]] = kappak[((param_index - 1) %% (N2 * N3)) %/% N3 + 1]
            tausq_list[index[j]] = tausqk[(param_index - 1) %% N3 + 1]
         }
      }
   }
   return(list(g_list = g_list, N_list = N_list, kappa_list = kappa_list, log_prob_N_list = log_prob_N_list))
}

sample.exact.onetime = function(X, Z, Nk, N.pr, kappak, kappa.pr, tausqk, tausq.pr,
                                beta = 2, mcmc, brn, sigsq = 0.01, seed = 1234){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa init should be 2 to match with GPI method
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(dim(X)[1]!=length(Z))
      stop("X and Z should be of same length!")
   n = length(Z)
   em = mcmc + brn # total number of sampling
   N_list = c()
   prob_list = c()
   g_list = vector("list", em)
   N1 = length(Nk)
   N2 = length(kappak)
   N3 = length(tausqk)
   log_prob_N_list = vector(length = N1 * N2 * N3)
   N_list = kappa_list = tausq_list = vector(length = em)
   mean_grid = list()
   prec_grid = list()
   chol_prec_grid = list()
   nu = beta - 1
   tausq0 = gamma(nu) / gamma(nu + 1) / (4*pi) / kappa^(2*nu)
   # 1. Sample from N | D
   for(k1 in 1:N1){
      for (k2 in 1:N2){
         for(k3 in 1:N3){
            index = (k1 - 1) * N2 * N3 + (k2 - 1) * N3 + k3
            N = Nk[k1]
            kappa = kappak[k2]
            tausq = tausqk[k3]
            Omega = Q2D(N, kappa, tausq, beta = 2)
            Phi = Phi_2D(X, N)
            # computation of the mean and the variance vector
            prec_grid[[index]] = Omega + t(Phi) %*% Phi / sigsq
            chol_prec_grid[[index]] = chol(prec_grid[[index]])
            mean_grid[[index]] = solve(prec_grid[[index]], t(Phi) %*% Z / sigsq)
            log_prob_N_list[index] = log(N.pr(N)) + log(kappa.pr(kappa)) + log(tausq.pr(tausq)) - 
               1/2 * log(det(prec_grid[[index]])) + 1/2 * log(det(Omega)) +
               1/2 * t(mean_grid[[index]]) %*% prec_grid[[index]] %*% mean_grid[[index]] - t(Z) %*% Z/(2*sigsq)
         }
      }
   }
   # sample from p(N|D)
   set.seed(seed)
   param_index_list = sample(1:(N1 * N2 * N3), size = em, replace = TRUE, prob = exp(log_prob_N_list - max(log_prob_N_list)))
   return(param_index_list = param_index_list, mean_grid = mean_grid, chol_prec_grid = chol_prec_grid)
}


sample.exact.iter = function(Nk, N.pr, kappak, kappa.pr, tausqk, tausq.pr,
                             beta = 2, mcmc, brn, sigsq = 0.01, param_index_list, 
                             chol_prec_grid, mean_grid, seed = 1234){
   g_list = vector("list", em)
   N1 = length(Nk)
   N2 = length(kappak)
   N3 = length(tausqk)
   log_prob_N_list = vector(length = N1 * N2 * N3)
   N_list = kappa_list = tausq_list = vector(length = em)
   for(param_index in 1:(N1 * N2 * N3)){
      index = which(param_index_list == param_index)
      if(length(index >= 1)){
         set.seed(seed * param_index)
         stdnorms = matrix(rnorm(length(index) * length(mean_grid[[param_index]])), nrow = length(mean_grid[[param_index]]))
         g_samples = Matrix::solve(t(chol_prec_grid[[param_index]]), stdnorms) + mean_grid[[param_index]]
         g_samples = t(g_samples)
         for(j in 1:length(index)){
            g_list[[(index[j])]] = g_samples[j, ]
            N_list[index[j]] = Nk[(param_index - 1) %/% (N2 * N3) + 1]
            kappa_list[index[j]] = kappak[((param_index - 1) %% (N2 * N3)) %/% N3 + 1]
            tausq_list[index[j]] = tausqk[(param_index - 1) %% N3 + 1]
         }
      }
   }
   return(list(g_list = g_list, N_list = N_list, kappa_list = kappa_list, log_prob_N_list = log_prob_N_list))
}
