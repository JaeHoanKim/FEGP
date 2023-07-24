# sampling g_N with N.fixed

sample.ESS.Nfixed = function(X, Y, kappa.pr = function(x){return(1)}, 
                             beta = 2, mcmc, brn, thin, sigsq, kappa.init, N.init, 
                             Temp = 1){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(length(X)!=length(Y))
      stop("X and Y should be of same length!")
   n = length(Y)
   em = mcmc + brn # total number of sampling
   N_list = c()
   g_list = list()
   N = N.init
   kappa = kappa.init
   mv = 1/(4*kappa^3)
   g.in = Sampling_N_new(N, kappa)/sqrt(mv) # when I bring this line inside the for loop, it does not seem to converge.
   for(i in 1:em){
      g_ESS = Sampling_N_new(N, kappa)/sqrt(mv)
      g.out = ESS_post_tempered(Y, X, g.in, g_ESS, sigsq, Temp)
      N.out = N
      ### track g.out and N.out ###
      g_list[[i]] = g.out
      N_list[i] = N.out
      ### preparation for iteration ###
      g.in = g.out
      N = N.out
      if (i%%1000 == 0){
         print(c("interation number:", i))
         print(c("N: ", N))
      }
   }
   return(list(g_list = g_list, N_list = N_list))
}

#### Using Parallel tempering
# Tk : not a function, but a vector
sample.PT.ESS = function(X, Y, kappa.pr = function(x){return(1)}, Nk,
                         Tk = c(1:length(Nk)), N.pr, sigsq, kappa.init,
                         mcmc, brn, thin){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   # temperature length check
   # Sort Nk in an increasing order
   Nk = sort(Nk)
   if(length(Nk) != length(Tk)){
      stop("The number of temperatures should be equal to the number of N s!")
   }
   if (Tk[1] != 1){
      stop("The temperature of the first element of Tk should be equal to 1!")
   }
   if(length(X)!=length(Y))
      stop("X and Y should be of same length!")
   n = length(Y)
   em = mcmc + brn # total number of sampling
   ## initialize the collection of g at first step
   g.in = list()
   g.out = list()
   g_list = list()
   N_list = c()
   kappa = kappa.init
   mv = 1/(4*kappa^3)
   # Nk: the N value (nrow(g) - 1) of the k-th chain!
   for(k in 1:length(Nk)){
      g.in[[k]] = Sampling_N_new(Nk[k], kappa)/sqrt(mv) * sqrt(Tk[k])
   }
   ## Parallel tempering
   for (i in 1:em){
      ## sample of g at i th step
      # k : chain order
      # g.out : sample from N(g;0, Sigma) L(g)^{1/T}
      for(k in 1:length(Nk)){
         # reflecting the changing dimension of g.in[[k]]
         g_ESS = Sampling_N_new(Nk[k], kappa)/sqrt(mv) * sqrt(Tk[k])
         g.out[[k]] = ESS_post_tempered(Y, X, g.in[[k]], g_ESS, sigsq, Temp = Tk[k])
      }
      # try multiple swaps after one within sampling
      # regulate the number of mixing
      for (kk in 1:(1*length(Nk))){
         ## swap states - Sambridge (2014)! randomly choose two chains and decide
         chains = sample(1:(length(Tk)), 2, replace = F)
         k1 = chains[1]
         k2 = chains[2]
         log.prob.swap = min(0, (Tk[k2] - Tk[k1])/(Tk[k1] * Tk[k2]) *
                                (loglik(Y, X, g.out[[k2]], sigsq) + log(N.pr(Nk[k2])) -
                                    loglik(Y, X, g.out[[k1]], sigsq) - log(N.pr(Nk[k1]))))
         u = runif(1)
         if (u < exp(log.prob.swap)){
            # swapping the sample
            x = g.out[[k1]]
            g.out[[k1]] = g.out[[k2]]
            g.out[[k2]] = x
            # swapping the N
            y = Nk[k1]
            Nk[k1] = Nk[k2]
            Nk[k2] = y
         }
      }
      g.in = g.out
      if (i %% 100 == 0){
         print(c("iteration number: ", i))
         print(c("N: ", nrow(g.out[[1]]) - 1))
         print(c("N mixing status :", Nk))
      }
      g_list[[i]] = g.in[[1]]
      N_list[i] = nrow(g.in[[1]]) - 1
   }
   return(list(g_list = g_list, N_list = N_list))
}

## Sample exact using the original GP

sample.exact = function(X, Y, kappa.pr = function(x){return(1)}, 
                        beta = 2, mcmc, brn, thin, sigsq, kappa.init, gridsize){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(length(X)!=length(Y))
      stop("X and Y should be of same length!")
   n = length(Y)
   if(missing(gridsize)){
      stop("grid should be provided for the exact Nfixed method!")
   }
   em = mcmc + brn # total number of sampling
   g_list = list()
   N = gridsize
   # starting with the initial N and g
   # setting up for iterative work
   # g.in = g.init
   kappa = kappa.init
   Omega = Q1D(N, kappa)
   Phi = Phi_1D(X, N)
   # computation of the mean and the variance vector
   var_grid = solve(Omega + t(Phi) %*% Phi / sigsq)
   mean_grid = var_grid %*% t(Phi) %*% Y / sigsq
   g_samples = rmvnorm(n = em, mean = mean_grid, sigma = var_grid, checkSymmetry = FALSE)
   for(i in 1:em){
      g_list[[i]] = g_samples[i, ]
   }
   return(list(g_list = g_list))
}


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


