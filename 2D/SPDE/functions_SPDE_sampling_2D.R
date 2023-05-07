
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
                         g.init = rep(0, N.init+1), pred = TRUE, Xtest = NULL, gridsize = 15){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## kappa.pr: prior distribution of kappa (function); by default, noninformative prior
   if(dim(X)[1]!=length(Z))
      stop("X and Z should be of same length!")
   n = length(Z)
   if (length(g.init) != N.init + 1){
      stop("provided g.init should have the length of N.init!")
   }
   if(missing(mcmc))
      mcmc=500
   if(missing(brn))
      brn=100
   if(missing(thin))
      thin=1
   em = mcmc + brn # total number of sampling
   ## compatibility check when predicting
   if(pred == TRUE){
      if(is.null(Xtest) == TRUE){
         stop("Xtest should be provided when performing prediction!")
      }
   }
   N_list = c()
   prob_list = c()
   g_list = list()
   kappa_list = list()
   ntest = length(Xtest)
   # starting with the initial N and g
   # setting up for iterative work
   # g.in = g.init
   N = N.init
   kappa = kappa.init
   g.in = as.vector(Sampling_N_new2D(N, kappa))  # (N+1)^2 vector
   for(i in 1:em){
      jump.prob = runif(1)
      g.in.mat = matrix(g.in, N+1, N+1, byrow = F)
      if ((N>=3 & jump.prob < 1/3) | (N==2 & jump.prob < 1/2)){
         # N increase
         g.cand.mat = matrix(0, N+2, N+2)
         g.cand.mat[1:(N+1), 1:(N+1)] = g.in.mat
         g.cand.mat[N+2, 1:(N+1)] = as.vector(g.in.mat[N+1, 1:(N+1)] + Sampling_N_new1D(N, kappa))
         g.cand.mat[1:(N+1), N+2] = as.vector(g.in.mat[1:(N+1), N+1] + Sampling_N_new1D(N, kappa))
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
         u = runif(1)
         if (exp(log.test.prob) > u){
            g.out = g.cand
            N.out = N+1
         } else{
            g.out = g.in
            N.out = N
         }
      } else if (N >=3 & jump.prob < 2/3){
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
      } else{
         # N stay
         v1.mat = matrix(Sampling_N_new2D(N, kappa), N+1, N+1, byrow = F)
         g.out = matrix(ESS_post2D(Z, X, g.in.mat, v1.mat, sigsq), ncol = 1, byrow = F)
         N.out = N
         log.test.prob = 2
      }
      ### track g.out and N.out 
      g_list[[i]] = g.out
      N_list[i] = N.out
      prob_list[i] = log.test.prob
      gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                      rep(c(0:gridsize)/gridsize, gridsize+ 1))
      ### preparation for iteration
      g.in = g.out
      N = N.out
      ## prediction
      if (i == brn+1){
         Z.fitted = matrix(0, (gridsize + 1)^2, 1) # (gridsize+1)^2 by 1 matrix
      }
      if(i > brn){
         Zhat = f_N_h_2D_multi(gridmat, matrix(g.out, N+1, N+1, byrow = F))
         Z.fitted = Z.fitted + Zhat 
      }
      
      # kappa = kappa.out
      if (i%%100 == 0){
         print(c("interation number:", i))
         print(c("N: ", N))
      }
   }
   Z.fitted = Z.fitted / mcmc
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list, prob_list = prob_list, Z.fitted = Z.fitted))
}



