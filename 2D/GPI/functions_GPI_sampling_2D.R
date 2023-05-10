library(Rcpp)
#########################################################################
########## Function for MCMC samples using ESS and WC algoritm ##########
#########################################################################

## knitting v (size N vector) and g_N (size N+1 vector) in an alternative order
sample.ESS.Nfixed2D = function(Z, X, l.in, nu.in, mcmc, brn, thin, sigsq, N.init, g.init, tausq, pred = FALSE, Xtest = NULL, gridsize = 500){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   ## l: range of the fields::Matern function
   if(length(Z) != dim(X)[1])
      stop("Z and X should have the same length!")
   ## follow the initial value given in Ray (2020)
   if(missing(nu.in))
      nu.in=1
   if(missing(l.in))
      l.in=0.5
   #l_est(nu.in,c(0,1),0.05)
   ## initial value of N
   if(missing(N.init)){
      N.init = c(10, 10)
   } else{
      if(length(N.init) != 2){
         stop("N.init should have two number which denotes the number of of grid edges!")
      }
   }
   ## initial value of g
   if(missing(g.init)){
      g.init = matrix(0, N.init + 1, N.init + 1)
   } else{
      if (length(g.init) != N.init){
         stop("provided g.init should have the length of N.init!")
      }
   }
   if(missing(sigsq)){
      sigsq = 1
   }
   if(missing(tausq)){
      tausq = 1
   }
   if(missing(mcmc))
      mcmc=5000
   if(missing(brn))
      brn=1000
   if(missing(thin))
      thin=1
   em = mcmc + brn # total number of sampling
   ## compatibility check when predicting
   N_list = list()
   prob_list = c()
   g_list = list()
   Zgrid_list = list()
   # starting with the initial N and g
   # setting up for iterative work
   g.in = g.init
   N = N.init
   # since N does not change, we can fix 'result'.
   result = eigvals_exact(ndim = N, nu = nu.in, lambda_g = l.in)
   if (pred == TRUE){
      Xgrid = cbind(rep(c(0:gridsize)/gridsize, gridsize + 1), 
                    rep(c(0:gridsize)/gridsize, each = gridsize + 1)) ## making grid points from (0, 0) to (1, 1) 
   }
   for(i in 1:em){
      ### update knot_N for updated N
      ## sample new g given N, D_n
      nu.ess = matrix(samp_from_grid(N, mdim = result$mvec, egs = result$egvals, nu, lambda_g), N + 1)
      g.out = ESS(g.in, nu.ess, z = Z, x = X, sigsq)
      N.out = N
      if (i>brn){
         g_list[[i-brn]] = g.out
         if (pred == TRUE){
            Zgrid = matrix(f_N_h_2D_multi(Xgrid, g.in), gridsize + 1, gridsize + 1)
            Zgrid_list[[i-brn]] = Zgrid
         }
         N_list[[i]] = N.out
      }
      ### preparation for iteration
      g.in = g.out
      N = N.out
      if (i%%100 == 0){
         print(c("interation number:", i))
      }
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list, Zgrid_list = Zgrid_list))
}



## goal: draw samples of (N, g^{(N)}, \theta) given the data D_n (X_n, Y_n)
## input: X, Y, sigma 
## output: mcmc samples after burn-in
## without considering the error term sigma

## knitting v (size N vector) and g_N (size N+1 vector) in an alternative order
sample.ESS.nested2D = function(Z, X, N.pr, mcmc, brn, thin, l.in = NULL, nu.in = NULL, sigsq, N.init, g.init, tausq, 
                 pred = FALSE, Xtest = NULL, gridsize = 500){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   if(length(Z) != dim(X)[1])
      stop("Z and X should have the same length!")
   ## follow the initial value given in Ray (2020)
   if(is.null(nu.in))
      nu.in=0.75
   if(is.null(l.in))
      l.in=0.5
   #l_est(nu.in,c(0,1),0.05)
   ## initial value of N
   if(missing(N.init)){
      N.init = c(2, 2)
   } else{
      if(length(N.init) != 2){
         stop("N.init should have two number which denotes the number of of grid edges!")
      }
   }
   ## initial value of g
   if(missing(g.init)){
      g.init = matrix(0, N.init + 1, N.init + 1)
   } else{
      if (length(g.init) != N.init){
         stop("provided g.init should have the length of N.init!")
      }
   }
   if(missing(sigsq)){
      sigsq = 1
   }
   if(missing(tausq)){
      tausq = 1
   }
   if(missing(mcmc))
      mcmc=5000
   if(missing(brn))
      brn=1000
   if(missing(thin))
      thin=1
   em = mcmc + brn # total number of sampling
   ## compatibility check when predicting
   N_list = list()
   prob_list = c()
   g_list = list()
   Zgrid_list = list()
   # starting with the initial N and g
   # setting up for iterative work
   g.in = g.init
   N = N.init
   # since N does not change, we can fix 'result'.
   result = eigvals_exact(ndim = N, nu = nu.in, lambda_g = l.in)
   if (pred == TRUE){
      Xgrid = cbind(rep(c(0:gridsize)/gridsize, gridsize + 1), 
                    rep(c(0:gridsize)/gridsize, each = gridsize + 1)) ## making grid points from (0, 0) to (1, 1) 
   }
   # store result_N and result_2N to prevent re-calculation of the inverse of BTTB
   result_N_list = list()
   for(i in 1:6){
      N1 = 2 * 2^(i-1)
      N2 = 2 * 2^(i-1)
      X_N = cbind(rep(c(0:N1)/N1, N2 + 1), 
                  rep(c(0:N2)/N2, each = N1 + 1)) # (N1+1)*(N2+1) by 2 matrix
      Sigma_N = thekernel(X_N, nu.in, lambda = l.in)
      result_N_list[[N1]] = inv_BTTB(Sigma_N, N1 + 1)
   }
   
   for(i in 1:em){
      ### update knot_N for updated N
      ## sample new g given N, D_n
      jump.prob = runif(1)
      log.test.prob = NULL
      if (min(N) == 2){
         u = 1; log.test.prob = -100
      } else{
         ## check if accept half N
         result_half = eigvals_exact(ndim = N / 2, nu = nu.in, lambda_g = l.in)
         g.cand.half = g.in[2 * c(0:(N[1]/2)) + 1, 2 * c(0:(N[2]/2)) + 1]
         log.test.prob = loglik(Z, X, g.cand.half, sigsq) - loglik(Z, X, g.in, sigsq) +
            log(N.pr(N[1]/2)) - log(N.pr(N[1]))
         u = runif(1)
         if (u < exp(log.test.prob)){
            g.out = g.cand.half
            N.out = N / 2
            result = result_half
         }
      }
      if (min(N) == 2 | u >= exp(log.test.prob)){
         ## if cutting in half is rejected,
         ## check if accept doubling or not
         ## 1. sample from posterior distribution of g, (2N+1) by (2N+1) matrix
         result2 = eigvals_exact(ndim = 2 * N, nu = nu.in, lambda_g = l.in)
         gdouble = matrix(samp_from_grid(2*N, mdim = result2$mvec, egs = result2$egvals, nu.in, lambda_g = l.in), 2*N + 1)
         gdouble_ess = matrix(samp_from_grid(2*N, mdim = result2$mvec, egs = result2$egvals, nu.in, lambda_g = l.in), 2*N + 1)
         g.cand.2 = ESS_double(gdouble, g.in, gdouble_ess, nu.in, l.in, 
                               result_N = result_N_list[[nrow(g.in) - 1]], 
                               result_2N = result_N_list[[nrow(gdouble) - 1]])
         ## 2. likelihood check
         log.test.prob = loglik(Z, X, g.cand.2, sigsq) - loglik(Z, X, g.in, sigsq) +
            log(N.pr(N[1]*2)) - log(N.pr(N[1])) 
         u = runif(1)
         if (u < exp(log.test.prob)){
            g.out = g.cand.2
            N.out = 2 * N
            result = result2
         }else{
            nu.ess = matrix(samp_from_grid(N, mdim = result$mvec, egs = result$egvals, nu, lambda_g = l.in), N + 1)
            g.out = ESS(g.in, nu.ess, z = Z, x = X, sigsq)
            N.out = N
         }
      }
      prob_list[[i]] = log.test.prob
      if (i>brn){
         g_list[[i-brn]] = g.out
         if (pred == TRUE){
            Zgrid = matrix(f_N_h_2D_multi(Xgrid, g.in), gridsize + 1, gridsize + 1)
            Zgrid_list[[i-brn]] = Zgrid
         }
         N_list[[i]] = N.out
      }
      ### preparation for iteration
      g.in = g.out
      N = N.out
      if (i%%10 == 0){
         print(c("interation number:", i))
         print(N)
      }
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list, Zgrid_list = Zgrid_list, prob_list = prob_list))
}


