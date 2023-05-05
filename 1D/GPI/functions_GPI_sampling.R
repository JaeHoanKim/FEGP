
#########################################################################
########## Function for MCMC samples using ESS and WC algoritm ##########
#########################################################################


## goal: draw samples of (N, g^{(N)}, \theta) given the data D_n (X_n, Y_n)
## input: X, Y, sigma 
## output: mcmc samples after burn-in
## without considering the error term sigma
## knitting v (size N vector) and g_N (size N+1 vector) in an alternative order

sample.ESS.nest = function(X, Y, N.pr, l.in = 0.5, nu.in = 0.75, mcmc, brn, thin, sigsq, N.init = 2, g.init, tausq = 1){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   if(length(X)!=length(Y))
      stop("X and Y should be of same length!")
   n = length(Y)
   ## initial value of N
   ## initial value of g
   g.init = rep(0, N.init + 1)
   em = mcmc + brn # total number of sampling
   ## initialize the collection of g at first step
   g.in = list()
   g.out = list()
   g_list = list()
   N_list = c()
   # Nk: the N value (nrow(g) - 1) of the k-th chain!
   for(k in 1:length(Nk)){
      g.in[[k]] = Sampling_N_new(Nk[k], kappa)/sqrt(mv)
   }
   # setting up for iterative work
   g.in = g.init
   N = N.init
   for(i in 1:em){
      ### update knot_N for updated N
      knot_N = c(0:N)/N
      ### which stage to run 
      ## if u < -2 : consider N/2; else if u < -1: consider staying; else, consider double N
      uu = ((N==2) - 3) * runif(1)
      if (uu < - 2){
         # consider half N
         g.cand.half = g.in[2*c(0:(N/2)) + 1]
         log.test.prob = min(0, log(N.pr(N/2) / N.pr(N)) + loglik(Y, X, g.cand.half, sigsq) - 
                                loglik(Y, X, g.in, sigsq))
         u = runif(1)
         if (exp(log.test.prob) > u){
            g.out = g.cand.half
            N.out = N/2
         } else{
            g.out = g.in
            N.out = N
         }
      } else if (uu < -1){
         # consider N staying
         ## sample new g given N, D_n
         nu.ess = samp.WC(knot_N, nu.in, l.in, tausq)
         g.out = ESS(g.in, nu.ess, Y, X, sigsq)
         N.out = N
         log.test.prob = 1
      } else{
         # consider double N
         ### applying ESS for sampling v : when g is given, make v to be used for N doubling procedure
         nu_ESS_v = samp.WC(knot_N, nu.in, l.in, tausq)[1:N]
         v.in = samp.WC(knot_N, nu.in, l.in, tausq)[1:N] 
         v.in = v_g_ESS(v.in, g.in, nu_ESS_v, nu.in, l.in)
         ### procedure to decide whether to double N or stay
         g.cand.2 = knit_vec(g.in, v.in)
         log.test.prob = min(0, log(N.pr(2*N) / N.pr(N)) +  loglik(Y, X, g.cand.2, sigsq) - 
                                loglik(Y, X, g.in, sigsq))
         u = runif(1)
         if (exp(log.test.prob) > u){
            ## if g is updated to 2N+1 vector, nu should start again from the first of the loop
            g.out = g.cand.2
            N.out = 2*N
         } else{
            g.out = g.in
            N.out = N
         }
      }
      ### g.out and N.out are under compliance
      g_list[[i]] = g.out
      N_list[i] = N.out
      ### preparation for iteration
      g.in = g.out
      N = N.out
      prob_list[i] = exp(log.test.prob)
      if (i%%1000 == 0){
         print(c("interation number:", i))
         print(c("N: ", N))
      }
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list, prob_list = prob_list))
}



sample.PTESS = function(X, Y, N.pr, Nk, Tk,
                        l.in = 0.5, nu.in = 0.75, mcmc, brn, thin, sigsq, tausq = 1){
   ## X, Y: given data
   ## N.pr: prior distribution of N (function)
   ## l.in, nu.in: initial value of l and nu (does not change throughout the simulation)
   if(length(X)!=length(Y))
      stop("X and Y should be of same length!")
   n = length(Y)
   Nk = sort(Nk)
   if(length(Nk) != length(Tk)){
      stop("The length of Nk and Tk should be the same!")
   }
   if (Tk[1] != 1){
      stop("The first element of Tk should be 1!")
   }
   ## initialize the collection of g at first step
   g.in = list()
   g.out = list()
   g_list = list()
   N_list = c()
   em = mcmc + brn # total number of sampling
   N_list = c()
   prob_list = c()
   g_list = list()
   # starting with the initial N and g
   # setting up for iterative work
   g.in = g.init
   N = N.init
   for(i in 1:em){
      ### update knot_N for updated N
      knot_N = c(0:N)/N
      ### which stage to run 
      ## if u < -2 : consider N/2; else if u < -1: consider staying; else, consider double N
      uu = ((N==2) - 3) * runif(1)
      if (uu < - 2){
         # consider half N
         g.cand.half = g.in[2*c(0:(N/2)) + 1]
         log.test.prob = min(0, log(N.pr(N/2) / N.pr(N)) + loglik(Y, X, g.cand.half, sigsq) - 
                                loglik(Y, X, g.in, sigsq))
         u = runif(1)
         if (exp(log.test.prob) > u){
            g.out = g.cand.half
            N.out = N/2
         } else{
            g.out = g.in
            N.out = N
         }
      } else if (uu < -1){
         # consider N staying
         ## sample new g given N, D_n
         nu.ess = samp.WC(knot_N, nu.in, l.in, tausq)
         g.out = ESS(g.in, nu.ess, Y, X, sigsq)
         N.out = N
         log.test.prob = 1
      } else{
         # consider double N
         ### applying ESS for sampling v : when g is given, make v to be used for N doubling procedure
         nu_ESS_v = samp.WC(knot_N, nu.in, l.in, tausq)[1:N]
         v.in = samp.WC(knot_N, nu.in, l.in, tausq)[1:N] 
         v.in = v_g_ESS(v.in, g.in, nu_ESS_v, nu.in, l.in)
         ### procedure to decide whether to double N or stay
         g.cand.2 = knit_vec(g.in, v.in)
         log.test.prob = min(0, log(N.pr(2*N) / N.pr(N)) +  loglik(Y, X, g.cand.2, sigsq) - 
                                loglik(Y, X, g.in, sigsq))
         u = runif(1)
         if (exp(log.test.prob) > u){
            ## if g is updated to 2N+1 vector, nu should start again from the first of the loop
            g.out = g.cand.2
            N.out = 2*N
         } else{
            g.out = g.in
            N.out = N
         }
      }
      ### g.out and N.out are under compliance
      g_list[[i]] = g.out
      N_list[i] = N.out
      ### preparation for iteration
      g.in = g.out
      N = N.out
      prob_list[i] = exp(log.test.prob)
      if (i%%1000 == 0){
         print(c("interation number:", i))
         print(c("N: ", N))
      }
   }
   ## when pred is FALSE, Ypred would be the zero matrix
   return(list(g_list = g_list, N_list = N_list, prob_list = prob_list))
}
