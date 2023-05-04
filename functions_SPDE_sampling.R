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
