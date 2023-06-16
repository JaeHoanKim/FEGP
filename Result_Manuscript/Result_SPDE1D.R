rm(list = ls())
gc()
library(grDevices)
library(ggplot2)
library(Matrix)
source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")

target = 2500
brn.ESS = 100
kappa = 2

const = function(x){
   return(1)
}

f0_1D = function(x){return(1*x^2+sin(8*x))}
alpha1 = 0.9 # coverage probability 1 (dark gray area)
alpha2 = 0.95 # coverage probability 2 (light gray area)
Nk = c(4, 6, 8, 10, 12)

M = 50
nlist = c(200, 500, 1000)
MSE_list = matrix(nrow = M, ncol = length(nlist))

for(a in 1:length(nlist)){
   n = nlist[a]
   filename = paste0("Result_Manuscript/obs_n1D", n, ".RData")
   load(filename)
   for(m in 1:M){
      # m th dataset among M = 50 dataset
      X = df$X[((m-1)*n+1):(m*n)]
      Y = df$Z[((m-1)*n+1):(m*n)]
      result = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = const,
                                kappa.init = kappa, mcmc = target, brn=brn, seed = 1234)
      g.plot = tail(result$g_list, target) # choosing last `target` samples
      grid.plot = c(0:1000)/1000
      obs = data.frame(X, Y)
      y.plot = glist_to_plotdf(g.plot, grid.plot, true = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
      MSE_list[m, a] = mean((y.plot$true - y.plot$mean)^2)
      print(m)
   }
}
## Plot 3. Fitting of the samples (credible intervals)

################## converting predicted values into usable form ##################
##################################################################################

obs = data.frame(X, Y)

ggplot(y.plot, aes(x = x)) +
   geom_line(aes(y=mean), colour="blue") + 
   geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) + 
   geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
   geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
   geom_point(data = obs, aes(X, Y), size = 0.3) +
   labs(title = paste0("CI 95% of ", target, " samples (n = ", n, ")"))+
   geom_point(data = obs, aes(X, Y), size = 0.3)


N_list = tail(result$N_list, target)
plot(N_list)
lines(N_list)
MSE = mean((y.plot$true - y.plot$mean)^2)
