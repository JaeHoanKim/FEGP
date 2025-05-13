## comparison between original GP and d-frGPloc SPDE
rm(list = ls())
gc()
library(gridExtra)
library(fields)
library(FastGP)
library(Rcpp)
library(ggplot2)
library(Matrix)
library(foreach)
library(rSPDE)
library(doParallel)
library(tidyverse)
# sourceCpp("1D/GPI/inv_chol.cpp")

### 1. true function setting & data generation
alpha_val = 7
f0_1D = function(x, alpha = alpha_val, trun = 500){
   # value = 0
   # for(j in 1:trun){
   #    value = value + sin(j) * cos(pi * (j - 1/2) * x) * j^(- alpha - 1)
   # }
   # return(value * sqrt(2))
   return (2^alpha * abs(0.5 - x)^alpha)
}
const = function(x){return(1)}

M = 50
nlist = c(200, 500, 1000)
df_1D = list(length = length(nlist))

for(i in 1:length(nlist)){
   set.seed(i)
   n = nlist[i]
   # 1D data generation
   X = runif(n*M)
   Z = f0_1D(X) + rnorm(n*M) * 0.1
   df_1D[[i]] = data.frame(X, Z)
   # 2D data generation
   X = matrix(runif(2*n*M), n*M)
}


target = 2500
brn = 0
brn.ESS = 100
# setting for the Matern parameters
# kappak = seq(0.2, 4, length.out = 10)
kappa_cand = seq(5, 10, length.out = 8)
tausqk = 1
Nk = c(4, 6, 8, 10, 12, 16, 20, 25, 30)
kappa.pr = tausq.pr = N.pr = const 
beta = 2

### 2. MSE calculation - 1D

source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")
source("1D/MaternGP1D.R")
grid.plot = c(0:300)/300

### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers/2)
registerDoParallel(cl)

for (kappak in kappa_cand){
   MSE_list_SPDE1D = matrix(nrow = M, ncol = length(nlist))
   MSE_list_Matern1D = matrix(nrow = M, ncol = length(nlist))
   for(a in 1:length(nlist)){
      n = nlist[a]
      df = df_1D[[a]]
      l = 1/kappak[1]
      v = alpha_val
      output <- foreach (m = 1:M, .packages = c("Matrix", "fields", "FastGP")) %dopar% {
         X = df$X[((m-1)*n+1):(m*n)]
         Y = df$Z[((m-1)*n+1):(m*n)]
         result = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr, beta = beta,
                                   kappak = kappak, kappa.pr = kappa.pr,
                                   tausqk = tausqk, tausq.pr = tausq.pr,
                                   mcmc = target, brn=0, seed = 1234)
         g.plot = tail(result$g_list, target) # choosing last `target` samples
         obs = data.frame(X, Y)
         y.plot = glist_to_plotdf(g.plot, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
         mean((y.plot$true - y.plot$mean)^2)
      }
      MSE_list_SPDE1D[, a] = simplify(output)
      ## MSE for original Matern
      output2 <- foreach (m = 1:M) %dopar% {
         X = df$X[((m-1)*n+1):(m*n)]
         Y = df$Z[((m-1)*n+1):(m*n)]
         tmp = sample_posterior(X, Y, grid.plot, l, v, sigma_n = 0.1, n_samples = 20)
         y.plot.Matern = matrix_to_plotdf(tmp, truefun = f0_1D, grid = grid.plot)
         mean((y.plot.Matern$true - y.plot.Matern$mean)^2)
      }
      MSE_list_Matern1D[, a] = simplify(output2)
      print(a)
   }
   aa = colMeans(MSE_list_SPDE1D)
   bb = colMeans(MSE_list_Matern1D)
   
   MSE_list_1D = rbind(MSE_list_SPDE1D, MSE_list_Matern1D)
   
   filename = paste0("MSE_comparison/MSE_list_SPDE_Matern_1D_kappa", round(kappak, 2), ".RData")
   save(MSE_list_1D, file = filename)  
   print(kappak)
}


stopCluster(cl)