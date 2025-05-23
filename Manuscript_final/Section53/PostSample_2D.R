rm(list = ls())
gc()
library(tidyverse)
library(fields)
library(FastGP)
library(ggplot2)
library(Matrix)
library(doParallel)
library(foreach)
library(doParallel)
library(tidyverse)
library(invgamma)
source("2D/FullGP/functions_FullGP.R")

### 1. true function setting & data generation
# alpha = 0.5
# f0_2D = function(x, y){return(x^2 + sqrt(abs(y-0.5)) + sin(8*x))}
# f0_2D = function(x, y){return(sin(5*x + 2*y) + 2*y^2)}
# For comparison between GPI and NNGP, we choose f0 such that GPI can do better.
# As \nu \infty in Matern makes Matern similar to RBF kernel, we use this.
f0_2D = function(x, y){return(sin(5*x + 2*y) + 2*y^2)}
M = 50
nlist = c(200, 500, 1000)
df_2D = list(length = length(nlist))
sigma = 0.1
for(i in 1:length(nlist)){
   set.seed(i)
   n = nlist[i]
   # 2D data generation
   X = matrix(runif(2*n*M), n*M)
   Z = f0_2D(X[, 1], X[, 2]) + rnorm(n*M) * sigma
   df_2D[[i]] = data.frame(X, Z)
}

# setting for the Matern parameters and sampling

beta = 2
d = 2
nu = beta - d/2

kappak = 30
tausqk = 1
kappa.pr = function(x){return(dgamma(x, 5, 1/5))}
kappa.sampler = function(){rgamma(1, 5, 1/5)}
tausq.pr = function(x){return(dgamma(x, 1, 1))}
tausq.sampler = function(){return (1)}
tausq = tausq.sampler()

MSE_list_FullGP2D = matrix(nrow = M, ncol = length(nlist))

target = 2500
brn = 1000
# brn.ESS = 100

### 3. MSE calculation - 2D

gridsize = 40
# gridmat is a (gridsize^2) by 2 matrix!
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

### Parallel computing ###
nworkers <- 5 # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)

# Exact posterior distribution of full GP
for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_2D[[a]]
   output <- foreach (m = 1:M, .packages = c("Matrix", "fields", "FastGP")) %dopar% {
      X = as.matrix(df[((m-1)*n+1):(m*n), c(1, 2)])
      Z = as.matrix(df$Z[((m-1)*n+1):(m*n)])
      PostMean <- sample_fullGP(X, as.vector(Z), sigma = sigma, nu = nu, 
                                kappa_pr = kappa.pr, kappa_sampler = kappa.sampler, 
                                tausq_pr = tausq.pr, tausq_sampler = tausq.sampler, grid_size = 40, target, brnin = brn, kappa_fixed = kappak)$PostMean
      truefun = f0_2D(gridmat[, 1], gridmat[, 2])
      mean((truefun - PostMean)^2)
   }
   MSE_list_FullGP2D[, a] = purrr::simplify(output)
   print(a)
}

stopCluster(cl)

############## GPI 2D ###############

MSE_list_GPI2D = matrix(nrow = M, ncol = length(nlist))

source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")

## Param check ##
N_supp = c(6, 10, 14, 18)
N.pr = function(x){return(rep(1, length(x)))}
kappa_supp = c(20, 30, 40)
tausq_supp = 1
param_check_list = vector("list", length = length(nlist))
for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_2D[[a]]
   X = df[1:n, c(1, 2)]
   Z = df$Z[1:n]
   param_check_list[[a]] = param_check_2D(X, Z, beta, N_supp, kappa_supp, tausq_supp,  N.pr, kappa.pr, tausq.pr, sigsq = sigma^2)
}
param_check_list[[1]]$N_list
param_check_list[[1]]$kappa_list
### Parallel computing ###
cl <- makeCluster(nworkers)
registerDoParallel(cl)

for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_2D[[a]]
   output <- foreach (m = 1:M, .packages = c("Matrix", "fields", "FastGP")) %dopar% {
      X = df[((m-1)*n+1):(m*n), c(1, 2)]
      Z = df$Z[((m-1)*n+1):(m*n)]
      # result = sample.RJESS2D.seq(Z = Z, X = X, Nk = Nk, N.pr = N.pr, 
      #                             kappak = kappak, kappa.pr = kappa.pr, 
      #                             tausqk = tausqk, tausq.pr = tausq.pr, sigsq = sigma^2, beta = beta,
      #                             mcmc = target, brn = 0, brn.ESS = brn.ESS)
      result = sample.GPI2D(Z = Z, X = X, Nk = N_supp, N.pr = N.pr, 
                            kappa.pr = kappa.pr, kappa.sampler = kappa.sampler,
                            tausq.pr = tausq.pr, tausq.sampler = tausq.sampler, 
                            sigsq = sigma^2, beta = beta, mcmc = target, brn = 0, brn.ESS = brn.ESS)
      g_list = result$g_list
      y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)
      mean((y.plot$truefun - y.plot$mean)^2)
   }
   print(a)
   print(Sys.time())
   MSE_list_GPI2D[, a] = purrr::simplify(output)
}

stopCluster(cl)
save(MSE_list_GPI2D, file = "MSE_comparison/MSE_list_generated_data_2D_GPI.RData")


MSE_list_2D = rbind(MSE_list_GPI2D, MSE_list_FullGP2D)
save(MSE_list_2D, file = "MSE_comparison/MSE_list_generated_data_2D_GPI_FullGP.RData")