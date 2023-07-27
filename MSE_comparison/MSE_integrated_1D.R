rm(list = ls())
gc()
library(tidyverse)
library(fields)
library(FastGP)
library(Rcpp)
library(ggplot2)
library(Matrix)
library(doParallel)
library(foreach)
library(rSPDE)
library(doParallel)
library(tidyverse)
library(spNNGP)
sourceCpp("1D/GPI/inv_chol.cpp")

### 1. true function setting & data generation

alpha = 0.5

f0_1D = function(x, trun = 500){
   value = 0
   for(j in 1:trun){
      value = value + sin(j) * cos(pi * (j - 1/2) * x) * j^(- alpha - 1)
   }
   return(value * sqrt(2))
}


M = 50
nlist = c(200, 500, 1000)
target = 2500
brn = 0
brn.ESS = 1000
# setting for the Matern parameters
kappa = 2
beta = 20
d = 1
nu = beta - d/2
l.in = 1/kappa

const = function(x){
   return(1)
}
Nk = c(4, 8, 12, 20, 30, 40, 50)
grid.plot = c(0:1000)/1000
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

### 2. MSE calculation - 1D

source("1D/GPI/functions_GPI.R")
source("1D/GPI/functions_GPI_sampling.R")
MSE_list_GPI1D = matrix(nrow = M, ncol = length(nlist))

### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)
##########################
for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   df = df_1D[[a]]
   output <- foreach (m = 1:M, .packages = c("Matrix", "fields", "FastGP")) %dopar% {  
      # m th dataset among M = 50 dataset
      X = df$X[((m-1)*n+1):(m*n)]
      Y = df$Z[((m-1)*n+1):(m*n)]
      result = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, nu.in = nu, l.in = 1/kappa,
                              N.pr = const,
                              mcmc = target, brn=0, brn.ESS = brn.ESS)
      g.plot = tail(result$g_list, target) # choosing last `target` samples
      obs = data.frame(X, Y)
      y.plot = glist_to_plotdf(g.plot, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
      print(m)
      mean((y.plot$true - y.plot$mean)^2)
   }
   MSE_list_GPI1D[, a] = simplify(output)
}
stopCluster(cl)


source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")
MSE_list_SPDE1D = matrix(nrow = M, ncol = length(nlist))
grid.plot = c(0:1000)/1000

### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)

for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_1D[[a]]
   output <- foreach (m = 1:M, .packages = c("Matrix", "fields", "FastGP")) %dopar% {
      X = df$X[((m-1)*n+1):(m*n)]
      Y = df$Z[((m-1)*n+1):(m*n)]
      result = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = const, beta = 2,
                                kappa.init = kappa, mcmc = target, brn=0, seed = 1234)
      g.plot = tail(result$g_list, target) # choosing last `target` samples
      obs = data.frame(X, Y)
      y.plot = glist_to_plotdf(g.plot, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
      mean((y.plot$true - y.plot$mean)^2)
   }
   MSE_list_SPDE1D[, a] = simplify(output)
}

stopCluster(cl)

MSE_list_1D = rbind(MSE_list_GPI1D, MSE_list_SPDE1D)
save(MSE_list_1D, file = "MSE_list_generated_data_1D.RData")