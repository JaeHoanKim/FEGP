rm(list = ls())
gc()
library(fields)
library(FastGP)
library(Rcpp)
library(ggplot2)
library(Matrix)
library(foreach)
library(doParallel)
library(tidyverse)
source("Result_Manuscript/GraphAesthetics.R")

### 1. true function setting & data generation

alpha = 1.4

f0_1D = function(x, trun = 500){
   value = 0
   for(j in 1:trun){
      value = value + sin(j) * cos(pi * (j - 1/2) * x) * j^(- alpha - 1)
   }
   return(value * sqrt(2))
}
# f0_1D = function(x){return(x^2 + sin(x))}

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
}


target = 2500
brn = 0
brn.ESS = 100
# setting for the Matern parameters
kappak = seq(1, 5, 0.5)
tausqk = 1
Nk = c(4, 6, 8, 10, 12)
kappa.pr = tausq.pr = N.pr = const 
beta = 4

const = function(x){
   return(1)
}

grid.plot = c(0:1000)/1000

### 2. MSE calculation - 1D

source("1D/GPI/functions_GPI.R")
source("1D/GPI/functions_GPI_sampling.R")
MSE_list_GPI1D = matrix(nrow = M, ncol = length(nlist))

### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers/2)
registerDoParallel(cl)
##########################
for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   df = df_1D[[a]]
   output <- foreach (m = 1:M, .packages = c("Matrix", "fields", "FastGP")) %dopar% {  
      # m th dataset among M = 50 dataset
      X = df$X[((m-1)*n+1):(m*n)]
      Y = df$Z[((m-1)*n+1):(m*n)]
      result = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr,
                              kappak = kappak, kappa.pr = kappa.pr, 
                              tausqk = tausqk, tausq.pr = tausq.pr, 
                              beta = beta, mcmc = target, brn=0, brn.ESS = brn.ESS)
      g.plot = tail(result$g_list, target) # choosing last `target` samples
      obs = data.frame(X, Y)
      y.plot = glist_to_plotdf(g.plot, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
      print(m)
      mean((y.plot$true - y.plot$mean)^2)
   }
   MSE_list_GPI1D[, a] = simplify(output)
   print(a)
}
stopCluster(cl)

save(MSE_list_GPI1D, file = "MSE_1D_GPI_1.RData")