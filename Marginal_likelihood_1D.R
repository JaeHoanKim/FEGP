## Marginal likelihood for 1D
# comparing GPI and SPDE according to parameters

rm(list = ls())
gc()
library(fields)
library(FastGP)
library(Rcpp)
library(ggplot2)
library(Matrix)
library(foreach)
library(rSPDE)
library(doParallel)
library(tidyverse)
sourceCpp("1D/GPI/inv_chol.cpp")

### 1. true function setting & data generation

alpha = 1.4

# f0_1D = function(x, trun = 500){
#    value = 0
#    for(j in 1:trun){
#       value = value + sin(j) * cos(pi * (j - 1/2) * x) * j^(- alpha - 1)
#    }
#    return(value * sqrt(2))
# }
f0_1D = function(x){return(x^2 + sin(x))}

const = function(x){return(1)}


M = 1
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
brn.ESS = 1000
# setting for the Matern parameters
kappak = seq(2, 5, 0.3)
tausqk = seq(1, 5, 0.5)
Nk = c(6, 8, 10, 14, 18, 22, 26, 30)
kappa.pr = tausq.pr = N.pr = const 
beta = 4

const = function(x){
   return(1)
}

grid.plot = c(0:1000)/1000

### 2. MSE calculation - 1D

source("1D/GPI/functions_GPI.R")
source("1D/GPI/functions_GPI_sampling.R")
m = 1
log_prob_list_GPI = list()
##########################
for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   df = df_1D[[a]]
   # m th dataset among M = 50 dataset
   X = df$X[((m-1)*n+1):(m*n)]
   Y = df$Z[((m-1)*n+1):(m*n)]
   result = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr,
                           kappak = kappak, kappa.pr = kappa.pr, 
                           tausqk = tausqk, tausq.pr = tausq.pr, 
                           beta = beta, mcmc = target, brn=0, brn.ESS = brn.ESS)
   log_prob_list_GPI[[a]] = result$log_prob_N_list
}

source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")
log_prob_list_SPDE = list()
grid.plot = c(0:1000)/1000
m = 1
for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_1D[[a]]
   X = df$X[((m-1)*n+1):(m*n)]
   Y = df$Z[((m-1)*n+1):(m*n)]
   result = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr, beta = beta,
                             kappak = kappak, kappa.pr = kappa.pr,
                             tausqk = tausqk, tausq.pr = tausq.pr,
                             mcmc = target, brn=0, seed = 1234)
   log_prob_list_SPDE[[a]] = result$log_prob_N_list
}

## calculate the marginal likelihood of N

stopCluster(cl)

MSE_list_1D = rbind(MSE_list_GPI1D, MSE_list_SPDE1D)
save(MSE_list_1D, file = "MSE_list_generated_data_1D_2.RData")