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
source("Result_Manuscript/GraphAesthetics.R")
sourceCpp("1D/GPI/inv_chol.cpp")

### 1. true function setting & data generation

alpha = 1.4
f0_1D = function(x, trun = 500){
   value = 0
   for(j in 1:trun){
      value = value + sin(j) * cos(pi * (j - 1/2) * x) * j^(- alpha - 1)
   }
   return(value * sqrt(2))
}

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
tausqk = 1
Nk = c(6, 8, 10, 14, 18, 22, 26, 30)
kappa.pr = tausq.pr = N.pr = const 
beta = 4

const = function(x){
   return(1)
}


### 2. MSE calculation - 1D

source("1D/GPI/functions_GPI.R")
source("1D/GPI/functions_GPI_sampling.R")
m = 1
log_prob_list_GPI = list()
prob_N_GPI = list()
prob_kappa_GPI = list()
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
   log_prob_mat = matrix(log_prob_list_GPI[[a]], nrow = length(kappak))
   log_prob_mat = log_prob_mat - max(log_prob_mat)
   prob_mat = exp(log_prob_mat)
   prob_N_GPI[[a]] = colSums(prob_mat) / sum(prob_mat)
   prob_kappa_GPI[[a]] = colSums(t(prob_mat)) / sum(prob_mat)
}

source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")
log_prob_list_SPDE = list()
prob_N_SPDE = list()
prob_kappa_SPDE = list()
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
   # calculate the marginal likelihood of N
   log_prob_mat = matrix(log_prob_list_SPDE[[a]], nrow = length(kappak))
   log_prob_mat = log_prob_mat - max(log_prob_mat)
   prob_mat = exp(log_prob_mat)
   prob_N_SPDE[[a]] = colSums(prob_mat) / sum(prob_mat)
   prob_kappa_SPDE[[a]] = colSums(t(prob_mat)) / sum(prob_mat)
}

prob_N_SPDE <- data.frame(prob_N = unlist(prob_N_SPDE), n = rep(nlist, each = length(Nk)), 
                        N = rep(Nk, times = length(nlist))) %>%
   mutate(n = factor(n))
ggplot(prob_N_SPDE, aes(x = N, y = prob_N)) + geom_point() + geom_line(aes(color = n)) + theme1
   
prob_kappa_SPDE <- data.frame(prob_kappa = unlist(prob_kappa_SPDE), n = rep(nlist, each = length(kappak)), 
                          kappa = rep(kappak, times = length(nlist))) %>%
   mutate(n = factor(n))
ggplot(prob_kappa_SPDE, aes(x = kappa, y = prob_kappa)) + geom_point() + geom_line(aes(color = n)) + theme1

