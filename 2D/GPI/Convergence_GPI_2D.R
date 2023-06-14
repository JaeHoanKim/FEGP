rm(list = ls())
gc()
library(plotly)
library(tidyverse)
library(ggpubr)
library(Rcpp)
library(ggplot2)
source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")

n = 200 # the number of observed data; 200, 500, 1000
filename = paste0("Result_Manuscript/obs_n", n, ".RData")
load(filename)
m = 1 # m th dataset among M = 50 dataset
X = df[((m-1)*n+1):(m*n), c(1, 2)]
Z = df$Z[((m-1)*n+1):(m*n)]


kappa = 2
N.init = 10
brnin = 100
target = 100
const = function(x){
   return(1)
}
dpoi5 = function(x){
   return(dpois(x, lambda = 5))
}
## convergence diagnosis with a fixed N
N = 6
result1 = sample.ESS.Nfixed2D(Z = Z, X = X, l.in = 1/kappa, nu.in = 1, 
                             mcmc = target, brn = brnin, sigsq = 0.1^2, N.init = c(N, N), tausq = 1)
result2 = sample.ESS.Nfixed2D(Z = Z, X = X, l.in = 1/kappa, nu.in = 1, 
                              mcmc = target, brn = brnin, sigsq = 0.1^2, N.init = c(N, N), tausq = 1)
result3 = sample.ESS.Nfixed2D(Z = Z, X = X, l.in = 1/kappa, nu.in = 1, 
                              mcmc = target, brn = brnin, sigsq = 0.1^2, N.init = c(N, N), tausq = 1)
# install.packages("stableGR")
library(stableGR)
chain1 = do.call(rbind, result1$g_list)
chain2 = do.call(rbind, result2$g_list)
chain3 = do.call(rbind, result3$g_list)
## simplify reuslt$g_list into matrices
CODA = stableGR::stable.GR(x = list(chain1, chain2, chain3))
