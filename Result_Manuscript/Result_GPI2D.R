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
## discretized version of 1 over exponential distribution - which satisfy the condition for prior theoretically
# N.pr = function(N){return (1/N^2 * 3 * exp(-3 / N))}

kappa = 2
N.init = 10
brnin = 1000
target = 2500
const = function(x){
   return(1)
}
dpoi5 = function(x){
   return(dpois(x, lambda = 5))
}
Nk = c(4, 6, 8, 10, 12)
result = sample.RJESS2D(Z = Z, X = X, N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                        mcmc = target, brn = brnin, nu.in = 1, l.in = 1/kappa)
