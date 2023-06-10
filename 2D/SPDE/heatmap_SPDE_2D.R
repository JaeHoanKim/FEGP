################################
### heatmap for SPDE 2D case ###
################################

rm(list = ls())
gc()
library(Matrix)
library(grDevices)
library(ggplot2)
library(rSPDE)
source("2D/SPDE/functions_SPDE_sampling_2D.R")
source("2D/SPDE/functions_SPDE_2D.R")

n = 1000 # the number of observed data
f0 = function(x, y){
   return(sin(5*x + 2*y) + 2*y^2)
} # true function
X = matrix(runif(2*n), n)
Z = f0(X[, 1], X[, 2]) + rnorm(n) * 0.1

algo = "RJexact"

# For the RJexact method, generate the Q matrix and \Phi^T \Phi matrix
N = 10
kappa = 2
Q = as.matrix(Q2D(N, kappa))
Phi = Phi_2D(X, N)
PhiTPhi = as.matrix(t(Phi) %*% Phi)
heatmap(Q, Colv= NA, Rowv = NA)
heatmap(PhiTPhi, Colv= NA, Rowv = NA)
