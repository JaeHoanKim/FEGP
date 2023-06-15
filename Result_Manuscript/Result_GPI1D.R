rm(list = ls())
gc()
library(tidyverse)
library(ggpubr)
library(Rcpp)
library(ggplot2)
library(Matrix)
sourceCpp("1D/GPI/inv_chol.cpp")
source("1D/GPI/functions_GPI.R")
source("1D/GPI/functions_GPI_sampling.R")

n = 200 # the number of observed data; 200, 500, 1000
filename = paste0("Result_Manuscript/obs_n1D", n, ".RData")
load(filename)
m = 1 # m th dataset among M = 50 dataset
X = df$X[((m-1)*n+1):(m*n)]
Y = df$Z[((m-1)*n+1):(m*n)]
target = 250
brn.ESS = 100
kappa = 2
dpois5 = function(x){
   return(dpois(x, lambda = 5))
}

Nk = c(4, 6, 8, 10, 12)
result = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                        N.pr = dpois5,
                        mcmc = target, brn=0, brn.ESS = brn.ESS)


## Plot 3. Fitting of the samples (credible intervals)

################## converting predicted values into usable form ##################
##################################################################################

f0_1D = function(x){return(1*x^2+sin(8*x))}
alpha1 = 0.9 # coverage probability 1 (dark gray area)
alpha2 = 0.95 # coverage probability 2 (light gray area)
sample.num = target # number of samples used for the credible intervals
g.plot = tail(result$g_list, sample.num) # choosing last `sample.num` samples
grid.plot = c(0:1000)/1000
obs = data.frame(X, Y)

y.plot = glist_to_plotdf(g.plot, grid.plot, true = f0_1D, alpha1 = 0.95, alpha2 = 0.9)

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