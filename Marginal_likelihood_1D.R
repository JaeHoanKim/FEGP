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


target = 0
brn = 0
brn.ESS = 0
# setting for the Matern parameters
kappak = seq(1, 5, 0.5)
tausqk = 1
Nk = c(4, 6, 8, 10, 14, 18)
const = function(x){
   return(1)
}

kappa.pr = tausq.pr = N.pr = const 
beta = 4



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

prob_N_plot <- data.frame(prob_N = unlist(prob_N_GPI), n = rep(nlist, each = length(Nk)), 
                          N = rep(Nk, times = length(nlist))) %>%
   mutate(n = factor(n))

plot_N_GPI <- ggplot(prob_N_plot, aes(x = N, y = prob_N)) + geom_line(aes(color = n),  linewidth = 1) + 
   geom_point() +
   labs(y = "p(N)") + theme1

prob_kappa_plot <- data.frame(prob_kappa = unlist(prob_kappa_GPI), n = rep(nlist, each = length(kappak)), 
                              kappa = rep(kappak, times = length(nlist))) %>%
   mutate(n = factor(n))

plot_kappa_GPI <- ggplot(prob_kappa_plot, aes(x = kappa, y = prob_kappa)) + geom_line(aes(color = n),  linewidth = 1) + 
   geom_point() +
   labs(y = "p(kappa)") + theme1


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

prob_N_plot <- data.frame(prob_N = unlist(prob_N_SPDE), n = rep(nlist, each = length(Nk)), 
                        N = rep(Nk, times = length(nlist))) %>%
   mutate(n = factor(n))

plot_N_SPDE <- ggplot(prob_N_plot, aes(x = N, y = prob_N)) + geom_line(aes(color = n),  linewidth = 1) + 
   geom_point() +
   labs(y = "p(N)") + theme1
   
prob_kappa_plot <- data.frame(prob_kappa = unlist(prob_kappa_SPDE), n = rep(nlist, each = length(kappak)), 
                          kappa = rep(kappak, times = length(nlist))) %>%
   mutate(n = factor(n))
plot_kappa_SPDE <- ggplot(prob_kappa_plot, aes(x = kappa, y = prob_kappa)) + geom_line(aes(color = n),  linewidth = 1) + 
   geom_point() +
   labs(y = "p(kappa)") + theme1

fileloc = "Result_Manuscript/marginal_likelihood/"
library(ggpubr)
plot_marg = ggarrange(plotlist = list(plot_N_GPI, plot_kappa_GPI,
                                      plot_N_SPDE, plot_kappa_SPDE), nrow = 2, ncol = 2)

pdf(file = paste0(fileloc, "marginal.pdf"), width = 10, height = 8)
print(plot_marg)
dev.off()

####################################
## For The 2D Marginal likelihood ##
####################################

f0_2D = function(x, y){return(sin(5*x + 2*y) + 2*y^2)}
M = 50
nlist = c(200, 500, 1000)
df_2D = list(length = length(nlist))

for(i in 1:length(nlist)){
   set.seed(i)
   n = nlist[i]
   # 2D data generation
   X = matrix(runif(2*n*M), n*M)
   Z = f0_2D(X[, 1], X[, 2]) + rnorm(n*M) * 0.1
   df_2D[[i]] = data.frame(X, Z)
}

# setting for the Matern parameters and sampling

beta = 4
d = 2
nu = beta - d/2

const = function(x){
   return(1)
}
Nk = c(6, 8, 10, 14, 18, 22, 26, 30)
N.pr = kappa.pr = tausq.pr = const
tausq.pr = function(x){return(invgamma::dinvgamma(x, 1, 1))}
kappa.pr = function(x){return(1/x^2)}
kappak = seq(1, 5, 1)
tausqk = 1

target = 0
brn = 0
brn.ESS = 0
gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))
m = 1
source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")

for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_2D[[a]]
   X = df[((m-1)*n+1):(m*n), c(1, 2)]
   Z = df$Z[((m-1)*n+1):(m*n)]
   result = sample.RJESS2D.seq(Z = Z, X = X, Nk = Nk, N.pr = N.pr, 
                               kappak = kappak, kappa.pr = kappa.pr, 
                               tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                               mcmc = target, brn = 0, brn.ESS = brn.ESS)
   log_prob_list_GPI[[a]] = result$log_prob_N_list
   log_prob_mat = matrix(log_prob_list_GPI[[a]], nrow = length(kappak))
   log_prob_mat = log_prob_mat - max(log_prob_mat)
   prob_mat = exp(log_prob_mat)
   prob_N_GPI[[a]] = colSums(prob_mat) / sum(prob_mat)
   prob_kappa_GPI[[a]] = colSums(t(prob_mat)) / sum(prob_mat)
   
}
source("2D/SPDE/functions_SPDE_sampling_2D.R")
source("2D/SPDE/functions_SPDE_2D.R")

for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_2D[[a]]
   X = df[((m-1)*n+1):(m*n), c(1, 2)]
   Z = df$Z[((m-1)*n+1):(m*n)]
   result = sample.exact2D.seq(X, Z, sigsq = 0.1^2,
                               Nk = Nk, N.pr = N.pr, 
                               kappak = kappak, kappa.pr = kappa.pr,
                               tausqk = tausqk, tausq.pr = tausq.pr,
                               beta = 2, mcmc = target, brn = brn, seed = 1234)
   log_prob_list_SPDE[[a]] = result$log_prob_N_list
   # calculate the marginal likelihood of N
   log_prob_mat = matrix(log_prob_list_SPDE[[a]], nrow = length(kappak))
   log_prob_mat = log_prob_mat - max(log_prob_mat)
   prob_mat = exp(log_prob_mat)
   prob_N_SPDE[[a]] = colSums(prob_mat) / sum(prob_mat)
   prob_kappa_SPDE[[a]] = colSums(t(prob_mat)) / sum(prob_mat)
}

## Plots


prob_N_plot <- data.frame(prob_N = unlist(prob_N_GPI), n = rep(nlist, each = length(Nk)), 
                          N = rep(Nk, times = length(nlist))) %>%
   mutate(n = factor(n))

plot_N_GPI <- ggplot(prob_N_plot, aes(x = N, y = prob_N)) + geom_line(aes(color = n),  linewidth = 1) + 
   geom_point() +
   labs(y = "p(N)") + theme1

prob_kappa_plot <- data.frame(prob_kappa = unlist(prob_kappa_GPI), n = rep(nlist, each = length(kappak)), 
                              kappa = rep(kappak, times = length(nlist))) %>%
   mutate(n = factor(n))

plot_kappa_GPI <- ggplot(prob_kappa_plot, aes(x = kappa, y = prob_kappa)) + geom_line(aes(color = n),  linewidth = 1) + 
   geom_point() +
   labs(y = "p(kappa)") + theme1


prob_N_plot <- data.frame(prob_N = unlist(prob_N_SPDE), n = rep(nlist, each = length(Nk)), 
                          N = rep(Nk, times = length(nlist))) %>%
   mutate(n = factor(n))

plot_N_SPDE <- ggplot(prob_N_plot, aes(x = N, y = prob_N)) + geom_line(aes(color = n),  linewidth = 1) + 
   geom_point() +
   labs(y = "p(N)") + theme1

prob_kappa_plot <- data.frame(prob_kappa = unlist(prob_kappa_SPDE), n = rep(nlist, each = length(kappak)), 
                              kappa = rep(kappak, times = length(nlist))) %>%
   mutate(n = factor(n))
plot_kappa_SPDE <- ggplot(prob_kappa_plot, aes(x = kappa, y = prob_kappa)) + geom_line(aes(color = n),  linewidth = 1) + 
   geom_point() +
   labs(y = "p(kappa)") + theme1

library(ggpubr)
plot_marg = ggarrange(plotlist = list(plot_N_GPI, plot_kappa_GPI,
                                      plot_N_SPDE, plot_kappa_SPDE), nrow = 2, ncol = 2)
plot_marg
