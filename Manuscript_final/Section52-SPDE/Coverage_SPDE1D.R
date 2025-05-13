## comparison between original GP and d-frGPloc SPDE
rm(list = ls())
gc()
library(gridExtra)
library(fields)
library(FastGP)
library(Rcpp)
library(ggplot2)
library(Matrix)
library(foreach)
library(rSPDE)
library(doParallel)
library(tidyverse)
# sourceCpp("1D/GPI/inv_chol.cpp")

source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")
source("1D/MaternGP1D.R")
### 1. true function setting & data generation
alpha_val = 7
f0_1D = function(x, alpha = alpha_val, trun = 500){
   # value = 0
   # for(j in 1:trun){
   #    value = value + sin(j) * cos(pi * (j - 1/2) * x) * j^(- alpha - 1)
   # }
   # return(value * sqrt(2))
   return (2^alpha * abs(0.5 - x)^alpha)
}
const = function(x){return(1)}

M = 50
nlist = c(200, 500, 1000)
df_1D = list(length = length(nlist))
for(i in 1:length(nlist)){
   set.seed(i)
   n = nlist[i]
   # 1D data generation
   X = runif(n)
   Z = f0_1D(X) + rnorm(n) * 0.1
   df_1D[[i]] = data.frame(X, Z)
}

target = 2500
brn = 0
brn.ESS = 100
# setting for the Matern parameters
kappak = 10
kappa_cand = seq(5, 10, length.out = 8)
tausqk = 1
Nk = c(4, 6, 8, 10, 12, 16, 20, 25, 30)
kappa.pr = tausq.pr = N.pr = const 
beta = 2
grid.plot = c(0:300)/300
##################################################
####### Coverage plot for a specific data ########
##################################################

m = 1
cover.plot.SPDE.list = list(length = length(nlist))
cover.plot.Matern.list = list(length = length(nlist))

## Save plots regarding SPDE

for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   df = df_1D[[a]]
   X = df$X
   Y = df$Z
   obs = data.frame(X, Y)
   result.SPDE = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr, beta = beta,
                                  kappak = kappak, kappa.pr = kappa.pr,
                                  tausqk = tausqk, tausq.pr = tausq.pr,
                                  mcmc = target, brn=0, seed = 1234)
   g.plot.SPDE = tail(result.SPDE$g_list, target) # choosing last `target` samples
   y.plot.SPDE = glist_to_plotdf(g.plot.SPDE, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
   
   l = 1/kappak[1]
   v = alpha_val
   tmp = sample_posterior(X, Y, grid.plot, l, v, sigma_n = 0.1, n_samples = 20)
   y.plot.Matern = matrix_to_plotdf(tmp, truefun = f0_1D, grid = grid.plot)
   
   ## Plots 
   y_min = min(min(y.plot.Matern), min(y.plot.SPDE))
   y_max = max(max(y.plot.Matern), max(y.plot.SPDE))
   
   cover.plot.SPDE.list[[a]] <- ggplot(y.plot.SPDE, aes(x = x)) +
      geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
      geom_line(aes(y=mean), colour="blue") +
      geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) +
      geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
      # geom_point(data = obs, aes(X, Y), size = 0.3) +
      labs(x = "x", y = "y")+
      scale_y_continuous(limits = c(y_min, y_max))+
      theme(plot.title = element_text(hjust = 0.5))
   
   cover.plot.Matern.list[[a]] <- ggplot(y.plot.Matern, aes(x = x)) +
      geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
      geom_line(aes(y=mean), colour="blue") +
      geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) +
      geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
      # geom_point(data = obs, aes(X, Y), size = 0.3) +
      labs(x = "x", y = "y")+
      scale_y_continuous(limits = c(y_min, y_max))+
      theme(plot.title = element_text(hjust = 0.5))
   print(a)
}

## Pair of the graphs

pdf(file = "Graphs/coverage_plot_SPDE_Matern.pdf", width = 12, height = 6)
grid.arrange(cover.plot.SPDE.list[[1]], cover.plot.Matern.list[[1]],
             cover.plot.SPDE.list[[2]], cover.plot.Matern.list[[2]],
             cover.plot.SPDE.list[[3]], cover.plot.Matern.list[[3]], nrow = 2, as.table = FALSE)
dev.off()
