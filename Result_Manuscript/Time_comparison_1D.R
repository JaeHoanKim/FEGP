rm(list = ls())
gc()
library(grDevices)
library(ggplot2)
library(Matrix)
library(microbenchmark)
library(tidyverse)
library(fields)
library(FastGP)
library(ggpubr)
library(Rcpp)
sourceCpp("1D/GPI/inv_chol.cpp")
source("Result_Manuscript/GraphAesthetics.R")
m = 1
nlist = c(200, 500, 1000)
Xlist = list(length = length(nlist))
Ylist = list(length = length(nlist))
for(a in 1:length(nlist)){
   n = nlist[a]
   filename = paste0("Result_Manuscript/obs_n1D", n, ".RData")
   load(filename)
   Xlist[[a]] = df$X[((m-1)*n+1):(m*n)]
   Ylist[[a]] = df$Z[((m-1)*n+1):(m*n)]
}
target = 500
Nk = c(4, 6, 8, 10, 12)
const = function(x){return(1)}
kappa = 2

# import time from SPDE
source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")
time_comparison_SPDE = microbenchmark(
   result1 = sample.exact.seq(Xlist[[1]], Ylist[[1]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                              kappa.init = kappa, mcmc = target, brn=0, seed = 1234),
   result2 = sample.exact.seq(Xlist[[2]], Ylist[[2]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                              kappa.init = kappa, mcmc = target, brn=0, seed = 1234),
   result3 = sample.exact.seq(Xlist[[3]], Ylist[[3]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                              kappa.init = kappa, mcmc = target, brn=0, seed = 1234),
   times = 10
)

# import time from GPI

source("1D/GPI/functions_GPI.R")
source("1D/GPI/functions_GPI_sampling.R")
brn.ESS = 100
time_comparison_GPI = microbenchmark(
   result1 = sample.ESS.seq(Xlist[[1]], Ylist[[1]], sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = target, brn=0, brn.ESS = brn.ESS),
   result2 = sample.ESS.seq(Xlist[[2]], Ylist[[2]], sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = target, brn=0, brn.ESS = brn.ESS),
   result3 = sample.ESS.seq(Xlist[[3]], Ylist[[3]], sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = target, brn=0, brn.ESS = brn.ESS),
   times = 10
)

time_comparison_unify = rbind(time_comparison_SPDE, time_comparison_GPI)
time_comparison_unify <- time_comparison_unify %>%
   mutate(n = factor(expr)) %>%
   mutate(log_time = log(time)) %>%
   mutate(method = c(rep("SPDE", nrow(time_comparison_SPDE)), rep("GPI", nrow(time_comparison_GPI)))) %>%
   mutate(method = factor(method))
levels(time_comparison_unify$n) <- nlist 
time_comparison_unify <- time_comparison_unify %>% mutate(n = as.numeric(as.character(n)))
   
ggplot(time_comparison_unify) +
   geom_boxplot(aes(x = factor(n), y = log_time, color = method)) +
   labs(title = paste0("time comparison when generating ", target, " samples")) + theme1
