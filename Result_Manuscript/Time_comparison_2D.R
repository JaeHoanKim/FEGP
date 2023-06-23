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
source("Result_Manuscript/GraphAesthetics.R")
m = 1
nlist = c(200, 500, 1000)
Xlist = list(length = length(nlist))
Zlist = list(length = length(nlist))
for(a in 1:length(nlist)){
   n = nlist[a]
   filename = paste0("Result_Manuscript/obs_n2D", n, ".RData")
   load(filename)
   Xlist[[a]] = df[((m-1)*n+1):(m*n), c(1, 2)]
   Zlist[[a]] = df$Z[((m-1)*n+1):(m*n)]
}
target = 500
brnin = 0
Nk = c(4, 6, 8, 10, 12)
const = function(x){return(1)}
kappa = 2
# import time from SPDE
source("2D/SPDE/functions_SPDE_sampling_2D.R")
source("2D/SPDE/functions_SPDE_2D.R")
time_comparison_SPDE = microbenchmark(
   result1 = sample.exact2D.seq(X = Xlist[[1]], Z = Zlist[[1]], sigsq = 0.1^2, # N.pr = function(x){return(1)},
                                N.pr = const,
                                Nk = Nk, kappa.init = kappa, mcmc = target, brn = brnin, seed = 1234),
   result2 = sample.exact2D.seq(X = Xlist[[2]], Z = Zlist[[2]], sigsq = 0.1^2, # N.pr = function(x){return(1)},
                                N.pr = const,
                                Nk = Nk, kappa.init = kappa, mcmc = target, brn = brnin, seed = 1234),
   result3 = sample.exact2D.seq(Xlist[[3]], Zlist[[3]], sigsq = 0.1^2, # N.pr = function(x){return(1)},
                                N.pr = const,
                                Nk = Nk, kappa.init = kappa, mcmc = target, brn = brnin, seed = 1234),
   times = 10
)

# import time from GPI

source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")
brn.ESS = 100
time_comparison_GPI = microbenchmark(
   result1 = sample.RJESS2D.seq(Z = Zlist[[1]], X = Xlist[[1]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                mcmc = target, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS),
   result2 = sample.RJESS2D.seq(Z = Zlist[[2]], X = Xlist[[2]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                mcmc = target, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS),
   result3 = sample.RJESS2D.seq(Z = Zlist[[3]], X = Xlist[[3]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                mcmc = target, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS),
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
