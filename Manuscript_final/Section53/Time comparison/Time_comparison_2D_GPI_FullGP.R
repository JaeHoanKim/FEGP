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
source("GraphAesthetics.R")

## data generation
nlist = c(200, 500, 1000, 2000)
Xlist = list(length = length(nlist))
Zlist = list(length = length(nlist))
f0_2D = function(x, y){return(sin(5*x + 2*y) + 2*y^2)}
Xlist = Zlist = list(length = length(nlist))
for(a in 1:length(nlist)){
   set.seed(a)
   n = nlist[a]
   # 2D data generation
   Xlist[[a]] = matrix(runif(2*n), n)
   Zlist[[a]] = f0_2D(Xlist[[a]][, 1], Xlist[[a]][, 2]) + rnorm(n) * 0.1
}

## setting for the sampling
target = 1; brn = 0
brnin = 0


Nk= c(6, 10, 14, 18)
N.pr = function(x){return(rep(1, length(x)))}
kappa.pr = function(x){return(dgamma(x, 5, 1/5))}
kappa.sampler = function(){rgamma(1, 5, 1/5)}
tausq.pr = function(x){return(dgamma(x, 1, 1))}
tausq.sampler = function(){return (1)}
tausq = tausq.sampler()

beta = 2
d = 2
nu = beta - d/2
sigma = 0.1

gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize + 1))

# import time from GPI

source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")
brn.ESS = 0
time_comparison_GPI = microbenchmark(
   result1 = {result = sample.GPI2D(Z = Zlist[[1]], X = Xlist[[1]], Nk = Nk, N.pr = N.pr, 
                                    kappa.pr = kappa.pr, kappa.sampler = kappa.sampler,
                                    tausq.pr = tausq.pr, tausq.sampler = tausq.sampler, sigsq = sigma^2, beta = beta,
                                    mcmc = target, brn = 0, brn.ESS = brn.ESS); g_list = result$g_list;
                                          y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)},
   result2 = {result = sample.GPI2D(Z = Zlist[[2]], X = Xlist[[2]], Nk = Nk, N.pr = N.pr, 
                                    kappa.pr = kappa.pr, kappa.sampler = kappa.sampler,
                                    tausq.pr = tausq.pr, tausq.sampler = tausq.sampler, sigsq = sigma^2, beta = beta,
                                    mcmc = target, brn = 0, brn.ESS = brn.ESS); g_list = result$g_list;
                                          y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)},
   result3 = {result = sample.GPI2D(Z = Zlist[[3]], X = Xlist[[3]], Nk = Nk, N.pr = N.pr, 
                                    kappa.pr = kappa.pr, kappa.sampler = kappa.sampler,
                                    tausq.pr = tausq.pr, tausq.sampler = tausq.sampler, sigsq = sigma^2, beta = beta,
                                    mcmc = target, brn = 0, brn.ESS = brn.ESS); g_list = result$g_list;
                                          y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)},
   result4 = {result = sample.GPI2D(Z = Zlist[[4]], X = Xlist[[4]], Nk = Nk, N.pr = N.pr, 
                                    kappa.pr = kappa.pr, kappa.sampler = kappa.sampler,
                                    tausq.pr = tausq.pr, tausq.sampler = tausq.sampler, sigsq = sigma^2, beta = beta,
                                    mcmc = target, brn = 0, brn.ESS = brn.ESS); g_list = result$g_list;
                                          y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)},
   times = 10
)

# import from FullGP
   
source("2D/FullGP/functions_FullGP.R")

time_comparison_FullGP = microbenchmark(
   result1 = {result = sample_fullGP(X = Xlist[[1]], Z = Zlist[[1]], sigma = sigma, nu = nu, 
                            kappa_pr = kappa.pr, kappa_sampler = kappa.sampler,
                            tausq_pr = tausq.pr, tausq_sampler = tausq.sampler, target = target)},
   result2 = {result = sample_fullGP(X = Xlist[[2]], Z = Zlist[[2]], sigma = sigma, nu = nu, 
                                     kappa_pr = kappa.pr, kappa_sampler = kappa.sampler,
                                     tausq_pr = tausq.pr, tausq_sampler = tausq.sampler, target = target)},
   result3 = {result = sample_fullGP(X = Xlist[[3]], Z = Zlist[[3]], sigma = sigma, nu = nu, 
                                     kappa_pr = kappa.pr, kappa_sampler = kappa.sampler,
                                     tausq_pr = tausq.pr, tausq_sampler = tausq.sampler, target = target)},
   result4 = {result = sample_fullGP(X = Xlist[[4]], Z = Zlist[[4]], sigma = sigma, nu = nu, 
                                     kappa_pr = kappa.pr, kappa_sampler = kappa.sampler,
                                     tausq_pr = tausq.pr, tausq_sampler = tausq.sampler, target = target)},
   times = 10
)

########### Time df unify ###########
time_comparison_unify = data.frame(rbind(as.matrix(time_comparison_GPI), as.matrix(time_comparison_FullGP)))
time_comparison_unify <- time_comparison_unify %>%
   mutate(n = factor(expr)) %>%
   mutate(log_time = log(as.numeric(time))) %>%
   mutate(method = c(rep("GPI", nrow(time_comparison_GPI)), 
                     rep("FullGP", nrow(time_comparison_FullGP)))) %>%
   mutate(method = factor(method))
levels(time_comparison_unify$n) <- nlist 
time_comparison_unify <- time_comparison_unify %>% mutate(n = as.numeric(as.character(n)))


filename = paste0("Manuscript_final/Section53/Time comparison/time_2D_GPI_FullGP_", target, ".RData")
save(time_comparison_unify, file = filename)


############## code for plots ################
# library(ggplot2)
# target = 2
# load("Result_Manuscript/Time_dataframe/time_2D_2.RData")
# time.plot.1 <- ggplot(time_comparison_unify) +
#    geom_boxplot(aes(x = factor(n), y = log_time, color = method)) +
#    labs(title = paste0(target, " samples"), x = "n", y = "log(time)") + theme1
# 
# target = 500
# load("Result_Manuscript/Time_dataframe/time_2D_500.RData")
# time.plot.2 <- ggplot(time_comparison_unify) +
#    geom_boxplot(aes(x = factor(n), y = log_time, color = method)) +
#    labs(title = paste0(target, " samples"), x = "n", y = "log(time)") + theme1
# 
# library(gridExtra)
# pdf(file = "Graphs/Time_whole_plot.pdf", width = 12, height = 4)
# grid.arrange(time.plot.1, time.plot.2, ncol = 2)
# dev.off()