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
library(spNNGP)
source("Result_Manuscript/GraphAesthetics.R")

## data generation
nlist = c(200, 500, 1000)
Xlist = list(length = length(nlist))
Zlist = list(length = length(nlist))
f0_2D = function(x, y){return(x^2 + sqrt(abs(y-0.5)) + sin(8*x))}
Xlist = Zlist = list(length = length(nlist))
for(a in 1:length(nlist)){
   set.seed(a)
   n = nlist[a]
   # 2D data generation
   Xlist[[a]] = matrix(runif(2*n), n)
   Zlist[[a]] = f0_2D(Xlist[[a]][, 1], Xlist[[a]][, 2]) + rnorm(n) * 0.1
}

## setting for the sampling
target = 500; brn = 0
brnin = 0
Nk = c(4, 6, 8, 10, 12)
const = function(x){return(1)}
# singleton vector for kappak and tausqk
kappak = 3
tausqk = 3
beta = 2
N.pr = kappa.pr = tausq.pr = const

gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize + 1))

# import time from SPDE
source("2D/SPDE/functions_SPDE_sampling_2D.R")
source("2D/SPDE/functions_SPDE_2D.R")

# import time from GPI

source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")


# setting for NNGP
gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

starting <- list("phi" = 1/kappak[1], "sigma.sq" = tausqk[1], "tau.sq"=0.01, "nu" = 1)
tuning <- list("phi"= 0, "sigma.sq"= 0, "tau.sq"= 0, "nu" = 0)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1), "nu.unif" = c(1, 1.001))
cov.model <- "matern"
n.report <- 10 # any small integer

################## Splitting the one time / iterative calculation ############################
########## for the fair comparison, fix the parameter throughout the calculation #############

## GPI

onetime.GPI = microbenchmark(
   result.onetime.1 = sample.RJESS2D.onetime(Z = Zlist[[1]], X = Xlist[[1]], Nk = Nk, N.pr = N.pr, 
                                          kappak = kappak, kappa.pr = kappa.pr, 
                                          tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                                          mcmc = target, brn = 0, brn.ESS = brn.ESS),
   result.onetime.2 = sample.RJESS2D.onetime(Z = Zlist[[2]], X = Xlist[[2]], Nk = Nk, N.pr = N.pr, 
                                             kappak = kappak, kappa.pr = kappa.pr, 
                                             tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                                             mcmc = target, brn = 0, brn.ESS = brn.ESS),
   result.onetime.3 = sample.RJESS2D.onetime(Z = Zlist[[3]], X = Xlist[[3]], Nk = Nk, N.pr = N.pr, 
                                             kappak = kappak, kappa.pr = kappa.pr, 
                                             tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                                             mcmc = target, brn = 0, brn.ESS = brn.ESS),
   times = 10
)

# save variables to perform time comparison in the iterative step 
result.onetime.GPI.1 = sample.RJESS2D.onetime(Z = Zlist[[1]], X = Xlist[[1]], Nk = Nk, N.pr = N.pr, 
                                          kappak = kappak, kappa.pr = kappa.pr, 
                                          tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                                          mcmc = target, brn = 0, brn.ESS = brn.ESS)
result.onetime.GPI.2 = sample.RJESS2D.onetime(Z = Zlist[[2]], X = Xlist[[2]], Nk = Nk, N.pr = N.pr, 
                                          kappak = kappak, kappa.pr = kappa.pr, 
                                          tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                                          mcmc = target, brn = 0, brn.ESS = brn.ESS)
result.onetime.GPI.3 = sample.RJESS2D.onetime(Z = Zlist[[3]], X = Xlist[[3]], Nk = Nk, N.pr = N.pr, 
                                          kappak = kappak, kappa.pr = kappa.pr, 
                                          tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                                          mcmc = target, brn = 0, brn.ESS = brn.ESS)

iter.GPI = microbenchmark(
   result.iter.1 = sample.RJESS2D.iter(Z = Zlist[[1]], X = Xlist[[1]], N.pr = function(x){return(1)}, Nk = Nk, 
                                       result_list = result.onetime.GPI.1$result_list, 
                                       param_index_list = result.onetime.GPI.1$param_index_list, brn.ESS = brn.ESS),
   result.iter.2 = sample.RJESS2D.iter(Z = Zlist[[2]], X = Xlist[[2]], N.pr = function(x){return(1)}, Nk = Nk, 
                                       result_list = result.onetime.GPI.2$result_list, 
                                       param_index_list = result.onetime.GPI.2$param_index_list, brn.ESS = brn.ESS),
   result.iter.3 = sample.RJESS2D.iter(Z = Zlist[[3]], X = Xlist[[3]], N.pr = function(x){return(1)}, Nk = Nk, 
                                       result_list = result.onetime.GPI.3$result_list, 
                                       param_index_list = result.onetime.GPI.3$param_index_list, brn.ESS = brn.ESS),
   times = 10
)

############################### SPDE #####################################
onetime.SPDE = microbenchmark(
   result.onetime.1 = sample.exact2D.seq(X = Xlist[[1]], Z = Zlist[[1]], sigsq = 0.1^2,
                               Nk = Nk, N.pr = N.pr, 
                               kappak = kappak, kappa.pr = kappa.pr,
                               tausqk = tausqk, tausq.pr = tausq.pr,
                               beta = 2, mcmc = target, brn = brn, seed = 1234),
   result.onetime.2 = sample.exact2D.seq(X = Xlist[[2]], Z = Zlist[[2]], sigsq = 0.1^2,
                                         Nk = Nk, N.pr = N.pr, 
                                         kappak = kappak, kappa.pr = kappa.pr,
                                         tausqk = tausqk, tausq.pr = tausq.pr,
                                         beta = 2, mcmc = target, brn = brn, seed = 1234),
   result.onetime.3 = sample.exact2D.seq(X = Xlist[[3]], Z = Zlist[[3]], sigsq = 0.1^2,
                                         Nk = Nk, N.pr = N.pr, 
                                         kappak = kappak, kappa.pr = kappa.pr,
                                         tausqk = tausqk, tausq.pr = tausq.pr,
                                         beta = 2, mcmc = target, brn = brn, seed = 1234),
   times = 10
)

result.onetime.SPDE.1 = sample.exact.onetime(X = Xlist[[1]], Z = Zlist[[1]], sigsq = 0.1^2,
                                      Nk = Nk, N.pr = N.pr, 
                                      kappak = kappak, kappa.pr = kappa.pr,
                                      tausqk = tausqk, tausq.pr = tausq.pr,
                                      beta = 2, mcmc = target, brn = brn, seed = 1234)
result.onetime.SPDE.2 = sample.exact.onetime(X = Xlist[[2]], Z = Zlist[[2]], sigsq = 0.1^2,
                                      Nk = Nk, N.pr = N.pr, 
                                      kappak = kappak, kappa.pr = kappa.pr,
                                      tausqk = tausqk, tausq.pr = tausq.pr,
                                      beta = 2, mcmc = target, brn = brn, seed = 1234)
result.onetime.SPDE.3 = sample.exact.onetime(X = Xlist[[3]], Z = Zlist[[3]], sigsq = 0.1^2,
                                      Nk = Nk, N.pr = N.pr, 
                                      kappak = kappak, kappa.pr = kappa.pr,
                                      tausqk = tausqk, tausq.pr = tausq.pr,
                                      beta = 2, mcmc = target, brn = brn, seed = 1234)

iter.SPDE = microbenchmark(
   result.iter.1 = sample.exact.iter(Nk = Nk, N.pr = N.pr, 
                                     kappak = kappak, kappa.pr = kappa.pr,
                                     tausqk = tausqk, tausq.pr = tausq.pr,
                                     beta = 2, mcmc, brn, sigsq = 0.01, 
                                     param_index_list = result.onetime.SPDE.1$param_index_list, 
                                     chol_prec_grid = result.onetime.SPDE.1$chol_prec_grid, 
                                     mean_grid = result.onetime.SPDE.1$mean_grid, seed = 1234),
   result.iter.2 = sample.exact.iter(Nk = Nk, N.pr = N.pr, 
                                     kappak = kappak, kappa.pr = kappa.pr,
                                     tausqk = tausqk, tausq.pr = tausq.pr,
                                     beta = 2, mcmc, brn, sigsq = 0.01, 
                                     param_index_list = result.onetime.SPDE.2$param_index_list, 
                                     chol_prec_grid = result.onetime.SPDE.2$chol_prec_grid, 
                                     mean_grid = result.onetime.SPDE.2$mean_grid, seed = 1234),
   result.iter.3 = sample.exact.iter(Nk = Nk, N.pr = N.pr, 
                                     kappak = kappak, kappa.pr = kappa.pr,
                                     tausqk = tausqk, tausq.pr = tausq.pr,
                                     beta = 2, mcmc, brn, sigsq = 0.01, 
                                     param_index_list = result.onetime.SPDE.3$param_index_list, 
                                     chol_prec_grid = result.onetime.SPDE.3$chol_prec_grid, 
                                     mean_grid = result.onetime.SPDE.3$mean_grid, seed = 1234),
   times = 10
)

starting <- list("phi" = 1/kappak[1], "sigma.sq" = tausqk[1], "tau.sq"=0.01, "nu" = 1)
tuning <- list("phi"= 0, "sigma.sq"= 0, "tau.sq"= 0, "nu" = 0)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1), "nu.unif" = c(1, 1.001))
cov.model <- "matern"
n.report <- 10 #any small integer
iter.NNGP = microbenchmark(
   result1 = {Z = Zlist[[1]]; X = as.matrix(Xlist[[1]]); m.r <- spNNGP(Z ~ X-1, coords=X, starting=starting, method="response", n.neighbors=10,
                                                                       tuning=tuning, priors=priors, cov.model=cov.model,
                                                                       n.samples=target, n.omp.threads=1, n.report=n.report); p.r <- predict(m.r, X.0 = gridmat, coords.0 = gridmat, n.omp.threads=1); p.r$p.y.0},
   result2 = {Z = Zlist[[2]]; X = as.matrix(Xlist[[2]]); m.r <- spNNGP(Z ~ X-1, coords=X, starting=starting, method="response", n.neighbors=10,
                                                                       tuning=tuning, priors=priors, cov.model=cov.model,
                                                                       n.samples=target, n.omp.threads=1, n.report=n.report); p.r <- predict(m.r, X.0 = gridmat, coords.0 = gridmat, n.omp.threads=1); p.r$p.y.0},
   result3 = {Z = Zlist[[3]]; X = as.matrix(Xlist[[3]]); m.r <- spNNGP(Z ~ X-1, coords=X, starting=starting, method="response", n.neighbors=10,
                                                                       tuning=tuning, priors=priors, cov.model=cov.model,
                                                                       n.samples=target, n.omp.threads=1, n.report=n.report); p.r <- predict(m.r, X.0 = gridmat, coords.0 = gridmat, n.omp.threads=1);p.r$p.y.0},
   times = 10
)

iter.NNGP <- iter.NNGP$time / targets
########### save time for one time ###########

time_comparison_onetime_unify = rbind(onetime.SPDE, onetime.GPI)
time_comparison_onetime_unify <- time_comparison_onetime_unify %>%
   mutate(n = factor(expr)) %>%
   mutate(log_time = log(time)) %>%
   mutate(method = c(rep("SPDE", nrow(onetime.SPDE)), rep("GPI", nrow(onetime.GPI)))) %>%
   mutate(method = factor(method))
levels(time_comparison_unify$n) <- nlist
time_comparison_unify <- time_comparison_unify %>% mutate(n = as.numeric(as.character(n)))


filename = paste0("Result_Manuscript/Time_dataframe/time_2D_onetime_", target, ".RData")
save(time_comparison_onetime_unify, file = filename)

########### save time for iteration ###########

time_comparison_iter_unify = rbind(iter.SPDE, iter.GPI, iter.NNGP)
time_comparison_iter_unify <- time_comparison_iter_unify %>%
   mutate(n = factor(expr)) %>%
   mutate(log_time = log(time)) %>%
   mutate(method = c(rep("SPDE", nrow(iter.SPDE)), rep("GPI", nrow(iter.GPI)), rep("NNGP", nrow(iter.NNGP)))) %>%
   mutate(method = factor(method))
levels(time_comparison_unify$n) <- nlist
time_comparison_unify <- time_comparison_unify %>% mutate(n = as.numeric(as.character(n)))


filename = paste0("Result_Manuscript/Time_dataframe/time_2D_iter_", target, ".RData")
save(time_comparison_iter_unify, file = filename)
