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
f0 = function(x, y){
   return(sin(5*x + 2*y) + 2*y^2)
}

gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize + 1))

# import time from SPDE
source("2D/SPDE/functions_SPDE_sampling_2D.R")
source("2D/SPDE/functions_SPDE_2D.R")
time_comparison_SPDE = microbenchmark(
   result1 = {result = sample.exact2D.seq(X = Xlist[[1]], Z = Zlist[[1]], sigsq = 0.1^2, # N.pr = function(x){return(1)},
                                N.pr = const,
                                Nk = Nk, kappa.init = kappa, mcmc = target, brn = brnin, seed = 1234); g_list = result$g_list;
                                y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0, alpha1 = 0.9, alpha2 = 0.95)},
   result2 = {result = sample.exact2D.seq(X = Xlist[[2]], Z = Zlist[[2]], sigsq = 0.1^2, # N.pr = function(x){return(1)},
                                N.pr = const,
                                Nk = Nk, kappa.init = kappa, mcmc = target, brn = brnin, seed = 1234); g_list = result$g_list;
                                y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0, alpha1 = 0.9, alpha2 = 0.95)},
   result3 = {result = sample.exact2D.seq(X = Xlist[[3]], Z = Zlist[[3]], sigsq = 0.1^2, # N.pr = function(x){return(1)},
                                N.pr = const,
                                Nk = Nk, kappa.init = kappa, mcmc = target, brn = brnin, seed = 1234); g_list = result$g_list;
                                y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0, alpha1 = 0.9, alpha2 = 0.95)},
   times = 10
)

# import time from GPI

source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")
brn.ESS = 100
time_comparison_GPI = microbenchmark(
   result1 = {result = sample.RJESS2D.seq(Z = Zlist[[1]], X = Xlist[[1]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                mcmc = target, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS); g_list = result$g_list;
                                y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0, alpha1 = 0.9, alpha2 = 0.95)},
   result2 = {result = sample.RJESS2D.seq(Z = Zlist[[2]], X = Xlist[[2]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                mcmc = target, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS); g_list = result$g_list;
                                y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0, alpha1 = 0.9, alpha2 = 0.95)},
   result3 = {result = sample.RJESS2D.seq(Z = Zlist[[3]], X = Xlist[[3]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                mcmc = target, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS); g_list = result$g_list;
                                y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0, alpha1 = 0.9, alpha2 = 0.95)},
   times = 10
)

# import time from NNGP
gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

starting <- list("phi" = 1/kappa, "sigma.sq" = 1, "tau.sq"=0.01, "nu" = 1)
tuning <- list("phi"= 0, "sigma.sq"= 0, "tau.sq"= 0, "nu" = 0)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1), "nu.unif" = c(1, 1.001))
cov.model <- "matern"
n.report <- 10 #any small integer
library(spNNGP)
time_comparison_NNGP = microbenchmark(
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


time_comparison_unify = rbind(time_comparison_SPDE, time_comparison_GPI, time_comparison_NNGP)
time_comparison_unify <- time_comparison_unify %>%
   mutate(n = factor(expr)) %>%
   mutate(log_time = log(time)) %>%
   mutate(method = c(rep("SPDE", nrow(time_comparison_SPDE)), rep("GPI", nrow(time_comparison_GPI)), 
                     rep("NNGP", nrow(time_comparison_NNGP)))) %>%
   mutate(method = factor(method))
levels(time_comparison_unify$n) <- nlist 
time_comparison_unify <- time_comparison_unify %>% mutate(n = as.numeric(as.character(n)))
   

filename = paste0("Result_Manuscript/Time_dataframe/time_2D_", target, ".RData")
save(time_comparison_unify, file = filename)


############## code for plots ################
library(ggplot2)
target = 2
load("Result_Manuscript/Time_dataframe/time_2D_2.Rdata")
time.plot.1 <- ggplot(time_comparison_unify) +
   geom_boxplot(aes(x = factor(n), y = log_time, color = method)) +
   labs(title = paste0(target, " samples"), x = "n", y = "log(time)") + theme1

target = 500
load("Result_Manuscript/Time_dataframe/time_2D_500.Rdata")
time.plot.2 <- ggplot(time_comparison_unify) +
   geom_boxplot(aes(x = factor(n), y = log_time, color = method)) +
   labs(title = paste0(target, " samples"), x = "n", y = "log(time)") + theme1

library(gridExtra)
pdf(file = "Graphs/Time_plot.pdf", width = 12, height = 4)
grid.arrange(time.plot.1, time.plot.2, ncol = 2)
dev.off()

################## Splitting the one time / iterative calculation ############################
########## for the fair comparison, fix the parameter throughout the calculation #############

## GPI

onetime_GPI = microbenchmark(
   result.onetime.1 = sample.RJESS2D.onetime(Z = Zlist[[1]], X = Xlist[[1]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                             mcmc = target, brn = 0, nu.in = 1, tausq = 1 ,l.in = 1/kappa),
   result.onetime.2 = sample.RJESS2D.onetime(Z = Zlist[[2]], X = Xlist[[2]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                             mcmc = target, brn = 0, nu.in = 1, tausq = 1 ,l.in = 1/kappa), 
   result.onetime.3 = sample.RJESS2D.onetime(Z = Zlist[[3]], X = Xlist[[3]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                             mcmc = target, brn = 0, nu.in = 1, tausq = 1 ,l.in = 1/kappa),
   times = 10
   
)

iter_GPI = microbenchmark(
   result.iter.1 = sample.RJESS2D.iter(Z = Zlist[[1]], X = Xlist[[1]], N.pr = function(x){return(1)}, Nk = Nk, 
                                       result_list = result.onetime.1$result_list, N_list = result.onetime.1$N_list, sigsq = 0.1^2, brn.ESS = brn.ESS),
   result.iter.2 = sample.RJESS2D.iter(Z = Zlist[[1]], X = Xlist[[1]], N.pr = function(x){return(1)}, Nk = Nk, 
                                       result_list = result.onetime.2$result_list, N_list = result.onetime.2$N_list, sigsq = 0.1^2, brn.ESS = brn.ESS),
   result.iter.3 = sample.RJESS2D.iter(Z = Zlist[[1]], X = Xlist[[1]], N.pr = function(x){return(1)}, Nk = Nk, 
                                       result_list = result.onetime.3$result_list, N_list = result.onetime.3$N_list, sigsq = 0.1^2, brn.ESS = brn.ESS),
   times = 10
)

## SPDE


