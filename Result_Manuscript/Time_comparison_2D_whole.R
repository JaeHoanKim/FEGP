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

## data generation
nlist = 100000
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
target = 2; brn = 0
brnin = 0
Nk = 10 # fix Nk to be a singleton vector
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
time_comparison_SPDE = microbenchmark(
   result1 = {result = sample.exact2D.seq(X = Xlist[[1]], Z = Zlist[[1]], sigsq = 0.1^2,
                                          Nk = Nk, N.pr = N.pr, 
                                          kappak = kappak, kappa.pr = kappa.pr,
                                          tausqk = tausqk, tausq.pr = tausq.pr,
                                          beta = 2, mcmc = target, brn = brn, seed = 1234); g_list = result$g_list;
                                          y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)},
   result2 = {result = sample.exact2D.seq(X = Xlist[[2]], Z = Zlist[[2]], sigsq = 0.1^2,
                                          Nk = Nk, N.pr = N.pr, 
                                          kappak = kappak, kappa.pr = kappa.pr,
                                          tausqk = tausqk, tausq.pr = tausq.pr,
                                          beta = 2, mcmc = target, brn = brn, seed = 1234); g_list = result$g_list;
                                          y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)},
   result3 = {result = sample.exact2D.seq(X = Xlist[[3]], Z = Zlist[[3]], sigsq = 0.1^2,
                                          Nk = Nk, N.pr = N.pr, 
                                          kappak = kappak, kappa.pr = kappa.pr,
                                          tausqk = tausqk, tausq.pr = tausq.pr,
                                          beta = 2, mcmc = target, brn = brn, seed = 1234); g_list = result$g_list;
                                          y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)},
   times = 10
)

# import time from GPI

source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")
brn.ESS = 100
time_comparison_GPI = microbenchmark(
   result1 = {result = sample.RJESS2D.seq(Z = Zlist[[1]], X = Xlist[[1]], Nk = Nk, N.pr = N.pr, 
                                          kappak = kappak, kappa.pr = kappa.pr, 
                                          tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                                          mcmc = target, brn = 0, brn.ESS = brn.ESS); g_list = result$g_list;
                                          y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)},
   result2 = {result = sample.RJESS2D.seq(Z = Zlist[[2]], X = Xlist[[2]], Nk = Nk, N.pr = N.pr, 
                                          kappak = kappak, kappa.pr = kappa.pr, 
                                          tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                                          mcmc = target, brn = 0, brn.ESS = brn.ESS); g_list = result$g_list;
                                          y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)},
   result3 = {result = sample.RJESS2D.seq(Z = Zlist[[3]], X = Xlist[[3]], Nk = Nk, N.pr = N.pr, 
                                          kappak = kappak, kappa.pr = kappa.pr, 
                                          tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                                          mcmc = target, brn = 0, brn.ESS = brn.ESS); g_list = result$g_list;
                                          y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)},
   times = 10
)

# import time from NNGP
gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

starting <- list("phi" = 1/kappak[1], "sigma.sq" = tausqk[1], "tau.sq"=0.01, "nu" = 1)
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