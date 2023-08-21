rm(list = ls())
gc()
library(tidyverse)
library(fields)
library(FastGP)
library(Rcpp)
library(ggplot2)
library(Matrix)
library(doParallel)
library(foreach)
library(rSPDE)
library(doParallel)
library(tidyverse)
library(spNNGP)
library(invgamma)
sourceCpp("1D/GPI/inv_chol.cpp")

### 1. true function setting & data generation

# alpha = 0.5

f0_2D = function(x, y){return(x^2 + sqrt(abs(y-0.5)) + sin(8*x))}

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
l.in = 1/kappa

const = function(x){
   return(1)
}
Nk = c(4, 6, 8, 10, 12)
N.pr = kappa.pr = tausq.pr = const
tausq.pr = function(x){return(invgamma::dinvgamma(x, 1, 1))}
kappa.pr = function(x){return(1/x^2)}
kappak = seq(2, 5, 1)
tausqk = seq(2, 5, 1)

target = 2500
brn = 0
brn.ESS = 1000

### 3. MSE calculation - 2D

gridsize = 40
# gridmat is a (gridsize^2) by 2 matrix!
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

###################################################################

starting <- list("phi" = 1/kappak[1], "sigma.sq" =  tausqk[1], "tau.sq" = 0.01, "nu" = beta - 1)
tuning <- list("phi"= 0.1, "sigma.sq"= 0.1, "tau.sq"= 0, "nu" = 0)
# shape and scale parameter
priors <- list("phi.Unif"=c(1/kappak[length(kappak)], 1/kappak[1]), "sigma.sq.IG"=c(1, 1), "tau.sq.IG"=c(0.1, 0.1), "nu.unif" = c(beta - 1, beta - 1 + 1e-3))
cov.model <- "matern"
n.report <- 10 

MSE_list_NNGP2D = matrix(nrow = M, ncol = length(nlist))

### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)

for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_2D[[a]]
   output <- foreach (m = 1:M, .packages = c("Matrix", "spNNGP")) %dopar% {
      X = as.matrix(df[((m-1)*n+1):(m*n), c(1, 2)])
      Z = as.matrix(df$Z[((m-1)*n+1):(m*n)])
      ## Response
      m.r <- spNNGP(Z ~ X-1, coords=X, starting=starting, method="response", n.neighbors=10,
                    tuning=tuning, priors=priors, cov.model=cov.model,
                    n.samples=target, n.omp.threads=1, n.report=n.report)
      p.r <- predict(m.r, X.0 = gridmat, coords.0 = gridmat, n.omp.threads=1)
      pred.grid <- p.r$p.y.0
      true.grid <- f0_2D(gridmat[, 1], gridmat[, 2])
      mean.grid <- apply(pred.grid, 1, mean)
      mean((true.grid - mean.grid)^2)
   }
   MSE_list_NNGP2D[, a] = purrr::simplify(output)
}
stopCluster(cl)

MSE_list_GPI2D = matrix(nrow = M, ncol = length(nlist))

source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")

### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)

for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_2D[[a]]
   output <- foreach (m = 1:M, .packages = c("Matrix", "fields", "FastGP")) %dopar% {
      X = df[((m-1)*n+1):(m*n), c(1, 2)]
      Z = df$Z[((m-1)*n+1):(m*n)]
      result = sample.RJESS2D.seq(Z = Z, X = X, Nk = Nk, N.pr = N.pr, 
                                  kappak = kappak, kappa.pr = kappa.pr, 
                                  tausqk = tausqk, tausq.pr = tausq.pr, sigsq = 0.1^2, beta = beta,
                                  mcmc = target, brn = 0, brn.ESS = brn.ESS)
      g_list = result$g_list
      y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)
      print(m)
      mean((y.plot$truefun - y.plot$mean)^2)
   }
   MSE_list_GPI2D[, a] = purrr::simplify(output)
}

stopCluster(cl)

MSE_list_SPDE2D = matrix(nrow = M, ncol = length(nlist))

source("2D/SPDE/functions_SPDE_sampling_2D.R")
source("2D/SPDE/functions_SPDE_2D.R")

### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)

############################################################
for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_2D[[a]]
   output <- foreach (m = 1:M, .packages = c("Matrix", "rSPDE")) %dopar% {
      X = df[((m-1)*n+1):(m*n), c(1, 2)]
      Z = df$Z[((m-1)*n+1):(m*n)]
      result = sample.exact2D.seq(X, Z, sigsq = 0.1^2,
                                  Nk = Nk, N.pr = N.pr, 
                                  kappak = kappak, kappa.pr = kappa.pr,
                                  tausqk = tausqk, tausq.pr = tausq.pr,
                                  beta = 2, mcmc = target, brn = brn, seed = 1234)
      g_list = result$g_list
      y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)
      print(m)
      mean((y.plot$truefun - y.plot$mean)^2)
   }
   MSE_list_SPDE2D[, a] = purrr::simplify(output)
}
stopCluster(cl)

MSE_list_2D = rbind(MSE_list_GPI2D, MSE_list_SPDE2D, MSE_list_NNGP2D)
save(MSE_list_2D, file = "MSE_list_generated_data_2D.RData")