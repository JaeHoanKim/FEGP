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
sourceCpp("1D/GPI/inv_chol.cpp")
### 1. true function setting & data generation

f0_1D = function(x){return (x^2 + sin(x))}
f0_2D = function(x, y){return(x^2 + sqrt(abs(y-0.5)) + sin(8*x))}
# next try if it still preserves the pattern: abs(y-0.5)


M = 50
nlist = c(200, 500, 1000)
target = 2500
brn = 0
brn.ESS = 100
kappa = 2
const = function(x){
   return(1)
}
Nk = c(4, 6, 8, 10, 12)
grid.plot = c(0:1000)/1000
df_1D = list(length = length(nlist))
df_2D = list(length = length(nlist))

for(i in 1:length(nlist)){
   set.seed(i)
   n = nlist[i]
   # 1D data generation
   X = runif(n*M)
   Z = f0_1D(X) + rnorm(n*M) * 0.1
   df_1D[[i]] = data.frame(X, Z)
   # 2D data generation
   X = matrix(runif(2*n*M), n*M)
   Z = f0_2D(X[, 1], X[, 2]) + rnorm(n*M) * 0.1
   df_2D[[i]] = data.frame(X, Z)
}

### 2. MSE calculation - 1D

source("1D/GPI/functions_GPI.R")
source("1D/GPI/functions_GPI_sampling.R")
MSE_list_GPI1D = matrix(nrow = M, ncol = length(nlist))

### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)
##########################
for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   df = df_1D[[a]]
   output <- foreach (m = 1:M, .packages = c("Matrix", "fields", "FastGP")) %dopar% {  
      # m th dataset among M = 50 dataset
      X = df$X[((m-1)*n+1):(m*n)]
      Y = df$Z[((m-1)*n+1):(m*n)]
      result = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                              N.pr = const,
                              mcmc = target, brn=0, brn.ESS = brn.ESS)
      g.plot = tail(result$g_list, target) # choosing last `target` samples
      obs = data.frame(X, Y)
      y.plot = glist_to_plotdf(g.plot, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
      print(m)
      mean((y.plot$true - y.plot$mean)^2)
   }
   MSE_list_GPI1D[, a] = simplify(output)
}
stopCluster(cl)


source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")
MSE_list_SPDE1D = matrix(nrow = M, ncol = length(nlist))
grid.plot = c(0:1000)/1000

### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)

for(a in 1:length(nlist)){
   n = nlist[a]
   df = df_1D[[a]]
   output <- foreach (m = 1:M, .packages = c("Matrix", "fields", "FastGP")) %dopar% {
      X = df$X[((m-1)*n+1):(m*n)]
      Y = df$Z[((m-1)*n+1):(m*n)]
      result = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = const,
                                kappa.init = kappa, mcmc = target, brn=0, seed = 1234)
      g.plot = tail(result$g_list, target) # choosing last `target` samples
      obs = data.frame(X, Y)
      y.plot = glist_to_plotdf(g.plot, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
      mean((y.plot$true - y.plot$mean)^2)
   }
   MSE_list_SPDE1D[, a] = simplify(output)
}

stopCluster(cl)

MSE_list_1D = rbind(MSE_list_GPI1D, MSE_list_SPDE1D)
save(MSE_list_1D, file = "MSE_list_generated_data_1D.RData")


### 3. MSE calculation - 2D

target = 2500
brn = 0
brn.ESS = 1000
kappa = 2

gridsize = 40
# gridmat is a (gridsize^2) by 2 matrix!
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

###################################################################

starting <- list("phi" = 1/kappa, "sigma.sq" = 1, "tau.sq" = 0.01, "nu" = 1)
tuning <- list("phi"= 0, "sigma.sq"= 0, "tau.sq"= 0, "nu" = 0)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1), "nu.unif" = c(1, 1.001))
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
      result = sample.RJESS2D.seq(Z = Z, X = X, N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                  mcmc = target, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS)
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
      result = sample.exact2D.seq(X, Z, sigsq = 0.1^2, # N.pr = function(x){return(1)},
                                  N.pr = const,
                                  Nk = Nk, kappa.init = kappa, mcmc = target, brn = brn, seed = 1234)
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