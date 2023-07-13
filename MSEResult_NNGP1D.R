rm(list = ls())
gc()
# install.packages("spNNGP")
library(spNNGP)
library(doParallel)
library(tidyverse)
source("Result_Manuscript/GraphAesthetics.R")


############################# our example ####################################

kappa = 2
brnin = 0
target = 2500
sigsq = 0.01

f0_1D = function(x){return(1*x^2+sin(8*x))}
grid.plot = c(0:1000)/1000

M = 50
nlist = c(200, 500, 1000)
MSE_list = matrix(nrow = M, ncol = length(nlist))
# gridmat is a (gridsize^2) by 2 matrix!

starting <- list("phi" = 1/kappa, "sigma.sq" = 1, "tau.sq"=0.01, "nu" = 1)
tuning <- list("phi"= 0, "sigma.sq"= 0, "tau.sq"= 0, "nu" = 0)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1), "nu.unif" = c(1, 1.001))
cov.model <- "matern"
n.report <- 10 #any small integer

#################### Parallel computing ###########################
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers/2)
registerDoParallel(cl)
###################################################################

for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   filename = paste0("Result_Manuscript/obs_n1D", n, ".RData")
   load(filename)
   for(m in 1:M){  
      # m th dataset among M = 50 dataset
      X = as.matrix(cbind(df$X[((m-1)*n+1):(m*n)], rep(0, n)))
      Y = as.matrix(df$Z[((m-1)*n+1):(m*n)])
      m.r <- spNNGP(Y ~ X-1, coords=X, starting=starting, method="response", n.neighbors=10,
                    tuning=tuning, priors=priors, cov.model=cov.model,
                    n.samples=target, n.omp.threads=1, n.report=n.report)
      p.r <- predict(m.r, X.0 = grid.plot, coords.0 = grid.plot, n.omp.threads=1)
      pred.grid <- p.r$p.y.0
      true.grid <- f0(grid.plot)
      mean.grid <- apply(pred.grid, 1, mean)
      mean((true.grid - mean.grid)^2)
   }
}

##################################################################################

## Plot 4. For the replicated result ##
MSE.df = data.frame(MSE_list)
colnames(MSE.df) = nlist
MSE.df.plot = gather(MSE.df, key = "n", value = "MSE")
MSE.df.plot$n <- factor(MSE.df.plot$n, levels = c(200, 500, 1000))

## save as Rdata
save(MSE.df.plot, file = "Result_Manuscript/MSE_dataframe/MSE_NNGP_1D.RData")