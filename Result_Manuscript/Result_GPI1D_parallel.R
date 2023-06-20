rm(list = ls())
gc()
library(tidyverse)
library(fields)
library(FastGP)
library(ggpubr)
library(Rcpp)
library(ggplot2)
library(Matrix)
library(doParallel)
library(foreach)
sourceCpp("1D/GPI/inv_chol.cpp")
source("1D/GPI/functions_GPI.R")
source("1D/GPI/functions_GPI_sampling.R")

target = 2500
brn.ESS = 100
kappa = 2
const = function(x){
   return(1)
}

f0_1D = function(x){return(1*x^2+sin(8*x))}
alpha1 = 0.9 # coverage probability 1 (dark gray area)
alpha2 = 0.95 # coverage probability 2 (light gray area)
Nk = c(4, 6, 8, 10, 12)

M = 50
nlist = c(200, 500, 1000)
MSE_list = matrix(nrow = M, ncol = length(nlist))


### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)
##########################
for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   filename = paste0("Result_Manuscript/obs_n1D", n, ".RData")
   load(filename)
   output <- foreach (m = 1:M, .packages = c("Matrix", "fields", "FastGP")) %dopar% {  
      # m th dataset among M = 50 dataset
      X = df$X[((m-1)*n+1):(m*n)]
      Y = df$Z[((m-1)*n+1):(m*n)]
      result = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                              N.pr = const,
                              mcmc = target, brn=0, brn.ESS = brn.ESS)
      g.plot = tail(result$g_list, target) # choosing last `target` samples
      grid.plot = c(0:1000)/1000
      obs = data.frame(X, Y)
      y.plot = glist_to_plotdf(g.plot, grid.plot, true = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
      print(m)
      mean((y.plot$true - y.plot$mean)^2)
   }
   MSE_list[, a] = simplify(output)
}

## Plot 3. Fitting of the samples (credible intervals)

################## converting predicted values into usable form ##################
##################################################################################

## Plot 4. For the replicated result ##
MSE.df = data.frame(MSE_list)
colnames(MSE.df) = nlist
MSE.df.plot = gather(MSE.df, key = "n", value = "MSE")
MSE.df.plot$n <- factor(MSE.df.plot$n, levels = c(200, 500, 1000))

## save as Rdata
save(MSE.df.plot, file = "Result_Manuscript/MSE_dataframe/MSE_GPI_1D.RData")