rm(list = ls())
gc()
library(Matrix)
library(grDevices)
library(ggplot2)
library(rSPDE)
library(doParallel)
source("2D/SPDE/functions_SPDE_sampling_2D.R")
source("2D/SPDE/functions_SPDE_2D.R")

kappa = 2
brnin = 0
target = 25

const = function(x){
   return(1)
}
dpoi5 = function(x){
   return(dpois(x, lambda = 5))
}
f0 = function(x, y){
   return(sin(5*x + 2*y) + 2*y^2)
}
Nk = c(3, 5, 8, 10, 15)
gridsize = 40
M = 50
nlist = c(200, 500, 1000)
# gridmat is a (gridsize^2) by 2 matrix!
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))
MSE_list = matrix(nrow = M, ncol = length(nlist))

#################### Parallel computing ###########################
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)
############################################################
for(a in 1:length(nlist)){
   n = nlist[a]
   filename = paste0("Result_Manuscript/obs_n2D", n, ".RData")
   load(filename)
   output <- foreach (m = 1:M, .packages = c("Matrix", "rSPDE")) %dopar% {
      X = df[((m-1)*n+1):(m*n), c(1, 2)]
      Z = df$Z[((m-1)*n+1):(m*n)]
      result = sample.exact2D.seq(X, Z, sigsq = 0.1^2, # N.pr = function(x){return(1)},
                                  N.pr = dpoi5,
                                  Nk = Nk, kappa.init = kappa, mcmc = target, brn = brnin, seed = 1234)
      g_list = result$g_list
      y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0, alpha1 = 0.9, alpha2 = 0.95)
      print(m)
      mean((y.plot$truefun - y.plot$mean)^2)
   }
   MSE_list[, a] = simplify(output)
}

## Plot 4. For the replicated result ##
library(tidyverse)
MSE.df = data.frame(MSE_list)
colnames(MSE.df) = nlist
MSE.df.plot = gather(MSE.df, key = "n", value = "MSE")
MSE.df.plot$n <- factor(MSE.df.plot$n, levels = nlist)
ggplot(MSE.df.plot) + 
   geom_boxplot(aes(x = n, y = MSE))
## save as Rdata
save(MSE.df.plot, file = "Result_Manuscript/MSE_dataframe/MSE_SPDE_2D.RData")
