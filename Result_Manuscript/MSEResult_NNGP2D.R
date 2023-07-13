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

f0 = function(x, y){
   return(sin(5*x + 2*y) + 2*y^2)
}

gridsize = 40
M = 50
nlist = c(200, 500, 1000)
MSE_list = matrix(nrow = M, ncol = length(nlist))
# gridmat is a (gridsize^2) by 2 matrix!
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

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
   n = nlist[a]
   filename = paste0("Result_Manuscript/obs_n2D", n, ".RData")
   load(filename)
   output <- foreach (m = 1:M, .packages = c("Matrix", "rSPDE", "spNNGP")) %dopar% {
      X = as.matrix(df[((m-1)*n+1):(m*n), c(1, 2)])
      Z = as.matrix(df$Z[((m-1)*n+1):(m*n)])
      ## Response
      m.r <- spNNGP(Z ~ X-1, coords=X, starting=starting, method="response", n.neighbors=10,
                    tuning=tuning, priors=priors, cov.model=cov.model,
                    n.samples=target, n.omp.threads=1, n.report=n.report)
      p.r <- predict(m.r, X.0 = gridmat, coords.0 = gridmat, n.omp.threads=1)
      pred.grid <- p.r$p.y.0
      true.grid <- f0(gridmat[, 1], gridmat[, 2])
      mean.grid <- apply(pred.grid, 1, mean)
      mean((true.grid - mean.grid)^2)
   }
   MSE_list[, a] = purrr::simplify(output)
}

#########################################
library(tidyverse)
MSE.df = data.frame(MSE_list)
colnames(MSE.df) = nlist
MSE.df.plot = gather(MSE.df, key = "n", value = "MSE")
MSE.df.plot$n <- factor(MSE.df.plot$n, levels = nlist)
ggplot(MSE.df.plot) + 
   geom_boxplot(aes(x = n, y = MSE))
## save as Rdata
save(MSE.df.plot, file = "Result_Manuscript/MSE_dataframe/MSE_NNGP_2D.RData")