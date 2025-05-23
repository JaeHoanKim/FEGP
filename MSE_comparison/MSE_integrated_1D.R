rm(list = ls())
gc()
library(fields)
library(FastGP)
library(Rcpp)
library(ggplot2)
library(Matrix)
library(foreach)
library(rSPDE)
library(doParallel)
library(tidyverse)
# sourceCpp("1D/GPI/inv_chol.cpp")

### 1. true function setting & data generation

alpha = 1.4

f0_1D = function(x, trun = 500){
   value = 0
   for(j in 1:trun){
      value = value + sin(j) * cos(pi * (j - 1/2) * x) * j^(- alpha - 1)
   }
   return(value * sqrt(2))
}
# f0_1D = function(x){return(x^2 + sin(x))}

const = function(x){return(1)}


M = 50
nlist = c(200, 500, 1000)
df_1D = list(length = length(nlist))

for(i in 1:length(nlist)){
   set.seed(i)
   n = nlist[i]
   # 1D data generation
   X = runif(n*M)
   Z = f0_1D(X) + rnorm(n*M) * 0.1
   df_1D[[i]] = data.frame(X, Z)
   # 2D data generation
   X = matrix(runif(2*n*M), n*M)
}


target = 2500
brn = 0
brn.ESS = 1000
# setting for the Matern parameters
kappak = seq(1, 5, 0.5)
tausqk = 1
Nk = c(4, 6, 8, 10, 12)
kappa.pr = tausq.pr = N.pr = const 
beta = 4

const = function(x){
   return(1)
}

grid.plot = c(0:1000)/1000

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
      result = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr,
                              kappak = kappak, kappa.pr = kappa.pr, 
                              tausqk = tausqk, tausq.pr = tausq.pr, 
                              beta = beta, mcmc = target, brn=0, brn.ESS = brn.ESS)
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
      result = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr, beta = beta,
                                kappak = kappak, kappa.pr = kappa.pr,
                                tausqk = tausqk, tausq.pr = tausq.pr,
                                mcmc = target, brn=0, seed = 1234)
      g.plot = tail(result$g_list, target) # choosing last `target` samples
      obs = data.frame(X, Y)
      y.plot = glist_to_plotdf(g.plot, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
      mean((y.plot$true - y.plot$mean)^2)
   }
   MSE_list_SPDE1D[, a] = simplify(output)
}

stopCluster(cl)

MSE_list_1D = rbind(MSE_list_GPI1D, MSE_list_SPDE1D)
save(MSE_list_1D, file = "MSE_list_generated_data_1D_1.RData")

##################################################
####### Coverage plot for a specific data ########
##################################################

m = 1
cover.plot.GPI.list = list(length = length(nlist))
cover.plot.GPI.beta2.list = list(length = length(nlist))
cover.plot.SPDE.list = list(length = length(nlist))
cover.plot.SPDE.beta2.list = list(length = length(nlist))
cover.plot.SPDE.kappa.list = list(length = length(nlist)) # cover plot list for fixed kappa

## Save plots regarding SPDE

for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   df = df_1D[[a]]
   X = df$X[((m-1)*n+1):(m*n)]
   Y = df$Z[((m-1)*n+1):(m*n)]
   obs = data.frame(X, Y)
   # result for SPDE
   result.SPDE = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr, beta = beta,
                                  kappak = kappak, kappa.pr = kappa.pr,
                                  tausqk = tausqk, tausq.pr = tausq.pr,
                                  mcmc = target, brn=0, seed = 1234)
   g.plot.SPDE = tail(result.SPDE$g_list, target) # choosing last `target` samples
   y.plot.SPDE = glist_to_plotdf(g.plot.SPDE, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)

   # result for SPDE with beta = 2
   result.SPDE.beta2 = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr, beta = 2,
                                        kappak = kappak, kappa.pr = kappa.pr,
                                        tausqk = tausqk, tausq.pr = tausq.pr,
                                        mcmc = target, brn=0, seed = 1234)
   g.plot.SPDE.beta2 = tail(result.SPDE.beta2$g_list, target) # choosing last `target` samples
   y.plot.SPDE.beta2 = glist_to_plotdf(g.plot.SPDE.beta2, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)

   # result for SPDE with fixed kappa with beta = 4

   result.SPDE.kappa = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr, beta = beta,
                                        kappak = 3.2, kappa.pr = kappa.pr,
                                        tausqk = tausqk, tausq.pr = tausq.pr,
                                        mcmc = target, brn=0, seed = 1234)
   g.plot.SPDE.kappa = tail(result.SPDE.kappa$g_list, target) # choosing last `target` samples
   y.plot.SPDE.kappa = glist_to_plotdf(g.plot.SPDE.kappa, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)

   # Coverage plot

   cover.plot.SPDE.list[[a]] <- ggplot(y.plot.SPDE, aes(x = x)) +
      geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
      geom_line(aes(y=mean), colour="blue") +
      geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) +
      geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
      # geom_point(data = obs, aes(X, Y), size = 0.3) +
      labs(title = paste0("SPDE (n = ", n, ", beta = ", beta, ")"), x = "x", y = "y")+
      theme(plot.title = element_text(hjust = 0.5))

   cover.plot.SPDE.kappa.list[[a]] <- ggplot(y.plot.SPDE.kappa, aes(x = x)) +
      geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
      geom_line(aes(y=mean), colour="blue") +
      geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) +
      geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
      labs(title = paste0("SPDE (n = ", n, ", kappa = 3.2)"), x = "x", y = "y")+
      theme(plot.title = element_text(hjust = 0.5))

   cover.plot.SPDE.beta2.list[[a]] <- ggplot(y.plot.SPDE.beta2, aes(x = x)) +
      geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
      geom_line(aes(y=mean), colour="blue") +
      geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) +
      geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
      labs(title = paste0("SPDE (n = ", n, ", beta = 2)"), x = "x", y = "y")+
      theme(plot.title = element_text(hjust = 0.5))
}




## Save plots regarding GPI
for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   df = df_1D[[a]]
   X = df$X[((m-1)*n+1):(m*n)]
   Y = df$Z[((m-1)*n+1):(m*n)]
   obs = data.frame(X, Y)
   # result for GPI
   result.GPI = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr,
                           kappak = kappak, kappa.pr = kappa.pr,
                           tausqk = tausqk, tausq.pr = tausq.pr,
                           beta = beta, mcmc = target, brn=0, brn.ESS = brn.ESS)
   g.plot.GPI = tail(result.GPI$g_list, target) # choosing last `target` samples
   y.plot.GPI = glist_to_plotdf(g.plot.GPI, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)

   # result for GPI with beta = 2
   result.GPI.beta2 = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = N.pr,
                               kappak = kappak, kappa.pr = kappa.pr,
                               tausqk = tausqk, tausq.pr = tausq.pr,
                               beta = beta, mcmc = target, brn=0, brn.ESS = brn.ESS)
   g.plot.GPI.beta2 = tail(result.GPI.beta2$g_list, target) # choosing last `target` samples
   y.plot.GPI.beta2 = glist_to_plotdf(g.plot.GPI.beta2, grid.plot, truefun = f0_1D, alpha1 = 0.95, alpha2 = 0.9)

   # Coverage plot

   cover.plot.GPI.list[[a]] <- ggplot(y.plot.GPI, aes(x = x)) +
      geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
      geom_line(aes(y=mean), colour="blue") +
      geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) +
      geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
      labs(title = paste0("GPI (n = ", n, ", beta = ", beta, ")"), x = "x", y = "y")+
      theme(plot.title = element_text(hjust = 0.5))

   cover.plot.GPI.beta2.list[[a]] <- ggplot(y.plot.GPI.beta2, aes(x = x)) +
      geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
      geom_line(aes(y=mean), colour="blue") +
      geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) +
      geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
      labs(title = paste0("GPI (n = ", n, ", beta = 2)"), x = "x", y = "y")+
      theme(plot.title = element_text(hjust = 0.5))
}


## Pair the graphs
library(gridExtra)
pdf(file = "Graphs/coverage_plot_1.pdf", width = 12, height = 12)
grid.arrange(cover.plot.GPI.beta2.list[[1]], cover.plot.GPI.list[[1]], cover.plot.SPDE.beta2.list[[1]], cover.plot.SPDE.list[[1]],
             cover.plot.GPI.beta2.list[[2]], cover.plot.GPI.list[[2]], cover.plot.SPDE.beta2.list[[2]], cover.plot.SPDE.list[[2]],
             cover.plot.GPI.beta2.list[[3]], cover.plot.GPI.list[[3]], cover.plot.SPDE.beta2.list[[3]], cover.plot.SPDE.list[[3]], nrow = 4, as.table = FALSE)
dev.off()

pdf(file = "Graphs/coverage_plot_2.pdf", width = 12, height = 6)
grid.arrange(cover.plot.SPDE.list[[1]], cover.plot.SPDE.kappa.list[[1]],
             cover.plot.SPDE.list[[2]], cover.plot.SPDE.kappa.list[[2]],
             cover.plot.SPDE.list[[3]], cover.plot.SPDE.kappa.list[[3]], nrow = 2, as.table = FALSE)
dev.off()


## Coverage plot for presentation slides
pdf(file = "Graphs/coverage_plot_simple_1000.pdf", width = 10, height = 4)
grid.arrange(cover.plot.GPI.beta2.list[[3]], cover.plot.SPDE.list[[3]], nrow = 1)
dev.off()


