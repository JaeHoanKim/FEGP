rm(list = ls())
gc()
library(grDevices)
library(ggplot2)
library(Matrix)
library(doParallel)
library(foreach)
library(rSPDE)
source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")

target = 2500
brn = 0
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
grid.plot = c(0:1000)/1000


### Parallel computing ###
nworkers <- detectCores() # Initialize the cluster
cl <- makeCluster(nworkers)
registerDoParallel(cl)
##########################
for(a in 1:length(nlist)){
   n = nlist[a]
   filename = paste0("Result_Manuscript/obs_n1D", n, ".RData")
   load(filename)
   output <- foreach (m = 1:M, .packages = c("Matrix", "rSPDE")) %dopar% {
      # m th dataset among M = 50 dataset
      X = df$X[((m-1)*n+1):(m*n)]
      Y = df$Z[((m-1)*n+1):(m*n)]
      result = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = const,
                                kappa.init = kappa, mcmc = target, brn=brn, seed = 1234)
      g.plot = tail(result$g_list, target) # choosing last `target` samples
      obs = data.frame(X, Y)
      y.plot = glist_to_plotdf(g.plot, grid.plot, true = f0_1D, alpha1 = 0.95, alpha2 = 0.9)
      mean((y.plot$true - y.plot$mean)^2)
   }
   MSE_list[, a] = simplify(output)
}
## Plot 3. Fitting of the samples (credible intervals)

################## converting predicted values into usable form ##################
##################################################################################

obs = data.frame(X, Y)

ggplot(y.plot, aes(x = x)) +
   geom_line(aes(y=mean), colour="blue") + 
   geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) + 
   geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
   geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
   geom_point(data = obs, aes(X, Y), size = 0.3) +
   labs(title = paste0("CI 95% of ", target, " samples (n = ", n, ")"))+
   geom_point(data = obs, aes(X, Y), size = 0.3)


N_list = tail(result$N_list, target)
plot(N_list)
lines(N_list)

## Plot 4. For the replicated result ##
library(tidyverse)
MSE.df = data.frame(MSE_list)
colnames(MSE.df) = nlist
MSE.df.plot = gather(MSE.df, key = "n", value = "MSE")
MSE.df.plot$n <- factor(MSE.df.plot$n, levels = c(200, 500, 1000))
ggplot(MSE.df.plot) + 
   geom_boxplot(aes(x = n, y = MSE))
## save as Rdata
# save(MSE.df.plot, file = "Result_Manuscript/MSE_dataframe/MSE_SPDE_1D.RData")

## Time comparison
library(microbenchmark)
m = 1
Xlist = list(length = length(nlist))
Ylist = list(length = length(nlist))
for(a in 1:length(nlist)){
   n = nlist[a]
   filename = paste0("Result_Manuscript/obs_n1D", n, ".RData")
   load(filename)
   Xlist[[a]] = df$X[((m-1)*n+1):(m*n)]
   Ylist[[a]] = df$Z[((m-1)*n+1):(m*n)]
}

# time comparison according to n
target = 500
time_comparison = microbenchmark(
   result1 = sample.exact.seq(Xlist[[1]], Ylist[[1]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                             kappa.init = kappa, mcmc = target, brn=brn, seed = 1234),
   result2 = sample.exact.seq(Xlist[[2]], Ylist[[2]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                             kappa.init = kappa, mcmc = target, brn=brn, seed = 1234),
   result3 = sample.exact.seq(Xlist[[3]], Ylist[[3]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                             kappa.init = kappa, mcmc = target, brn=brn, seed = 1234),
   times = 10
)

time_comparison$n = time_comparison$expr
levels(time_comparison$n) <- nlist
time_comparison$n <- as.numeric(as.character(time_comparison$n))
ggplot(time_comparison, aes(x = n, y = time)) +
   geom_point() +
   geom_smooth() +
   labs(title = "SPDE 1D time comparison")

# time comparison according to s
s_list = c(100, 500, 1000, 5000, 10000)
microbenchmark(
   result1 = sample.exact.seq(Xlist[[1]], Ylist[[1]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                              kappa.init = kappa, mcmc = s_list[1], brn=brn, seed = 1234),
   result2 = sample.exact.seq(Xlist[[1]], Ylist[[1]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                              kappa.init = kappa, mcmc = s_list[2], brn=brn, seed = 1234),
   result3 = sample.exact.seq(Xlist[[1]], Ylist[[1]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                              kappa.init = kappa, mcmc = s_list[3], brn=brn, seed = 1234),
   result4 = sample.exact.seq(Xlist[[1]], Ylist[[1]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                              kappa.init = kappa, mcmc = s_list[4], brn=brn, seed = 1234),
   result5 = sample.exact.seq(Xlist[[1]], Ylist[[1]], sigsq = 0.1^2, Nk = Nk, N.pr = const,
                              kappa.init = kappa, mcmc = s_list[5], brn=brn, seed = 1234),
   times = 10
)



# time comparison according to N
microbenchmark(
   result1 = sample.exact.seq(Xlist[[2]], Ylist[[2]], sigsq = 0.1^2, Nk = c(5, 6, 7), N.pr = const,
                              kappa.init = kappa, mcmc = 1000, brn=brn, seed = 1234),
   result2 = sample.exact.seq(Xlist[[2]], Ylist[[2]], sigsq = 0.1^2, Nk = c(50, 60, 70), N.pr = const,
                              kappa.init = kappa, mcmc = 1000, brn=brn, seed = 1234),
   result3 = sample.exact.seq(Xlist[[2]], Ylist[[2]], sigsq = 0.1^2, Nk = c(200, 240, 280), N.pr = const,
                              kappa.init = kappa, mcmc = 1000, brn=brn, seed = 1234),
   times = 10
)

# R profiling to see detailed time cost
Rprof()
invisible(sample.exact.seq(Xlist[[1]], Ylist[[1]], sigsq = 0.1^2, Nk = c(50, 60, 70), N.pr = const,
                           kappa.init = kappa, mcmc = 100, brn=brn, seed = 1234))
Rprof(NULL)
summaryRprof()
