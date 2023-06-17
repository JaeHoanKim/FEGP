rm(list = ls())
gc()
library(tidyverse)
library(fields)
library(FastGP)
library(ggpubr)
library(Rcpp)
library(ggplot2)
library(Matrix)
sourceCpp("1D/GPI/inv_chol.cpp")
source("1D/GPI/functions_GPI.R")
source("1D/GPI/functions_GPI_sampling.R")

target = 2500
brn.ESS = 100
kappa = 2
dpois5 = function(x){
   return(dpois(x, lambda = 5))
}
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

for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   filename = paste0("Result_Manuscript/obs_n1D", n, ".RData")
   load(filename)
   for(m in 1:M){  
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
      MSE_list[m, a] = mean((y.plot$true - y.plot$mean)^2)
      print(m)
   }
}

## Plot 3. Fitting of the samples (credible intervals)

################## converting predicted values into usable form ##################
##################################################################################



ggplot(y.plot, aes(x = x)) +
   geom_line(aes(y=mean), colour="blue") + 
   geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) + 
   geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
   geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
   # geom_point(data = obs, aes(X, Y), size = 0.3) +
   labs(title = paste0("CI 95% of ", target, " samples (n = ", n, ")"))


N_list = tail(result$N_list, target)
plot(N_list)
lines(N_list)

## Plot 4. For the replicated result ##
MSE.df = data.frame(MSE_list)
colnames(MSE.df) = nlist
MSE.df.plot = gather(MSE.df, key = "n", value = "MSE")
MSE.df.plot$n <- factor(MSE.df.plot$n, levels = c(200, 500, 1000))
   

## save as Rdata
# save(MSE.df.plot, file = "Result_Manuscript/MSE_dataframe/MSE_GPI_1D.RData")

ggplot(MSE.df.plot) + 
   geom_boxplot(aes(x = n, y = MSE)) + themegg


# time comparison according to n
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

target = 200
# time comparison according to n
microbenchmark(
   result1 = sample.ESS.seq(Xlist[[1]], Ylist[[1]], sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = target, brn=0, brn.ESS = brn.ESS),
   result2 = sample.ESS.seq(Xlist[[2]], Ylist[[2]], sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = target, brn=0, brn.ESS = brn.ESS),
   result3 = sample.ESS.seq(Xlist[[3]], Ylist[[3]], sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = target, brn=0, brn.ESS = brn.ESS),
   times = 10
)

# time comparison according to s
microbenchmark(
   result1 = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = 100, brn=0, brn.ESS = brn.ESS),
   result2 = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = 500, brn=0, brn.ESS = brn.ESS),
   result3 = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = Nk, nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = 1000, brn=0, brn.ESS = brn.ESS),
   times = 10
)

# time comparison according to N
microbenchmark(
   result1 = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = c(4, 6), nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = target, brn=0, brn.ESS = brn.ESS),
   result2 = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = c(8, 12), nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = target, brn=0, brn.ESS = brn.ESS),
   result3 = sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = c(16, 24), nu.in = 1, l.in = 1/kappa,
                            N.pr = const, mcmc = target, brn=0, brn.ESS = brn.ESS),
   times = 10
)

# R profiling to see detailed time cost
Rprof()
invisible(sample.ESS.seq(X, Y, sigsq = 0.1^2, Nk = c(4, 6), nu.in = 1, l.in = 1/kappa,
                         N.pr = const, mcmc = target, brn=0, brn.ESS = brn.ESS))
Rprof(NULL)
summaryRprof()