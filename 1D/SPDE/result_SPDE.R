rm(list = ls())
gc()
library(grDevices)
library(ggplot2)
library(Matrix)
source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")
nlist = 300
plabslist = list()
index = 1
n = nlist[index]
X = c(runif(n * 0.8) * 0.2, 0.2 + runif (n * 0.2) * 0.8) 
f0 = function(x){return(1*x^2+sin(8*x))}
Y = f0(X) + rnorm(n) * 0.1
obs = data.frame(X, Y)
## discretized version of 1 over exponential distribution - hich satisfy the condition for prior theoretically
# N.pr = function(N){return (1/N^2 * 3 * exp(-3 / N))}
target = 100
brn = 100
kappa = 3
N.fixed = 20
algo = "exact"
if (algo == "ESS.Nfixed"){
   Temperature = 1
   ## MAKE SURE IF TEMP IS SET TO BE 1!
   result = sample.ESS.Nfixed(X, Y, sigsq = 0.1^2, kappa.init = kappa, N.init = N.fixed,
                              mcmc = target, brn=brn, thin = 1, Temp = Temperature)
} else if (algo == "PT.ESS"){
   ################ For the PTESS algorithm ####################
   N.pr.poi = 5
   Nk = c(3, 5, 10, 15, 20)
   result = sample.PT.ESS(X, Y, sigsq = 0.1^2, kappa.init = 5, Nk = Nk,
                          N.pr = function(x){return(dpois(x, lambda = N.pr.poi))},
                          Tk = c(1, 3, 10, 30, 100),
                          # Tk = c(1:(K - N.min + 1))^4,
                          # N.pr = function(x){return (1)},
                          mcmc = target, brn=brn, thin = 1, pred = FALSE)
} else if(algo == "exact"){
   result = sample.exact(X, Y, sigsq = 0.1^2, kappa.init = kappa, gridsize = N.fixed,
                         mcmc = target, brn=brn, thin = 1)
}




## Plot 3. Fitting of the samples (credible intervals)

#################################################################################
alpha1 = 0.9 # coverage probability 1 (dark gray area)
alpha2 = 0.95 # coverage probability 2 (light gray area)
sample.num = target # number of samples used for the credible intervals
g.plot = tail(result$g_list, sample.num) # choosing last `sample.num` samples
grid.plot = c(0:1000)/1000

y.plot = glist_to_plotdf(g.plot, grid.plot, true = f0, alpha1 = 0.95, alpha2 = 0.9)

ggplot(y.plot, aes(x = x)) +
   geom_line(aes(y=med), colour="blue") + 
   geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) + 
   geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
   geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
   geom_point(data = obs, aes(X, Y), size = 0.3) +
   labs(title = paste0("CI 95% of ", sample.num, " samples using ESS after ", target+brn - sample.num, " burn ins"))+
   geom_point(data = obs, aes(X, Y), size = 0.3)

N_list = tail(result$N_list, target)
plot(N_list)