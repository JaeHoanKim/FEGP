rm(list = ls())
gc()
library(grDevices)
library(ggplot2)
library(Matrix)
source("1D/SPDE/functions_SPDE_sampling.R")
source("1D/SPDE/functions_SPDE.R")

n = 200 # the number of observed data; 200, 500, 1000
filename = paste0("Result_Manuscript/obs_n1D", n, ".RData")
load(filename)
m = 1 # m th dataset among M = 50 dataset
X = df$X[((m-1)*n+1):(m*n)]
Y = df$Z[((m-1)*n+1):(m*n)]
target = 250
brn.ESS = 100
kappa = 2
dpois5 = function(x){
   return(dpois(x, lambda = 5))
}

## discretized version of 1 over exponential distribution - hich satisfy the condition for prior theoretically
# N.pr = function(N){return (1/N^2 * 3 * exp(-3 / N))}
target = 100
brn = 100
kappa = 2
Nk = c(3, 5, 10, 15, 20)
result = sample.exact.seq(X, Y, sigsq = 0.1^2, Nk = Nk, N.pr = dpois5,
                          kappa.init = kappa, mcmc = target, brn=brn)

## Plot 3. Fitting of the samples (credible intervals)

################## converting predicted values into usable form ##################
##################################################################################

f0_1D = function(x){return(1*x^2+sin(8*x))}
alpha1 = 0.9 # coverage probability 1 (dark gray area)
alpha2 = 0.95 # coverage probability 2 (light gray area)
sample.num = target # number of samples used for the credible intervals
g.plot = tail(result$g_list, sample.num) # choosing last `sample.num` samples
grid.plot = c(0:1000)/1000
obs = data.frame(X, Y)

y.plot = glist_to_plotdf(g.plot, grid.plot, true = f0_1D, alpha1 = 0.95, alpha2 = 0.9)

ggplot(y.plot, aes(x = x)) +
   geom_line(aes(y=med), colour="blue") + 
   geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) + 
   geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
   geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
   geom_point(data = obs, aes(X, Y), size = 0.3) +
   labs(title = paste0("CI 95% of ", target, " samples (n = ", n, ")"))+
   geom_point(data = obs, aes(X, Y), size = 0.3)


N_list = tail(result$N_list, target)
plot(N_list)
lines(N_list)
