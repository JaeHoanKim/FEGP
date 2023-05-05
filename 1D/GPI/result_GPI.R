rm(list = ls())
gc()
library(grDevices)
library(ggplot2)
require(Rcpp)
sourceCpp("inv_chol.cpp")
source("functions_GPI_sampling.R")
source("functions_GPI.R")

n = 50
# N.max = 100 # for uniform prior - maximum number of knots
plabslist = list()
X = runif(n)
f0 = function(x){return(1*x^2+sin(8*x))}
Y = f0(X) + rnorm(n) * 0.1
obs = data.frame(X, Y)
ntest = 100
Xtest = runif(ntest) # not used
Ytest = sin(5*pi*Xtest) + rnorm(ntest) * 0.1 # not used
## discretized version of 1 over exponential distribution - hich satisfy the condition for prior theoretically
# N.pr = function(N){return (1/N^2 * 3 * exp(-3 / N))}
# N.pr = function(N){return (ifelse(N <= N.max, 1, 0))}
N.fixed = 15
target = 9800
brn = 200
result = sample.ESS.nest(X, Y, N.pr = function(x){return(dpois(x, lambda = 5))}, N.init = 2, 
                     sigsq = 0.1^2, mcmc = target, brn=brn, thin = 1, l.in = 0.5)


################ For the PTESS algorithm ####################
N.pr.poi = 5
Nk = c(3, 5, 10, 15, 20)
result = sample.PT.ESS(X, Y, sigsq = 0.1^2, kappa.init = 5, Nk = Nk,
                       N.pr = function(x){return(dpois(x, lambda = N.pr.poi))},
                       Tk = c(1, 3, 10, 30, 100),
                       # Tk = c(1:(K - N.min + 1))^4,
                       # N.pr = function(x){return (1)},
                       mcmc = target, brn=brn, thin = 1, pred = FALSE)

## Plot 3. Fitting of the samples (credible intervals)

################## converting predicted values into usable form ##################
##################################################################################



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
lines(N_list)
