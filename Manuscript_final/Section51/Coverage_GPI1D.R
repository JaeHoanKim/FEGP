rm(list = ls())
gc()
library(fields)
library(FastGP)
library(Rcpp)
library(ggplot2)
library(Matrix)
library(foreach)
library(doParallel)
library(tidyverse)
source("Result_Manuscript/GraphAesthetics.R")

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

nlist = c(200, 500, 1000)
df_1D = list(length = length(nlist))

target = 2500
brn = 0
brn.ESS = 100
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

##################################################
####### Coverage plot for a specific data ########
##################################################

nlist = c(200, 500, 1000, 10000)

for(i in 1:length(nlist)){
   set.seed(i)
   n = nlist[i]
   # 1D data generation
   X = runif(n)
   Z = f0_1D(X) + rnorm(n) * 0.1
   df_1D[[i]] = data.frame(X, Z)
}

cover.plot.GPI.list = list(length = length(nlist))
cover.plot.GPI.beta2.list = list(length = length(nlist))
## Save plots regarding SPDE

for(a in 1:length(nlist)){
   n = nlist[a] # the number of observed data; 200, 500, 1000
   df = df_1D[[a]]
   X = df$X
   Y = df$Z
   obs = data.frame(X, Y)
}


result.GPI.list = vector(length = length(nlist), "list")

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
   result.GPI.list[[a]] = result.GPI
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
      labs(y = "y") + 
      theme(plot.title = element_text(hjust = 0.5))
   
   cover.plot.GPI.beta2.list[[a]] <- ggplot(y.plot.GPI.beta2, aes(x = x)) +
      geom_point(aes(x = x, y = true), col = 'red', size = 0.5) +
      geom_line(aes(y=mean), colour="blue") +
      geom_ribbon(aes(ymin=low1, ymax=upp1),  alpha=0.4, show.legend=TRUE) +
      geom_ribbon(aes(ymin=lowsup, ymax=uppsup),  alpha=0.2, show.legend=TRUE) +
      labs(y = "y") + 
      theme(plot.title = element_text(hjust = 0.5))
   print(a)
}


library(gridExtra)
pdf(file = "Graphs/coverage_plot_GPI.pdf", width = 12, height = 6)
grid.arrange(cover.plot.GPI.beta2.list[[1]], cover.plot.GPI.list[[1]],
             cover.plot.GPI.beta2.list[[2]], cover.plot.GPI.list[[2]], 
             cover.plot.GPI.beta2.list[[3]], cover.plot.GPI.list[[3]], nrow = 2, as.table = FALSE)
dev.off()


##########################
#### Marginal N, kappa ###
##########################


df_N = data.frame()
df_kappa = data.frame()
for (i in 1:length(nlist)){
   N_list = tail(result.GPI.list[[i]]$N_list, target)
   kappa_list = tail(result.GPI.list[[i]]$kappa_list, target)
   df_now = data.frame(table(factor(N_list, levels = Nk))/target)
   colnames(df_now)[1] = "N_list"
   df_now$n = nlist[i]
   df_N = rbind(df_N, df_now)
   df_now = data.frame(table(factor(kappa_list, levels = kappak))/target)
   colnames(df_now)[1] = "kappa_list"
   df_now$n = nlist[i]
   df_kappa = rbind(df_kappa, df_now)
}
df_N$n = factor(df_N$n)
df_kappa$n = factor(df_kappa$n)
plot_N <- ggplot(df_N) + geom_point(aes(x = N_list, y = Freq, color = n)) + 
   geom_line(aes(x = N_list, y = Freq, group = n, color = n)) + labs(x = "N", y = "frequency") + 
   theme1
plot_kappa <- ggplot(df_kappa) + geom_point(aes(x = kappa_list, y = Freq, color = n)) + 
   geom_line(aes(x = kappa_list, y = Freq, group = n, color = n)) + labs(x = expression(kappa), y = "frequency") + 
   theme1

pdf(file = "Graphs/Marginal_prob_N_kappa_GPI.pdf", width = 10, height = 4)
grid.arrange(plot_N, plot_kappa, nrow = 1)
dev.off()




# plot_N_large <- ggplot(df_N %>% filter(n == 1000 | n == 10000)) + geom_point(aes(x = N_list, y = Freq, color = n)) + 
#    geom_line(aes(x = N_list, y = Freq, group = n, color = n)) + labs(x = "N", y = "frequency") + 
#    theme1
# plot_kappa_large <- ggplot(df_kappa %>% filter(n == 1000 | n == 10000)) + geom_point(aes(x = kappa_list, y = Freq, color = n)) + 
#    geom_line(aes(x = kappa_list, y = Freq, group = n, color = n)) + labs(x = expression(kappa), y = "frequency") + 
#    theme1
# 
# pdf(file = "Graphs/Marginal_prob_N_kappa_GPI_large.pdf", width = 10, height = 4)
# grid.arrange(plot_N_large, plot_kappa_large, nrow = 1)
# dev.off()

