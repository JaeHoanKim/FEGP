rm(list = ls())
gc()
library(tidyverse)
library(ggpubr)
library(Rcpp)
library(ggplot2)
library(Matrix)
library(fields)
library(FastGP)
source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")
source("GraphAesthetics.R")

target = 2500
nlist = 500
M = 1
f0_2D = function(x, y){return(sin(5*x + 2*y) + 2*y^2)}

N_supp = c(6, 10, 14, 18)
N.pr = function(x){return(rep(1, length(x)))}
kappa.pr = function(x){return(dgamma(x, 3, 1/3))}
kappa.sampler = function(){rgamma(1, 3, 1/3)}
tausq.pr = function(x){return(dgamma(x, 1, 1))}
tausq.sampler = function(){rgamma(1, 1, 1)}
beta = 2
sigma = 0.1
gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

############################################################
## specify n
for(i in 1:length(nlist)){
   set.seed(i)
   n = nlist[i]
   # 2D data generation
   X = matrix(runif(2*n*M), n*M)
   Z = f0_2D(X[, 1], X[, 2]) + rnorm(n*M) * 0.1
}

result = sample.GPI2D(Z = Z, X = X, Nk = N_supp, N.pr = N.pr, 
                      kappa.pr = kappa.pr, kappa.sampler = kappa.sampler,
                      tausq.pr = tausq.pr, tausq.sampler = tausq.sampler, 
                      sigsq = sigma^2, beta = beta, mcmc = target, brn = 1000, brn.ESS = brn.ESS)
g_list = result$g_list
y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0_2D, alpha1 = 0.9, alpha2 = 0.95)

################## plot ###################
library(ggpubr)

grandmin <- round(min(y.plot$truefun) - 2, 2)
grandmax <- round(max(y.plot$truefun) + 2, 2)
mybreaks <- seq(grandmin, grandmax, length.out = 11)
mycolors<- function(x) {
   colors<-colorRampPalette(c("yellow", "blue"))( 10 )
   colors[1:x]
}
#Function to create labels for legend
breaklabel <- function(x){
   labels<- paste0(mybreaks[1:10], "-", mybreaks[2:11])
   labels[1:x]
}

plot_true <- ggplot(y.plot, aes(x1, x2)) +
   geom_contour_filled(aes(z = truefun), breaks = mybreaks,
                       show.legend = TRUE) + labs(x = expression(x[1]), y = expression(x[2])) +  themegg

plot_mean <- ggplot(y.plot, aes(x1, x2)) +
   geom_contour_filled(aes(z = mean), breaks = mybreaks, 
                       show.legend = TRUE) + labs(x = expression(x[1]), y = expression(x[2])) +  themegg

plot_low2 <- ggplot(y.plot, aes(x1, x2)) +
   geom_contour_filled(aes(z = low2), breaks = mybreaks, 
                       show.legend = TRUE) + labs(x = expression(x[1]), y = expression(x[2])) +  themegg

plot_upp2 <- ggplot(y.plot, aes(x1, x2)) +
   geom_contour_filled(aes(z = upp2), breaks = mybreaks, 
                       show.legend = TRUE) + labs(x = expression(x[1]), y = expression(x[2])) +  themegg

## plot arrangement : list(1, 2, 3, 4) => 1 2 // 3 4
final_plot = ggarrange(plotlist = list(plot_true, plot_mean, 
                                       plot_low2, plot_upp2), nrow = 2, ncol = 2)

fileloc = "Graphs/"
pdf(file = paste0(fileloc, "heatmap_GPI_2D_", n, ".pdf"))
print(final_plot)
dev.off()

N_list = tail(result$N_list, target)
par(mfrow = c(1, 2))
plot(N_list, xlab = "Index", ylab = "N")
lines(N_list)
N_list <- factor(N_list)
barplot(table(N_list))