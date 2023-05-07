rm(list = ls())
gc()
library(Matrix)
library(grDevices)
library(ggplot2)
source("2D/SPDE/functions_SPDE_sampling_2D.R")
source("2D/SPDE/functions_SPDE_2D.R")

n = 300 # the number of observed data
f0 = function(x, y){
   return(sin(11*x + 2*y) + 2*y^2)
} # true function
X = matrix(runif(2*n), n)
Z = f0(X[, 1], X[, 2]) + rnorm(n) * 0.1
## discretized version of 1 over exponential distribution - which satisfy the condition for prior theoretically
# N.pr = function(N){return (1/N^2 * 3 * exp(-3 / N))}

kappa = 4
N.init = 10
brnin = 1000
target = 500
algo = "ESS"
if (algo == "ESS.Nfixed"){
   result = sample.ESS.Nfixed2D(X, Z, sigsq = 0.1^2, kappa.init = kappa, N.init = N.init,
                          mcmc = target, brn=brnin, thin = 1)
}else if (algo == "ESS"){
   result = sample.ESS.2D(X, Z, sigsq = 0.1^2, kappa.init = kappa, N.init = N.init,
                          mcmc = target, brn=brnin, thin = 1)
}


g_list = tail(result$g_list, target)

##################### plot ###################
library(ggpubr)

gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

# gridmat is a (gridsize^2) by 2 matrix!
y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0, alpha1 = 0.9, alpha2 = 0.95)

grandmin <- round(min(y.plot$truefun) - 0.3, 2)
grandmax <- round(max(y.plot$truefun) + 0.3, 2)
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
                       show.legend = TRUE) +  themegg
  
plot_mean <- ggplot(y.plot, aes(x1, x2)) +
   geom_contour_filled(aes(z = mean), breaks = mybreaks, 
                       show.legend = TRUE) + themegg

plot_low2 <- ggplot(y.plot, aes(x1, x2)) +
   geom_contour_filled(aes(z = low2), breaks = mybreaks, 
                       show.legend = TRUE) + themegg

plot_upp2 <- ggplot(y.plot, aes(x1, x2)) +
   geom_contour_filled(aes(z = upp2), breaks = mybreaks, 
                       show.legend = TRUE) + themegg


final_plot = ggarrange(plotlist = list(plot_true, plot_mean, 
                                       plot_low2, plot_upp2), nrow = 2, ncol = 2)
final_plot

plot(result$N_list)
lines(result$N_list)
MSE = mean((y.plot$truefun - y.plot$mean)^2)
MSE