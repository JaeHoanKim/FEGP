rm(list = ls())
gc()
library(plotly)
library(tidyverse)
library(ggpubr)
library(Rcpp)
library(ggplot2)
source("2D/GPI/functions_GPI_2D.R")
source("2D/GPI/functions_GPI_sampling_2D.R")

kappa = 2
brn.ESS = 100
target = 250
nlist = c(200, 500, 1000)
M = 500
const = function(x){
   return(1)
}
dpoi5 = function(x){
   return(dpois(x, lambda = 5))
}
Nk = c(4, 6, 8, 10)
f0 = function(x, y){
   return(sin(5*x + 2*y) + 2*y^2)
}
gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

############################################################
## specify n
a =  2
n = nlist[a]
filename = paste0("Result_Manuscript/obs_n2D", n, ".RData")
load(filename)
## specify dataset
m = 1
X = df[((m-1)*n+1):(m*n), c(1, 2)]
Z = df$Z[((m-1)*n+1):(m*n)]
result = sample.RJESS2D.seq(Z = Z, X = X, N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                            mcmc = target, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS)
g_list = result$g_list
################## plot ###################
library(ggpubr)

# gridmat is a (gridsize^2) by 2 matrix!
y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0, alpha1 = 0.9, alpha2 = 0.95)

grandmin <- round(min(y.plot$truefun) - 0.5, 2)
grandmax <- round(max(y.plot$truefun) + 0.5, 2)
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

## plot arrangement : list(1, 2, 3, 4) => 1 2 // 3 4
final_plot = ggarrange(plotlist = list(plot_true, plot_mean, 
                                       plot_low2, plot_upp2), nrow = 2, ncol = 2)

fileloc = "Result_Manuscript/heatmap/"
pdf(file = paste0(fileloc, "heatmap_GPI_2D.pdf"))
print(final_plot)
dev.off()

N_list = tail(result$N_list, target)
par(mfrow = c(1, 2))
plot(N_list, xlab = "Index", ylab = "N")
lines(N_list)
N_list <- factor(N_list)
barplot(table(N_list))


####  time comparison according to n
library(microbenchmark)
m = 1
nlist = c(200, 500, 1000)
Xlist = list(length = length(nlist))
Zlist = list(length = length(nlist))
for(a in 1:length(nlist)){
   n = nlist[a]
   filename = paste0("Result_Manuscript/obs_n2D", n, ".RData")
   load(filename)
   Xlist[[a]] = df[((m-1)*n+1):(m*n), c(1, 2)]
   Zlist[[a]] = df$Z[((m-1)*n+1):(m*n)]
}

Nk = c(4, 6, 8, 10, 12)
brn.ESS = 100

# time comparison according to s
microbenchmark(
   result1 = sample.RJESS2D.seq(Z = Zlist[[1]], X = Xlist[[1]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                mcmc = 100, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS),
   result2 = sample.RJESS2D.seq(Z = Zlist[[1]], X = Xlist[[1]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                mcmc = 500, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS),
   result3 = sample.RJESS2D.seq(Z = Zlist[[1]], X = Xlist[[1]], N.pr = function(x){return(1)}, Nk = Nk, sigsq = 0.1^2,
                                mcmc = 1000, brn = 0, nu.in = 1, l.in = 1/kappa, brn.ESS = brn.ESS),
   times = 10
)


