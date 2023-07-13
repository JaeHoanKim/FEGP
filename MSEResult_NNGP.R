rm(list = ls())
gc()
# install.packages("spNNGP")
library(spNNGP)
source("Result_Manuscript/GraphAesthetics.R")
############################# our example ####################################

kappa = 2
brnin = 0
target = 2500
sigsq = 0.01

Nk = c(3, 5, 8, 10, 15)
gridsize = 40
M = 50
nlist = c(200, 500, 1000)
# gridmat is a (gridsize^2) by 2 matrix!
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

########################################################
a = 2
n = nlist[a]
filename = paste0("Result_Manuscript/obs_n2D", n, ".RData")
load(filename)
m = 1
X = as.matrix(df[((m-1)*n+1):(m*n), c(1, 2)])
Z = as.matrix(df$Z[((m-1)*n+1):(m*n)])



starting <- list("phi" = 1/kappa, "sigma.sq" = 1, "tau.sq"=0.01, "nu" = 1)
tuning <- list("phi"= 0, "sigma.sq"= 0, "tau.sq"= 0, "nu" = 0)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1), "nu.unif" = c(1, 1.001))
cov.model <- "matern"
n.report <- 10 #any small integer


## Response
m.r <- spNNGP(Z ~ X-1, coords=X, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=target, n.omp.threads=1, n.report=n.report)
p.r <- predict(m.r, X.0 = gridmat, coords.0 = gridmat, n.omp.threads=1)

pred.grid <- p.r$p.y.0


## data for plot
f0 = function(x, y){
   return(sin(5*x + 2*y) + 2*y^2)
}
mean.grid <- apply(pred.grid, 1, mean)
upp.grid <- apply(pred.grid, 1, FUN = function(x){return(quantile(x, 0.975))})
low.grid <- apply(pred.grid, 1, FUN = function(x){return(quantile(x, 0.025))})
true.grid <- f0(gridmat[, 1], gridmat[, 2])

y.plot = as.data.frame(cbind(mean = mean.grid, truefun = true.grid,
                         upp = upp.grid, low = low.grid, 
                         x1 = gridmat[, 1], x2 = gridmat[, 2]))

## plots

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
   geom_contour_filled(aes(z = low), breaks = mybreaks, 
                       show.legend = TRUE) + themegg

plot_upp2 <- ggplot(y.plot, aes(x1, x2)) +
   geom_contour_filled(aes(z = upp), breaks = mybreaks, 
                       show.legend = TRUE) + themegg


final_plot = ggarrange(plotlist = list(plot_true, plot_mean, 
                                       plot_low2, plot_upp2), nrow = 2, ncol = 2)
final_plot
