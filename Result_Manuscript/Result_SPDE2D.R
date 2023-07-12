rm(list = ls())
gc()
library(Matrix)
library(ggplot2)
library(rSPDE)
library(ggpubr)
source("2D/SPDE/functions_SPDE_sampling_2D.R")
source("2D/SPDE/functions_SPDE_2D.R")
source("Result_Manuscript/GraphAesthetics.R")

kappa = 2
brnin = 0
target = 2500

const = function(x){
   return(1)
}
dpoi5 = function(x){
   return(dpois(x, lambda = 5))
}
f0 = function(x, y){
   return(sin(5*x + 2*y) + 2*y^2)
}
Nk = c(3, 5, 8, 10, 15)
gridsize = 40
M = 50
nlist = c(200, 500, 1000)
# gridmat is a (gridsize^2) by 2 matrix!
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

############################################################
a = 2
n = nlist[a]
filename = paste0("Result_Manuscript/obs_n2D", n, ".RData")
load(filename)
m = 1
X = df[((m-1)*n+1):(m*n), c(1, 2)]
Z = df$Z[((m-1)*n+1):(m*n)]
result = sample.exact2D.seq(X, Z, sigsq = 0.1^2, # N.pr = function(x){return(1)},
                            N.pr = dpoi5,
                            Nk = Nk, kappa.init = kappa, mcmc = target, brn = brnin, seed = 1234)
g_list = result$g_list
y.plot = glist_to_plotdf_2D(g_list, gridmat, truefun = f0, alpha1 = 0.9, alpha2 = 0.95)

################## plot ###################
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


final_plot = ggarrange(plotlist = list(plot_true, plot_mean, 
                                       plot_low2, plot_upp2), nrow = 2, ncol = 2)


fileloc = "Result_Manuscript/heatmap/"
pdf(file = paste0(fileloc, "heatmap_SPDE_2D.pdf"))
print(final_plot)
dev.off()


N_list = tail(result$N_list, target)
par(mfrow = c(1, 2))
plot(N_list, xlab = "Index", ylab = "N")
lines(N_list)
N_list <- factor(N_list)
barplot(table(N_list))
