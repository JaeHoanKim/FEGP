rm(list = ls())
gc()
library(Matrix)
library(grDevices)
library(ggplot2)
source("Sampling_functions_ESS.R")
source("Sampling_GMRF_2D.R")

n = 300 # the number of observed data
f0 = function(x, y){
   return(sin(11*x + 2*y) + 2*y^2)
} # true function
X = matrix(runif(2*n), n)
Z = f0(X[, 1], X[, 2]) + rnorm(n) * 0.1
## discretized version of 1 over exponential distribution - which satisfy the condition for prior theoretically
# N.pr = function(N){return (1/N^2 * 3 * exp(-3 / N))}

kappa = 2
N.init = 10
brnin = 1000
target = 500
result = sample.ESS.Nfixed2D(X, Z, sigsq = 0.1^2, kappa.init = kappa, N.init = N.init,
                             mcmc = target, brn=brnin, thin = 1)

gridsize = 40
gridmat = cbind(rep(c(0:gridsize)/gridsize, each = gridsize + 1),
                rep(c(0:gridsize)/gridsize, gridsize+ 1))

glist = tail(result$g_list, target)
ghat = Reduce("+", glist) / length(glist) 

Ztrue = f0(gridmat[, 1], gridmat[, 2])
Zhat = f_N_h_2D_multi(gridmat, ghat)
data_plot = data.frame(cbind(gridmat, Zhat, Ztrue))
colnames(data_plot) = c('x', 'y', 'Zhat', 'Ztrue')

grandmin <- round(min(Ztrue), 2)
grandmax <- round(max(Ztrue), 2)
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


MSE = mean((Zhat - Ztrue)^2)
MSE

##################### plot ###################

themegg = theme(
   # LABLES APPEARANCE
   panel.grid.major = element_blank(), 
   panel.grid.minor = element_blank(),
   panel.background = element_rect(fill = "transparent",colour = NA),
   plot.background = element_rect(fill = "transparent",colour = NA),
   plot.title = element_text(hjust = 0.5, size=14, face= "bold", colour= "black" ),
   axis.title.x = element_text(size=20, face="bold", colour = "black"),    
   axis.title.y = element_text(size=20, face="bold", colour = "black"),    
   axis.text.x = element_text(size=18, colour = "black"), 
   axis.text.y = element_text(size=18, colour = "black"),
   strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
   strip.text.y = element_text(size = 12, face="bold", colour = "black"),
   strip.background =  element_rect(fill = "transparent",colour = NA),
   axis.line.x = element_line(color="black", linewidth = 0.2),
   axis.line.y = element_line(color="black", linewidth =  0.2),
   panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.2),
   legend.title=element_blank(),
   legend.text=element_text(size=12, face="bold", colour = "black"),
   legend.position="none"
)

plot_true <- ggplot(data_plot, aes(x, y)) +
   geom_contour_filled(aes(z = Ztrue), breaks = mybreaks, show.legend = TRUE) +
   themegg

plot_mean <- ggplot(data_plot, aes(x, y)) +
   geom_contour_filled(aes(z = Zhat), breaks = mybreaks, show.legend = FALSE) +
   themegg

library(ggpubr)

plot_true
plot_mean
final_plot = ggarrange(plotlist = list(plot_mean, plot_true), nrow = 1)
final_plot

pdf(file = "poster_SPDE_comparison2DD.pdf", width = 16, height = 8)
final_plot
dev.off()                      

# ggplot(data_plot, aes(x, y)) +
#    geom_contour_filled(aes(z = abs(Zhat - Ztrue)), breaks = mybreaks, show.legend = TRUE)+
#    labs(title = paste0("difference plot between true function and estimated mean, MSE = ", MSE))