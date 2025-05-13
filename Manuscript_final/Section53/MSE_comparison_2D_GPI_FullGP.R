library(ggplot2)
library(tidyverse)
source("Result_Manuscript/GraphAesthetics.R")
fileloc = "MSE_comparison/"
nlist = c(200, 500, 1000, 2000)
load(paste0(fileloc, "MSE_list_generated_data_2D_GPI.RData"))
MSE.df.GPI <- as.data.frame(MSE_list_GPI2D)
colnames(MSE.df.GPI) <- nlist
MSE.df.GPI$method = "GPI"
load(paste0(fileloc, "MSE_list_generated_data_2D_FullGP.RData"))
MSE.df.FullGP <- as.data.frame(MSE_list_FullGP2D)
colnames(MSE.df.FullGP) <- nlist
MSE.df.FullGP$method = "Matern"
# load(paste0(fileloc, "MSE_SPDE_2D.RData"))
# MSE.df.SPDE <- MSE.df.plot
# MSE.df.SPDE$method = "SPDE"
# load(paste0(fileloc, "MSE_NNGP_2D.RData"))
# MSE.df.NNGP <- MSE.df.plot
# MSE.df.NNGP$method = "NNGP"

MSE.df = bind_rows(MSE.df.GPI, MSE.df.FullGP)
MSE.df = pivot_longer(MSE.df, cols = as.character(nlist), names_to = "n", values_to = "MSE")
MSE.df$n = factor(MSE.df$n, levels = as.character(nlist))
MSE.plot.2D = ggplot(MSE.df) +
   geom_boxplot(aes(x = n, y = MSE, color = method)) + labs(y = "AMSE") +
   theme1
## Save plots ##
library(gridExtra)
# MSE.plot = grid.arrange(MSE.plot.1D, MSE.plot.2D, ncol = 2)
# ggsave("Graphs/MSE_plot.pdf", MSE.plot)


pdf(file = "Graphs/MSE_plot_2D_GPI_FullGP.pdf", width = 6, height = 4)
grid.arrange(MSE.plot.2D, ncol = 1)
dev.off()
