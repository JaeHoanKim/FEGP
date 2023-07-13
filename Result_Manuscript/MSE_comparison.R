library(ggplot2)
source("Result_Manuscript/GraphAesthetics.R")
fileloc = "Result_Manuscript/MSE_dataframe/"
load(paste0(fileloc, "MSE_GPI_1D.RData"))
MSE.df.GPI <- MSE.df.plot
MSE.df.GPI$method = "GPI"
load(paste0(fileloc, "MSE_SPDE_1D.RData"))
MSE.df.SPDE <- MSE.df.plot
MSE.df.SPDE$method = "SPDE"

library(tidyverse)
MSE.df = bind_rows(MSE.df.GPI, MSE.df.SPDE)
MSE.plot.1D = ggplot(MSE.df) +
   geom_boxplot(aes(x = n, y = MSE, color = method)) + labs(title = "MSE comparison plot - 1D") +
   theme1

load(paste0(fileloc, "MSE_GPI_2D.RData"))
MSE.df.GPI <- MSE.df.plot
MSE.df.GPI$method = "GPI"
load(paste0(fileloc, "MSE_SPDE_2D.RData"))
MSE.df.SPDE <- MSE.df.plot
MSE.df.SPDE$method = "SPDE"
load(paste0(fileloc, "MSE_NNGP_2D.RData"))
MSE.df.NNGP <- MSE.df.plot
MSE.df.NNGP$method = "NNGP"

MSE.df = bind_rows(MSE.df.GPI, MSE.df.SPDE, MSE.df.NNGP)
MSE.plot.2D = ggplot(MSE.df) +
   geom_boxplot(aes(x = n, y = MSE, color = method)) + labs(title = "MSE comparison plot - 2D") +
   theme1


## Save plots ##
library(gridExtra)
# MSE.plot = grid.arrange(MSE.plot.1D, MSE.plot.2D, ncol = 2)
# ggsave("Graphs/MSE_plot.pdf", MSE.plot)


pdf(file = "Graphs/MSE_plot.pdf", width = 12, height = 4)
grid.arrange(MSE.plot.1D, MSE.plot.2D, ncol = 2)
dev.off()
