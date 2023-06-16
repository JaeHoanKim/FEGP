fileloc = "Result_Manuscript/MSE_dataframe/"
load(paste0(fileloc, "MSE_GPI_1D.RData"))
MSE.df.GPI <- MSE.df.plot
MSE.df.GPI$method = "GPI"
load(paste0(fileloc, "MSE_SPDE_1D.RData"))
MSE.df.SPDE <- MSE.df.plot
MSE.df.SPDE$method = "SPDE"

library(tidyverse)
MSE.df = bind_rows(MSE.df.GPI, MSE.df.SPDE)
ggplot(MSE.df) +
   geom_boxplot(aes(x = n, y = MSE, color = method))
