rm(list = ls())
library(ggplot2)
library(tidyverse)
source("GraphAesthetics.R")
kappa_cand = seq(5, 10, length.out = 8)
for (i in 1:length(kappa_cand)){
   kappa_cand[i] = round(kappa_cand[i], 2)
}
MSE_plot_1D = list()
for(i in c(1:length(kappa_cand))){
   filename = paste0("MSE_comparison/MSE_list_SPDE_Matern_1D_kappa", kappa_cand[i], ".RData")
   load(filename)
   MSE_list_1D = data.frame(MSE_list_1D)
   M = 50
   n_list = c(200, 500, 1000)
   MSE_list_1D$method = c(rep("SPDE", M), rep("Matern", M))
   colnames(MSE_list_1D) = c(n_list, "method")
   MSE_list_1D <- MSE_list_1D %>%
      gather(key = "n", value = "MSE", `200`:`1000`) %>%
      mutate(n = factor(n, levels = c("200", "500", "1000")))
   MSE_plot_1D[[i]] <- ggplot(MSE_list_1D) +
      geom_boxplot(aes(x = factor(n), y = log(MSE), color = method)) + 
      labs(
         x = "n", y = "log(AMSE)"
      ) +
      theme1 + ylim(min(log(MSE_list_1D$MSE)), max(log(MSE_list_1D$MSE)))
   
}
## Save plots ##
library(gridExtra)

pdf(file = "Graphs/log_AMSE_plot_SPDE_Matern_kappa.pdf", width = 24, height = 12)
grid.arrange(grobs = MSE_plot_1D, ncol = 4, nrow = 2)
dev.off()
