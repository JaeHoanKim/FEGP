rm(list = ls())
library(ggplot2)
library(tidyverse)
source("Result_Manuscript/GraphAesthetics.R")
filename = "MSE_comparison/comparison_trueftn/MSE_list_generated_data"
MSE_plot_1D = list()
M = 50
n_list = c(200, 500, 1000)
## Save plots ##
library(gridExtra)

## Result only for GPI
## save plots for 1D
for(index in c(1:2)){
   load(paste0(filename, "_1D_", index, ".RData"))
   MSE_list_1D = data.frame(MSE_list_1D)
   M = 50
   n_list = c(200, 500, 1000)
   MSE_list_1D$method = c(rep("GPI", M), rep("SPDE", M))
   colnames(MSE_list_1D) = c(n_list, "method")
   MSE_list_1D <- MSE_list_1D %>%
      gather(key = "n", value = "MSE", `200`:`1000`) %>%
      mutate(n = factor(n, levels = c("200", "500", "1000"))) %>%
      filter(method == "GPI")
   MSE_plot_1D[[index]] <- ggplot(MSE_list_1D) +
      geom_boxplot(aes(x = factor(n), y = MSE)) + labs(x = "n", y = "AMSE") +
      theme1 + ylim(c(0, max(MSE_list_1D$MSE))) + 
      scale_y_continuous(labels = function(x) sprintf("%.1e", x))
}
pdf(file = "Graphs/MSE_plot_1D_GPI.pdf", width = 10, height = 4)
grid.arrange(MSE_plot_1D[[1]], MSE_plot_1D[[2]], ncol = 2)
dev.off()