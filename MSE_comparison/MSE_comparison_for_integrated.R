rm(list = ls())
library(ggplot2)
library(tidyverse)
source("Result_Manuscript/GraphAesthetics.R")
filename = "MSE_comparison/comparison_trueftn/MSE_list_generated_data"
MSE_plot_1D = list()
MSE_plot_2D = list()
for(index in c(1:2)){
   ## save plots for 1D
   load(paste0(filename, "_1D_", index, ".RData"))
   MSE_list_1D = data.frame(MSE_list_1D)
   M = 50
   n_list = c(200, 500, 1000)
   MSE_list_1D$method = c(rep("GPI", M), rep("SPDE", M))
   colnames(MSE_list_1D) = c(n_list, "method")
   MSE_list_1D <- MSE_list_1D %>%
      gather(key = "n", value = "MSE", `200`:`1000`) %>%
      mutate(n = factor(n, levels = c("200", "500", "1000")))
   MSE_plot_1D[[index]] <- ggplot(MSE_list_1D) +
      geom_boxplot(aes(x = factor(n), y = MSE, color = method)) + labs(title = "MSE comparison plot - 1D", x = "n") +
      theme1 + ylim(c(0, max(MSE_list_1D$MSE)))
   
   ## save plots for 2D
   load(paste0(filename, "_2D_", index, ".RData"))
   MSE_list_2D = data.frame(MSE_list_2D)
   M = 50
   n_list = c(200, 500, 1000)
   MSE_list_2D$method = c(rep("GPI", M), rep("SPDE", M), rep("NNGP", M))
   colnames(MSE_list_2D) = c(n_list, "method")
   MSE_list_2D <- MSE_list_2D %>%
      gather(key = "n", value = "MSE", `200`:`1000`) %>%
      mutate(n = factor(n, levels = c("200", "500", "1000")))
   MSE_plot_2D[[index]] <- ggplot(MSE_list_2D) +
      geom_boxplot(aes(x = factor(n), y = MSE, color = method)) + labs(title = "MSE comparison plot - 2D", x = "n") +
      theme1
}


## Save plots ##
library(gridExtra)
# MSE.plot = grid.arrange(MSE.plot.1D, MSE.plot.2D, ncol = 2)
# ggsave("Graphs/MSE_plot.pdf", MSE.plot)


# pdf(file = "Graphs/MSE_plot_1D.pdf", width = 12, height = 4)
# grid.arrange(MSE_plot_1D[[1]], MSE_plot_1D[[2]], ncol = 2)
# dev.off()

# pdf(file = "Graphs/MSE_plot_2D.pdf", width = 12, height = 4)
# grid.arrange(MSE_plot_2D[[1]], MSE_plot_2D[[2]], ncol = 2)
# dev.off()

### MSE as the number of burnin differs
library(gridExtra)
filename = "MSE_comparison/comparison_burnin/brn1000/MSE_list_generated_data"

load(paste0(filename, "_2D_1.RData"))
MSE_list_2D = data.frame(MSE_list_2D)

M = 50
n_list = c(200, 500, 1000)
MSE_list_2D$method = c(rep("GPI", M), rep("SPDE", M), rep("NNGP", M))
colnames(MSE_list_2D) = c(n_list, "method")
MSE_list_2D <- MSE_list_2D %>%
   gather(key = "n", value = "MSE", `200`:`1000`) %>%
   mutate(n = factor(n, levels = c("200", "500", "1000")))

MSE_plot_2D_1000 <- ggplot(MSE_list_2D) +
   geom_boxplot(aes(x = factor(n), y = MSE, color = method)) + labs(title = "MSE comparison plot - 2D", x = "n") +
   theme1


pdf(file = "Graphs/MSE_plot_burnin.pdf", width = 12, height = 4)
grid.arrange(MSE_plot_2D[[1]], MSE_plot_2D_1000, ncol = 2)
dev.off()
