rm(list = ls())
library(ggplot2)
library(tidyverse)
source("Result_Manuscript/GraphAesthetics.R")
filename = "MSE_list_generated_data"
load(paste0(filename, "_1D.RData"))
MSE_list_1D = data.frame(MSE_list_1D)
M = 50
n_list = c(200, 500, 1000)
MSE_list_1D$method = c(rep("GPI", M), rep("SPDE", M))
colnames(MSE_list_1D) = c(n_list, "method")
MSE_list_1D <- MSE_list_1D %>%
   gather(key = "n", value = "MSE", `200`:`1000`) %>%
   mutate(n = factor(n, levels = c("200", "500", "1000")))
   
ggplot(MSE_list_1D) +
   geom_boxplot(aes(x = factor(n), y = MSE, color = method)) + labs(title = "MSE comparison plot - 1D") +
   theme1
