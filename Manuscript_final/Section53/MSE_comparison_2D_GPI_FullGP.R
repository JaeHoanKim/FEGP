rm(list = ls())
library(ggplot2)
library(tidyverse)
source("GraphAesthetics.R")

fileloc <- "MSE_comparison/"
nlist <- c(200, 500, 1000, 2000)

# Load GPI
load(paste0(fileloc, "MSE_list_generated_data_2D_GPI.RData"))
MSE.df.GPI <- as.data.frame(MSE_list_GPI2D)
colnames(MSE.df.GPI) <- nlist
MSE.df.GPI$method <- "GPI"

# Load FullGP
load(paste0(fileloc, "MSE_list_generated_data_2D_FullGP.RData"))
MSE.df.FullGP <- as.data.frame(MSE_list_FullGP2D)
colnames(MSE.df.FullGP) <- nlist
MSE.df.FullGP$method <- "Matern"

# Combine and reshape
MSE.df <- bind_rows(MSE.df.GPI, MSE.df.FullGP) %>%
   pivot_longer(cols = as.character(nlist), names_to = "n", values_to = "MSE") %>%
   mutate(
      n = factor(n, levels = as.character(nlist)),
   )

# Plot with log(n) x-axis labels
MSE.plot.2D <- ggplot(MSE.df, aes(x = n, y = log(MSE), color = method)) +
   geom_boxplot() +
   labs(
      x = "n",
      y = "log(AMSE)"
   ) +
   themegg

# Save
pdf(file = "Graphs/MSE_plot_2D_GPI_FullGP.pdf", width = 6, height = 4)
grid.arrange(MSE.plot.2D, ncol = 1)
dev.off()
