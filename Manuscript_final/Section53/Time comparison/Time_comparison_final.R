library(ggplot2)
library(tidyverse)
source("GraphAesthetics.R")

# Load and preprocess
target <- 10
fileloc <- "Manuscript_final/Section53/Time comparison/time_2D_GPI_FullGP_"
load(paste0(fileloc, target, ".RData"))

time_comparison_unify <- time_comparison_unify %>%
   filter(method %in% c("GPI", "FullGP")) %>%
   mutate(
      method = recode(method, "FullGP" = "GP"),
      time = as.numeric(time) / target / 1e9
   )

# Compute log(n) labels
log_n_labels <- time_comparison_unify %>%
   distinct(n) %>%
   arrange(n) %>%
   mutate(log_n = round(log(n), 2)) %>%
   pull(log_n)

# Plot: Boxplot grouped by factor(n), x-axis labeled with log(n)
time.plot <- ggplot(time_comparison_unify, aes(x = factor(n), y = log(time), color = method)) +
   geom_boxplot() +
   labs(
      title = "",
      x = "n",
      y = "log time(s)"
   ) +
   theme1

# Display and save
print(time.plot)
pdf(file = "Graphs/Time_plot_GPI_FullGP.pdf", width = 5, height = 4)
time.plot
dev.off()
