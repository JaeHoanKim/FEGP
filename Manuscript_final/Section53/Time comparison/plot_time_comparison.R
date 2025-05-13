rm(list = ls())
gc()

library(ggplot2)
source("Result_Manuscript/GraphAesthetics.R")

load("Result_Manuscript/Time_dataframe/time_2D_iter_500.Rdata")

iter.plot <- ggplot(time_comparison_iter_unify) +
   geom_boxplot(aes(x = factor(n), y = log_time, color = method)) +
   labs(title = paste0(500, " samples"), x = "n", y = "log(time)") + theme1

load("Result_Manuscript/Time_dataframe/time_2D_onetime_500.Rdata")

onetime.plot <- ggplot(time_comparison_onetime_unify) +
   geom_boxplot(aes(x = factor(n), y = log_time, color = method)) +
   labs(title = paste0(500, " samples"), x = "n", y = "log(time)") + theme1

# library(gridExtra)
pdf(file = "Graphs/Time_iter_plot.pdf", width = 6, height = 4)
# grid.arrange(time.plot.1, time.plot.2, ncol = 2)
iter.plot
dev.off()
   
