########################################
# 1. generate plots for the moderate n #
########################################

library(ggplot2)
target = 2
fileloc = "Result_Manuscript/Time_dataframe/time_2D_"
load(paste0(fileloc, target, ".RData"))
time_comparison_unify$time = time_comparison_unify$time / 1e9
time.plot.2 <- ggplot(time_comparison_unify) +
   geom_boxplot(aes(x = factor(n), y = time, color = method)) +
   labs(title = paste0(target, " samples"), x = "n", y = "time(s)") + theme1

target = 500
fileloc = "Result_Manuscript/Time_dataframe/time_2D_"
load(paste0(fileloc, target, ".RData"))
time_comparison_unify$time = time_comparison_unify$time / 1e9
time.plot.500 <- ggplot(time_comparison_unify) +
   geom_boxplot(aes(x = factor(n), y = time, color = method)) +
   labs(title = paste0(target, " samples"), x = "n", y = "time(s)") + theme1


library(gridExtra)
pdf(file = "Graphs/Time_whole_plot_updated.pdf", width = 12, height = 4)
grid.arrange(time.plot.2, time.plot.500, ncol = 2)
dev.off()