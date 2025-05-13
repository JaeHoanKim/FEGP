########################################
# 1. generate plots for the moderate n #
########################################

library(ggplot2)
source("GraphAesthetics.R")
target = 10
fileloc = "Manuscript_final/Section53/Time comparison/Time_dataframe/time_2D_GPI_FullGP_"
load(paste0(fileloc, target, ".RData"))
time_comparison_unify = time_comparison_unify %>%
   filter(method == "GPI" | method == "FullGP")
# Change column name
time_comparison_unify$method <- as.character(time_comparison_unify$method)  # Convert factor to character
time_comparison_unify$method[time_comparison_unify$method == "FullGP"] <- "Matern"  # Replace values
time_comparison_unify$method <- as.factor(time_comparison_unify$method)  # Convert back to factor

time_comparison_unify$time = as.numeric(time_comparison_unify$time) / target / 1e9
time.plot <- ggplot(time_comparison_unify) +
   geom_boxplot(aes(x = factor(n), y = time, color = method)) +
   labs(title = paste0(""), x = "n", y = "time(s)") + theme1 +
   ylim(0, NA)  
time.plot
pdf(file = "Graphs/Time_plot_GPI_FullGP.pdf", width = 5, height = 4)
time.plot
dev.off()

target = 500
load(paste0(fileloc, target, ".RData"))
time_comparison_unify = time_comparison_unify %>%
   filter(method == "GPI" | method == "FullGP")
time_comparison_unify$time = as.numeric(time_comparison_unify$time) / 1e9
time.plot.500 <- ggplot(time_comparison_unify) +
   geom_boxplot(aes(x = factor(n), y = time, color = method)) +
   labs(title = paste0(target, " samples"), x = "n", y = "time(s)") + theme1 +
   ylim(0, NA)  


library(gridExtra)
pdf(file = "Graphs/Time_whole_plot_mod_n_GPI_FullGP.pdf", width = 10, height = 4)
grid.arrange(time.plot.2, time.plot.500, ncol = 2)
dev.off()

#####################################
# 2. generate plots for the large n #
#####################################

target = 2
fileloc = "Result_Manuscript/Time_dataframe/time_2D_large_"
load(paste0(fileloc, target, ".RData"))
time_comparison_large_n = time_comparison_unify %>%
   filter(method == "GPI" | method == "NNGP")
time_comparison_large_n$time = time_comparison_large_n$time / 1e9 
time.plot.2 <- ggplot(time_comparison_large_n) +
   geom_boxplot(aes(x = factor(n), y = time, color = method)) +
   labs(title = paste0(target, " samples"), x = "n", y = "time(s)") + theme1

target = 500
fileloc = "Result_Manuscript/Time_dataframe/time_2D_large_"
load(paste0(fileloc, target, ".RData"))
time_comparison_large_n = time_comparison_unify %>%
   filter(method == "GPI" | method == "NNGP")
time_comparison_large_n$time = time_comparison_large_n$time / 1e9 
time.plot.500 <- ggplot(time_comparison_large_n) +
   geom_boxplot(aes(x = factor(n), y = time, color = method)) +
   labs(title = paste0(target, " samples"), x = "n", y = "time(s)") + theme1


library(gridExtra)
pdf(file = "Graphs/Time_whole_plot_large_n_GPI_NNGP.pdf", width = 12, height = 4)
grid.arrange(time.plot.2, time.plot.500, ncol = 2)
dev.off()