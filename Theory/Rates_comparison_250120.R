library(ggplot2)
beta = 4
alphas = seq(1, 10, 0.1)
lower_bound = pmin(alphas, beta - 1/2) / (2 * pmin(alphas, beta - 1/2) + 1)
# From Information Rates (2009) paper
upper_bound = pmin(alphas, beta - 1/2) / (2 * beta)
optimal = alphas / (2*alphas + 1)
ggplot(df = data.frame(alphas, lower_bound, upper_bound, optimal)) + 
   geom_line(aes(alphas, lower_bound), col = 'red') + 
   geom_line(aes(alphas, upper_bound), col = 'blue')
   