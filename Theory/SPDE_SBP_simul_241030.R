M = 1000
epsilon = 0.0001
kappa = 0
# N = floor(epsilon^(-2/3)) + 1
N = floor((pi / 2 / epsilon / epsilon)^{ 1/3}) + 1
N = floor(0.01/sqrt(epsilon)) + 1
a = replicate(M, {
   w = runif(N + 2, -epsilon, epsilon)
   - 2 * N^3 * sum((w[3:(N+2)] - (2 + kappa^2 / N^2) * w[2:(N+1)] + w[1:N])^2) + (N+1) * log(2 * epsilon) + (1.5 * N - 1.5) * log(N) - 
       (N+1)/2 * log(2*pi)})
   # (N+1) * log(2 * epsilon) + (1.5 * N - 1.5) * log(N) - 
   #   (N+1)/2 * log(2*pi)})

a
mean(a)

# N = 2 case
N = 10
M = 1000
epsilon_vec = 10^seq(-2, -6, length.out = 30)
Eexp_vec = vector(length = length(epsilon_vec))
for(i in c(1:length(epsilon_vec))){
   epsilon = epsilon_vec[i]
   a = replicate(M, {
      w = runif(3, -epsilon, epsilon)
      w = runif(N + 1, -epsilon, epsilon)
      t = - 2 * N^3 * sum((w[3:(N+1)] - (2 + kappa^2 / N^2) * w[2:N] + w[1:(N-1)])^2) + (N+1) * log(2 * epsilon) + (1.5 * N - 1.5) * log(N) - 
         (N+1)/2 * log(2*pi)
      exp(t)
   })
   Eexp_vec[i] = mean(a)
}
plot(epsilon_vec^3, Eexp_vec)
plot(epsilon_vec, log(Eexp_vec))
plot(log(epsilon_vec), log(Eexp_vec))
max(log(Eexp_vec)) - min(log(Eexp_vec))
max(log(epsilon_vec)) - min(log(epsilon_vec))
