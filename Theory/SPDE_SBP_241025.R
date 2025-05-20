M = 100
N = 5000
a = replicate(M, {
   w = vector(length = N + 1)
   w[1] = 0.01
   w[2] = 0.010001
   kappa = 1
   for (i in 2:N){
      w[i + 1] = (2 + kappa^2 / N^2) * w[i] - w[i - 1] + rnorm(1) / N / sqrt(N)
   }
   max(abs(w))
})
epsilon = 0.1
mean(a < epsilon)
hist(a)
plot(w)

## integrated Brownian motion
M = 100
N = 1000
vec = vector(length = N+1)
vec[1] = 0
vec[2] = 0
for(i in 3:(N+1)){
   vec[i] = 2 * vec[i-1] - vec[i-2] + rnorm(1) / N / sqrt(N)
}
plot(vec)
