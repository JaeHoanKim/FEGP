rm(list = ls())
gc()
beta = 2
i = 6
a = vector(length = beta + 1)
for (j in 1:(beta+1)){
   a[j] = (-1)^j * choose(beta, j) * j^i
}
print(a)
sum(a)
library(mvtnorm)
N = 30
kappa = 1
epsilon = 0.01
G = matrix(0, N+1, N+1)
for (i in 2:N){
   G[i, i] = 2 * N
   G[i, i-1] = -N
   G[i, i+1] = -N
}
C = diag(N+1) / N
C[1, 1] = 1/(2 * N)
C[N + 1, N + 1] = 1/(2 * N)
B = kappa^2 * C +G
Q2 = t(B) %*% solve(C) %*% B
Q4 = t(B) %*% solve(C) %*% Q2 %*% solve(C) %*% B
Q4inv = solve(Q4[2:N, 2:N])

A = kappa^2 * diag(N+1) + N * G
Ainv = solve(A)
Varbeta = C
for(i in 1:(beta/2)){
   Varbeta = Ainv %*% Varbeta %*% t(Ainv)
}

y = Varbeta[N/2 + 1, ]
x = Varbeta[c(1:(beta /2 + 1), (N-beta/2 + 1):(N + 1)), ]
df = data.frame(y, t(x))
model = lm(y~. -1, data =df)
predicted_values <- predict(model)
# Calculate the sum of squared errors (SSE)
sse <- sum((y - predicted_values)^2)
print(sse)


a = solve(B, c(0, rep(1, N-1), 0))
max(a)
B_sub = B[2:N, 2:N]
Bsub_inv = solve(B_sub)
A = Bsub_inv %*% Bsub_inv
eg = eigen(G)
max(A)
max(solve(B)[2:N, 2:N])
VarU = t(G) %*% solve(B) %*% solve(C) %*% solve(t(B)) %*% G
VarUinv = solve(VarU[2:N, 2:N])
plot(diag(VarUinv))

K = solve(B) %*% solve(C) %*% solve(t(B))
K
Q = t(B) %*% solve(C) %*% B
Q2 = t(B) %*% solve(C) %*% Q %*% solve(C) %*% B
nsamp = 10
Sigma = solve(Q2)
Sigma = Sigma + t(Sigma)
x = rmvnorm(nsamp, sigma = Sigma)

eg = eigen(Q)
1 / sqrt(eg$values)
nsamp = 10
x = rmvnorm(nsamp, sigma = solve(Q))
plot(x[1,])
for(i in 1:nsamp){
   lines(x[i, ])
}
Qinv = solve(Q)

A = matrix(c(1, 0, 0, 1), 2, 2)
Sigma2 = A %*% Qinv[1:2, 1:2] %*% t(A)
min(eigen(Sigma2)$values)


Btinv = solve(t(B))
Vw = t(Btinv) %*% C %*% Btinv
Q = solve(Vw)
sum(eigen(Q)$values)
tmp = -3
for(i in 1:N){
   j = i+1
   print(Vw[i, i] + Vw[j, j] - Vw[i, j] - Vw[j, i])
}
a = rmvnorm(1, sigma = Vw)[,]
plot(a, type = "l")
A = matrix(0, N-1, N+1)
for (i in 1:(N-1)){
   A[i, i] = 1
   A[i, i+1] = -2
   A[i, i+2] = 1
}
Fmat = diag(N+1) -kappa^2 * Btinv %*% C
V = t(Fmat) %*% C %*% Fmat / N^2
V = V[2:N, 2:N]
a = rmvnorm(1, sigma = V)
a = a[,]
plot(a, type = "l")
P = diag(N-1)
for (i in 1:(N-1)){
   for(j in 1:i){
      P[i, j] = 1
   }
}
b = P %*% a
plot(cumsum(b))
max(diag(V))
max(V)
min(V)
heatmap(V, Rowv= NA, Colv = NA)
BTB = t(B) %*% B
eigen(BTB)$values

Binv1 = solve(B)
## explicit calculation - theta[i] denotes theta_{i-1} in notation
theta = vector(length = N+2)
theta[1] = 1 # theta0
theta[2] = kappa^2 / (2*N)
theta[3] = theta[2] * (2 * N + kappa^2 / N)
for(i in 4:(N+1)){
   theta[i] = (2*N + kappa^2 / N) * theta[i-1] - N^2 * theta[i-2]
}
theta[N+2] = kappa^2 / (2*N) * theta[N + 1]
Binv2 = matrix(0, N+1, N+1)
for (i in 2:N){
   for(j in 2:N){
      if (i <= j){
         Binv2[i, j] = N^abs(j-i) * theta[i] * theta[N + 2 - j] / theta[N+2]  
      } else{
         Binv2[i, j] = N^abs(j-i) * theta[j] * theta[N + 2 - i] / theta[N+2] 
      }
   }
}
for (i in 1:(N+1)){
   if (i != 1){
      Binv2[i, N+1] = N^{N+1 - i} * theta[i]/theta[N+2]
   }
   if (i != N+1){
      Binv2[i, 1] = N^{i-1} * theta[N+2-i]/theta[N+2]
   }
}
BTB = t(B) %*% B
eigen(BTB)$values
## Property of Q
C = diag(N+1) / N
C[1, 1] = 1/2/N
C[N+1, N+1] = 1/2/N
Q1 = t(B) %*% solve(C) %*% B
min(diag(solve(Q1)))

###### pdf plotting
M = 100
w0 = seq(-epsilon, epsilon, length = M)
w1 = seq(-epsilon, epsilon, length = M)
w2 = seq(-epsilon, epsilon, length = M)





tmp = Binv2 %*% t(Binv2) / (2*N)
solve(Q1) - tmp
# solve(Q1)
# t(Binv2)
epsilon = 0.1
wvec = runif(N+1, -epsilon, epsilon)
# wvec = rep(epsilon, N+1)
# sum(wvec * (abs(Q1) %*% wvec))
# 32 * N^4 * epsilon^2
library(mvtnorm)
b = rmvnorm(100, sigma = solve(Q1))
var(apply(b, 1, max))
## MC estimate for the lower bound of exp(-w^T Q w / 2)
M = 1000
a = replicate(M, {
   wvec = runif(N+1, -epsilon, epsilon)
   exp(-sum(wvec * (abs(Q1) %*% wvec)) / 2)
})
mean(a)



## RBF Sigma
kappa = 2
cov_GPI = matrix(0, N+1, N+1)
for(i in 1:(N+1)){
   for (j in 1:(N+1)){
      cov_GPI[i, j] = exp(-kappa^2 * (i-j)^2 / N^2)
   }
}
max(eigen(solve(cov_GPI))$values)

