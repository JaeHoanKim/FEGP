N = 20
kappa = 0.2
beta = 6
r1 = 1 + kappa^2/(2*N^2) + kappa/N * sqrt(1 + kappa^2 / 4 / N^2) 
r2 = 1 + kappa^2/(2*N^2) - kappa/N * sqrt(1 + kappa^2 / 4 / N^2) 
a = vector(length = N)
for (i in 1:N){
   a[i] = (r1^i - r2^i)/(r1-r2)
}
X = matrix(0, N-1, N-1)
for(i in 1:(N-1)){
   for(j in 1:(i-1)){
      X[i, j] = a[j] * a[N-i]/a[N]
   }
   for(j in i:(N-1)){
      X[i, j] = a[i] * a[N-j]/a[N]
   }
}
Z = cbind(rev(a[1:(N-1)]), a[1:(N-1)])/a[N]

## Post mean
G = matrix(0, N+1, N+1)
for (i in 2:N){
   G[i, i] = 2 * N
   G[i, i-1] = -N
   G[i, i+1] = -N
}
C = diag(N+1) / N
C[1, 1] = 1/2/N
C[N+1, N+1] = 1/2/N
A = kappa^2 * diag(N+1) + N * G
Ainv = solve(A)
Varbeta = solve(C)
for(i in 1:(beta/2)){
   Varbeta = Ainv %*% Varbeta %*% t(Ainv)
}
Qbeta = C
for(i in 1:(beta/2)){
   Qbeta = t(A) %*% Qbeta %*% A
}
solve(Qbeta)[1, 1]
idx_end = c(1, N+1)
idx_mid = c(1:(N+1))[-idx_end]


idx_end2 = sort(1 + c(0:(beta/2), N - c(0:(beta/2))))
idx_mid2 = c(1:(N+1))[-idx_end2]

a11 = -solve(Qbeta[idx_mid2, idx_mid2]) %*% Qbeta[idx_mid2, idx_end2]
a11
# conditional mean form exact calculation
a1 = -solve(Qbeta[idx_mid, idx_mid]) %*% Qbeta[idx_mid, idx_end]
# conditional mean from our formula 
a2 = Z
a2 = Z + kappa^2 / N^2 * X %*% Z + kappa^4 / N^4 * X %*% X %*% Z
max(abs(a1-a2))

Y = X %*% t(X)
# plot(Y[N/2, ])

for(i in 1:(beta/2 - 1)){
   Y = X %*% Y %*% t(X) 
}
plot(Y[N/2, ])
Y = Y/N^(2*beta - 1)
a = vector(length = N-2)
for(i in c(1:(N-2))){
   a[i] = Y[i+1, i+1] + Y[i, i] - 2 * Y[i, i+1]
}

plot(a)

for(i in c((N/4+1):(N/2))){
   b[i-N/4] = Y[i+1, i+1] + Y[N/4, N/4] - 2 * Y[i+1, N/4]
}

plot(b)

1/2^(4 * beta - 3) / exp(3 * kappa*(beta - 1))

## Correlation for Y
dYinvhalf = diag(1/sqrt(diag(Y)))
Corr_Y = dYinvhalf %*% Y %*% dYinvhalf
min(Corr_Y)

## Check the explicit formula

G = matrix(0, N+1, N+1)
for (i in 2:N){
   G[i, i] = 2 * N
   G[i, i-1] = -N
   G[i, i+1] = -N
}
C = diag(N+1) / N
C[1, 1] = 1/2/N
C[N+1, N+1] = 1/2/N
A = kappa^2 * diag(N+1) + N * G
Ainv = solve(A)
Varbeta = solve(C)
for(i in 1:(beta/2)){
   Varbeta = Ainv %*% Varbeta %*% t(Ainv)
}
Qbeta = C
for(i in 1:(beta/2)){
   Qbeta = t(A) %*% Qbeta %*% A
}
idx1 = c(c(0:beta), c((N-beta):N))
idx2 = c(0:N)[-(idx1+1)]
Qbeta[idx1, idx2]
idx_center = 1 + (N/2-beta/2):(N/2+beta/2)
idx_end = 1 + c((0:(beta/2)), (N-beta/2):(N))
idx_not_end = c(1:(N+1))[-idx_end]

Qbeta[idx_not_end, idx_end]

Q11 = Qbeta[idx_not_end, idx_not_end]
Q12 = Qbeta[idx_not_end, idx_end]
solve(Q11) %*% Q12



solve(Qbeta[idx])


Qbeta[(beta/2+1):(N-beta/2-1), ]
Tmat = matrix(0, N-beta+1, N-beta+1)
for(i in c(1:(N-beta+1))){
   for(j in c(1:(N-beta+1))){
      if (i==j){
         Tmat[i, j] = 2+kappa^2/N^2
      }
      if(abs(i-j) == 1){
         Tmat[i, j] = -1
      }
   }
}
Tmatbeta = Tmat
Tmatbeta_inv = solve(Tmat)
for(j in 1:(beta - 1)){
   Tmatbeta = Tmatbeta %*% Tmat
   Tmatbeta_inv = Tmatbeta_inv %*% solve(Tmat)
}
Tfin = N^(2*beta - 1) * Tmatbeta

Tmatbeta_inv %*% Qbeta[idx_not_end, idx_end]



idx_range = (beta/2+1):(N-beta/2 + 1)
Qbeta0 = Qbeta[idx_range, idx_range]



Tfin
Qbeta0
evals = eigen(Qbeta0)$values
evals = evals / N^(2*beta - 1)
plot(1/evals^(1/beta))
solve(Qbeta)[c(1, N+1), c(1, N+1)]


solve(Qbeta)[3:(N-1), 3:(N-1)]

Varbeta[3:(N-1), 3:(N-1)]
max(abs(solve(Qbeta)[2:N, 2:N] - Varbeta[2:N, 2:N]))


solve(Qbeta[(beta/2+1):(N-beta/2), (beta/2+1):(N-beta/2)])

Y[(beta/2):(N-beta/2-1), (beta/2):(N-beta/2-1)]
2*N*(1/kappa^4)^(beta/2)

-solve(Qbeta[2:N, 2:N]) %*% Qbeta[2:N, c(1, N+1)]

Z + kappa^2 / N^2 * X %*% Z + kappa^4 / N^4 * X %*% X %*% Z



## Cofficient for the conditional expectation


N^2 / (4 * kappa^2)

first = 1
mid = 2:N
last = N+1
# Well-coded for var(w_{1:N-1} | w_{0, N})
mat22 = rbind(cbind(Varbeta[first, first], Varbeta[first, last]), cbind(Varbeta[last, first], Varbeta[last, last]))
tmp = Varbeta[mid, mid] -  Varbeta[mid, c(first, last)] %*% solve(mat22) %*%  Varbeta[c(first, last), mid]
tmp[(nrow(tmp) + 1)/2, (ncol(tmp) + 1)/2] # Var(w_{N/2} | w_0, w_N)

# Calculating var(w_{1:N-1} | w_{0, N}, w_{1:beta/2, N-beta/2:N-1})

first = 1:(beta/2)
mid = (beta/2 + 1):(N - beta/2-1)
last = (N - beta/2):(N-1)

mat22 = rbind(cbind(tmp[first, first], tmp[first, last]), cbind(tmp[last, first], tmp[last, last]))
tmp2 = tmp[mid, mid] -  tmp[mid, c(first, last)] %*% solve(mat22) %*%  tmp[c(first, last), mid]
tmp2[(nrow(tmp2) + 1)/2, (ncol(tmp2) + 1)/2] # Var(w_{N/2} | w_0, w_N)





## Check
tmp[1, 1]
1/N^3/a[N]^2 * sum(a[1:(N-1)]^2)
tmp[(nrow(tmp) + 1)/2, (ncol(tmp) + 1)/2]

X = X / N^
beta = 4
tmp = X
for(i in 1:(beta-1)){
   tmp = tmp %*% X
}
tmp = tmp * N
tmp
# tmp = tmp / (2 * N)^{2*beta - 1}
# tmp
y = tmp[N/2 + 1, ]
x = tmp[c(1:(beta /2), (N-beta/2):(N-1)), ]
df = data.frame(y, t(x))
model = lm(y~. -1, data =df)
predicted_values <- predict(model)
# Calculate the sum of squared errors (SSE)
sse <- sum((y - predicted_values)^2)
print(sse)
