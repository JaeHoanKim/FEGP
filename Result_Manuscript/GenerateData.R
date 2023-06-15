## generate the observation for both methods
## 50 sets of data
nlist = c(200, 500, 1000)
f0 = function(x, y){
   return(sin(5*x + 2*y) + 2*y^2)
}
M = 50
for(i in 1:length(nlist)){
   set.seed(i)
   n = nlist[i]
   X = matrix(runif(2*n*M), n*M)
   Z = f0(X[, 1], X[, 2]) + rnorm(n*M) * 0.1
   df = data.frame(X, Z)
   filename = paste0("Result_Manuscript/obs_n2D", n, ".RData")
   save(df, file = filename)
}

f0_1D = function(x){return(1*x^2+sin(8*x))}

for(i in 1:length(nlist)){
   set.seed(i)
   n = nlist[i]
   X = runif(n*M)
   Z = f0_1D(X) + rnorm(n*M) * 0.1
   df = data.frame(X, Z)
   filename = paste0("Result_Manuscript/obs_n1D", n, ".RData")
   save(df, file = filename)
}

