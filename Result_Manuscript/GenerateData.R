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
   X = matrix(runif(2*n), n)
   Z = f0(X[, 1], X[, 2]) + rnorm(n) * 0.1
   df = data.frame(X, Z)
   filename = paste0("Result_Manuscript/obs_n", n, ".RData")
   save(df, file = filename)
}