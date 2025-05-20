M = 1234
x = c(1:M)/M
N = 10000
alpha = 10
f = function(x, alpha){
   return(x^alpha)
}


f_N_vec = function(x, g){
   N = length(g) - 1
   knot_N = c(0:N)/N
   i = pmin(1 + floor(x*N), N)
   vec1 = (1 - abs(x - knot_N[i]) * N) * g[i]
   vec2 = (1 - abs(x - knot_N[i + 1]) * N) * g[i+1]
   return(vec1 + vec2)
}

w = f(c(0:N)/N, alpha)
y = f(x, alpha)
yhat = f_N_vec(x, w)
plot(x, y)
plot(x, yhat)
plot(x, y - yhat, type = "l")
print(max(abs(y - yhat)))
   
