True function


1D
1: f0_1D = function(x, trun = 500){
   value = 0
   for(j in 1:trun){
      value = value + sin(j) * cos(pi * (j - 1/2) * x) * j^(- alpha - 1)
   }
   return(value * sqrt(2))
}
2: f0_1D = function(x){return(x^2 + sin(x))}

2D

1: f0_2D = function(x, y){return(x^2 + sqrt(abs(y-0.5)) + sin(8*x))}
2: f0_2D = function(x, y){return(sin(5*x + 2*y) + 2*y^2)}

################ MSE Setting #################
1D:

kappak = seq(1, 5, 0.5)
tausqk = 1
Nk = c(4, 6, 8, 10, 12)
kappa.pr = tausq.pr = N.pr = const 
beta = 4

2D:

beta = 2
d = 2
nu = beta - d/2

const = function(x){
   return(1)
}
Nk = c(6, 8, 10, 14, 18, 22, 26, 30)
N.pr = kappa.pr = tausq.pr = const
tausq.pr = function(x){return(invgamma::dinvgamma(x, 1, 1))}
kappa.pr = function(x){return(1/x^2)}
kappak = seq(1, 6, 0.5)
tausqk = 1
