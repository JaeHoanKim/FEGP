True function


1D
1: f0_1D = function(x, trun = 500){
   value = 0
   for(j in 1:trun){
      value = value + sin(j) * cos(pi * (j - 1/2) * x) * j^(- alpha - 1)
   }
   return(value * sqrt(2))
}
2D

1: f0_2D = function(x, y){return(x^2 + sqrt(abs(y-0.5)) + sin(8*x))}
2: f0_2D = function(x, y){return(sin(5*x + 2*y) + 2*y^2)}