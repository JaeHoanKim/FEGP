For the RData files in brn100 and 1000, the regression functions are as follows:

f0_1D = function(x){return (x^2 + sin(x))}
f0_2D = function(x, y){return(x^2 + sqrt(abs(y-0.5)) + sin(8*x))}