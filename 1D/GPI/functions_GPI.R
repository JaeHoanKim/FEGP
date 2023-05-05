library(fields)
library(FastGP)

h_j = function(x, my_knot, delta_N){
   N = length(my_knot) - 1 # 0 / N to N / N
   i = pmin(1 + floor(x * N), N)
   k = rep(0, N+1) # store the value of h_0(x) to h_N(x)
   k[i] = 1 - abs(x - my_knot[i]) / delta_N
   k[i+1] = 1 - abs(x - my_knot[i+1]) / delta_N
   return(k)
}


# f_N - with tent basis
f_N_h = function(x, g){
   N = length(g) - 1
   knot_N = c(0:N)/N
   i = min(1 + floor(x*N), N)
   value = (1 - abs(x - knot_N[i]) * N) * g[i] +
      (1 - abs(x - knot_N[i+1]) * N) * g[i+1]
   return(value)
}

## returning f_N_h for vector x
f_N_h_vec = function(x, g){
   N = length(g) - 1
   knot_N = c(0:N)/N
   i = pmin(1 + floor(x*N), N)
   vec1 = (1 - abs(x - knot_N[i]) * N) * g[i]
   vec2 = (1 - abs(x - knot_N[i + 1]) * N) * g[i+1]
   return(vec1 + vec2)
}

#Given a \nu (smoothness parameter of matern kernel) finding a value of 
# l (length-scale parameter) such that the correlation between the 
# maximum seperation is some small value, say 0.05

## Matern kernel with smoothness nu and length-scale l:
MK = function(x, y ,l, nu){
   ifelse(abs(x-y)>0, (sqrt(2*nu)*abs(x-y)/l)^nu/(2^(nu-1)*gamma(nu))*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu), 1.0)
}



# Covariance matrix
covmat=function(knot,nu,l){
   return(MK(rdist(knot), 0, l, nu))
}


# Order of the circulant matrix:
# minimum value of g and m so that G can be embedded into C
min_m=function(knot){
   N=length(knot)
   g=ceiling(log(2*N,2))   #m=2^g and m>=2(n-1) : Wood & Chan notation; 
   m = 2^g
   #since we are going upto n and not stopping at (n-1), the condition is modified!
   return("m" = m)
}

# forming the circulant matrix:
circulant=function(x){
   n = length(x)
   mat = matrix(0, n, n)
   for (j in 1:n) {
      mat[j, ] <- c(x[-(1:(n+1-j))], x[1:(n+1-j)])
   }
   return(mat)
}

# Function for forming the vector of circulant matrix - MK:
circ_vec=function(knot,m,nu,l,tausq){
   delta_N=1/(length(knot)-1)
   cj=integer()
   for(j in 1:m){
      if(j<=(m/2))
         cj[j]=(j-1)*delta_N
      else
         cj[j]=(m-(j-1))*delta_N
   }
   x=(tausq*MK(cj,0,l,nu))
   return(x)
}

# Function for finding a m such that C is nnd:
eig.eval=function(knot,g,nu,l,tausq){
   vec=circ_vec(knot,g,nu,l,tausq)
   C=circulant(vec)
   ev=min(eigen(C)$values)
   return(list("vec" = vec, "min.eig.val" = ev))
}

# Function for finding a m such that C is nnd:
# without forming the circulant matrix and without computing eigen values:
C.eval=function(knot, m, nu, l, tausq){
   vec=circ_vec(knot, m, nu, l, tausq)
   val=fft(vec) # eigenvalues will be real as the circulant matrix formed by the 
   # vector is by construction is symmetric!
   ev=min(Re(val))
   return(list("vec" = vec, "min.eig.val" = ev))
}

nnd_C=function(knot,m,nu,l,tausq){
   C.vec=C.eval(knot,m,nu,l,tausq)$vec
   eval=C.eval(knot,m,nu,l,tausq)$min.eig.val
   if(eval>0)
      return(list("cj" = C.vec,"m" = m))
   else{
      m = m + 2
      nnd_C(knot,m,nu,l,tausq)
   }
}

# computing the eigen values of C using FFT:
eigval=function(knot,nu,l,tausq){
   m = min_m(knot)
   c.j=nnd_C(knot,m,nu,l,tausq)$cj
   lambda=Re(fft(c.j))
   if(min(lambda)>0)
      return(lambda)
   else
      stop("nnd condition is NOT satisfied!!!")
}

#################################################################
########## Samples drawn using Wood and Chan Algorithm ##########
#################################################################

samp.WC=function(knot,nu,l,tausq){
   N=length(knot)
   lambda=eigval(knot,nu,l,tausq)
   m=length(lambda)
   samp.vec=rep(0,N)
   a=rep(0,m)
   a[1]=sqrt(lambda[1])*rnorm(1)/sqrt(m)
   a[(m/2)+1]=sqrt(lambda[(m/2)+1])*rnorm(1)/sqrt(m)
   i=sqrt(as.complex(-1))
   for(j in 2:(m/2)){
      uj=rnorm(1); vj=rnorm(1)
      a[j]=(sqrt(lambda[j])*(uj + i*vj))/(sqrt(2*m))
      a[m+2-j]=(sqrt(lambda[j])*(uj - i*vj))/(sqrt(2*m))
   }
   samp=fft(a)
   samp.vec=Re(samp[1:N])
   return(samp.vec)
}



## Defining the loglik function to be used in ESS:
## loglik calculates the log of the likelihood:

loglik = function(y, x, g, sigsq){
   ## X: n by 1 vector, y: n by 1 vector
   n = length(y)
   f_N_x = vector(length = n)
   for (i in 1:n){
      f_N_x[i] = f_N_h(x[i], g)
   }
   loglik_value = - sum((y - f_N_x)^2) / (2 * sigsq) - n/2 * log((2*pi) * sigsq) 
   return(loglik_value)
}


loglik_v = function(v, g, nu, l){
   if (length(v) + 1 != length(g)){
      stop("Error: the length of v should be equal to length(g) - 1")
   }
   N = length(v)
   knot_N = c(0:N)/N
   knot_2N = c(0:(2*N))/(2*N)
   Sigma_N = covmat(knot_N, nu, l)[1:N, 1:N]
   Sigma_2N = covmat(knot_2N, nu, l)
   S_N = inv_chol(Sigma_N)
   S_2N = inv_chol(Sigma_2N)
   gv = knit_vec(g, v)
   eta_N = t(S_N) %*% v # N vector
   eta_2N = t(S_2N) %*% gv # 2N+1 vector
   result = -0.5 * sum(eta_2N^2) + 0.5 * sum(eta_N^2) # + sum(log(diag(S_2N))) - sum(log(diag(S_N)))  
   return(result)
}

knit_vec = function(a, b){
   if (length(a) != length(b)+1){
      stop("error: length(a) should be equal to 1 + length(b)!")
   }
   else{
      N = length(b)
      vec = rep(0, 2 * N + 1)
      vec[2 * c(1:N)] = b
      vec[(2 * c(1:(N+1))) - 1] = a
   }
   return(vec)
}

#############################################
########## Functions for using ESS ##########
#############################################
ESS = function(g, nu_ess, y, x, sigsq){
   thetamin = 0; 
   thetamax = 2*pi;
   
   u = runif(1)
   logy = loglik(y, x, g, sigsq) + log(u); 
   
   theta = runif(1,thetamin,thetamax); 
   thetamin = theta - 2*pi; 
   thetamax = theta;
   gprime = g*cos(theta) + nu_ess*sin(theta);
   
   while(loglik(y, x, gprime, sigsq) <= logy){
      if(theta < 0)
         thetamin = theta
      else
         thetamax = theta
      theta = runif(1,thetamin,thetamax)
      gprime = g*cos(theta) + nu_ess*sin(theta)
   }
   return(gprime)       
}

v_g_ESS = function(v, g, nu_ess, nu, l){
   thetamin = 0; 
   thetamax = 2*pi;
   u = runif(1)
   logy = loglik_v(v, g, nu, l) + log(u); 
   
   theta = runif(1,thetamin,thetamax); 
   thetamin = theta - 2*pi; 
   thetamax = theta;
   vprime = v*cos(theta) + nu_ess*sin(theta);
   
   while(loglik_v(vprime, g, nu, l) <= logy){
      if(theta < 0)
         thetamin = theta
      else
         thetamax = theta
      theta = runif(1,thetamin,thetamax)
      vprime = v*cos(theta) + nu_ess*sin(theta)
   }
   return(vprime)
}

# 1D Levinson algorithm
Levinson = function (gamma, y){
   gamma_mat = circulant(gamma)
   ## direct solution  
   sol1 = solve(gamma_mat, y)
   ## Levinson algorithm
   M = length(gamma) - 2 # suppose M = 3
   if (M <= -1){
      stop("gamma should have length of at least 2!")
   }
   beta = gamma[2] / gamma[1] # initial vector
   alpha = y[1] / gamma[1]
   alpha_end = (y[2] - alpha * gamma[2]) / (gamma[1] - beta * gamma[2])
   alpha = c((alpha - beta * alpha_end), alpha_end) 
   # one more iteration
   if (M >= 1){
      for(m in 1:M){
         # at m stage, beta should be m+1 length vector & alpha: m+2
         beta0 = (gamma[m+2] - sum(beta * gamma[2:(m+1)])) / (gamma[1] - sum(beta * gamma[(m+1):2]))
         beta = c(beta0, (beta - beta0 * beta[m:1])) # length of m+1
         alpha_end = (y[m+2] - sum(alpha * gamma[(m+2):2])) / (gamma[1] - sum(beta * gamma[(m+2):2])) 
         alpha = c((alpha - beta * alpha_end), alpha_end)
      } 
   }
   return(list(sol1, alpha))
}



glist_to_plotdf = function(g, grid, true = f0, alpha1 = 0.95, alpha2 = 0.9){
   Nsample = length(g)
   N = length(grid)
   y.tot = matrix(nrow = Nsample, ncol = N)
   for(i in 1:Nsample){
      y.tot[i, ] = f_N_h_vec(grid, g[[i]])
   }
   y.plot = matrix(nrow = length(grid), ncol = 9)
   y.plot[, 1] = grid
   for(i in 1:N){
      y.plot[i, 2:6] = quantile(y.tot[, i], probs = c((1-alpha2)/2, (1-alpha1)/2, 0.5, (1+alpha1)/2, (1+alpha2)/2))
   }
   ## replace median to mean
   y.plot[, 4] = colMeans(y.tot)
   y.plot[, 7] = f0(grid)
   ## add sup norm ribbon bands
   diff.tot = t(abs(t(y.tot) - y.plot[, 4])) # difference with the posterior mean at each gridpoints
   index = max.col(diff.tot)
   supnorms = vector(length = Nsample)
   for(i in 1:Nsample){
      supnorms[i] = diff.tot[i, index[i]]
   }
   supband = quantile(supnorms, probs = alpha1)
   y.plot[, 8] = y.plot[, 4] - supband
   y.plot[, 9] = y.plot[, 4] + supband
   y.plot = data.frame(y.plot)
   colnames(y.plot) <- c('x', 'low2', 'low1', 'med', 'upp1', 'upp2', 'true', 'lowsup', 'uppsup')
   return(y.plot)
}