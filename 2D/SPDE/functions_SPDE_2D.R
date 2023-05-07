# fast multivariate sampling from the GMRF, conditioned on the data (Using ESS)

Q1D = function(N, kappa, beta = 2){
   # generating C as a sparse matrix (N=1 by N+1 matrix)
   i1 = j1 = c(1:(N+1))
   x = 1 / (2*N) * c(1, rep(2, N-1), 1)
   C = sparseMatrix(i = i1, j = j1, x = x)
   i2 = rep(c(2:N), each = 3)
   # generating G as a sparse matrix (N=1 by N+1 matrix)
   j2 = i2 + c(-1, 0, 1)
   x2 = N*rep(c(-1, 2, -1), times = N-1)   
   G = sparseMatrix(i = i2, j = j2, x = x2, dims = c(N+1, N+1))
   if (beta == 2){
      Q = kappa^4*C + kappa^2 * (G+t(G)) + N*t(G)%*%G
   }
   return(Q)
}


Q2D = function(N, kappa, beta = 2){
   # generating C as a sparse matrix (N=1 by N+1 matrix)
   i1 = j1 = c(1:(N+1))
   x = 1 / (2*N) * c(1, rep(2, N-1), 1)
   C = sparseMatrix(i = i1, j = j1, x = x)
   i2 = rep(c(2:N), each = 3)
   # generating G as a sparse matrix (N=1 by N+1 matrix)
   j2 = i2 + c(-1, 0, 1)
   x2 = N*rep(c(-1, 2, -1), times = N-1)   
   G = sparseMatrix(i = i2, j = j2, x = x2, dims = c(N+1, N+1))
   I = sparseMatrix(i = i1, j = j1, x = rep(1, N+1))
   A = kappa^2 * kronecker(I, I) + N * kronecker(I, G) + N * kronecker(G, I)
   Q = crossprod(A, kronecker(C, C)) %*% A
   return(Q)
}

# when the vector w is sampled from N(0, Q), calculate its log likelihood in 2D, input should be (N+1)^2 long vector
loglik_w2D = function(N, kappa, w){
   Q = Q2D(N, kappa)
   L = chol(Q)
   x = L %*% w
   return(sum(log(diag(L))) - (N+1)^2/2*log(2*pi) - 1/2* sum(x^2))
}


# sampling from N(0, Q^{-1}) using the fast Cholesky decomposition in 2D - returns (N+1)^2 vector
Sampling_N_new2D = function(N, kappa, beta = 2){
   Q = Q2D(N, kappa)
   L = chol(Q)
   x = rnorm((N+1)^2)
   vec = solve(L, x)
   return(vec)
}


f_N_h_2D_multi = function(x, g){
   n = dim(x)[1]
   out = vector(length = n)
   N1 = dim(g)[1] - 1
   N2 = dim(g)[2] - 1
   knot_N1 = c(0:N1)/N1
   knot_N2 = c(0:N2)/N2
   for(index in 1:n){
      ## x, y coordinate of index'th point
      xx = x[index, 1]
      yy = x[index, 2]
      i = min(1 + floor(xx * N1), N1)
      j = min(1 + floor(yy * N2), N2)
      wx1 = (1 - abs(xx - knot_N1[i]) * N1)
      wx2 = (1 - abs(xx - knot_N1[i+1]) * N1)
      wy1 = (1 - abs(yy - knot_N2[j]) * N2)
      wy2 = (1 - abs(yy - knot_N2[j+1]) * N2)
      value =  wx1 * wy1 * g[i, j] +
         wx1 * wy2 * g[i, j+1] +
         wx2 * wy1 * g[i+1, j] +
         wx2 * wy2 * g[i+1, j+1]
      out[index] = value 
   }
   return(out)
}


loglik2D = function(z, x, g, sigsq){
   f_N_x = f_N_h_2D_multi(x, g) # n by 1 vector when x is n by 2 vector, g is a square matrix
   loglik_value = - sum((z - f_N_x)^2) / (2 * sigsq) - n / 2 *log((2*pi) * sigsq)
   return(loglik_value)
}

# run ESS algorithm

ESS_post2D = function(z, x, g, g_ESS, sigsq){
   thetamin = 0; 
   thetamax = 2*pi;
   u = runif(1)
   logy = loglik2D(z, x, g, sigsq) + log(u); 
   theta = runif(1,thetamin,thetamax); 
   thetamin = theta - 2*pi 
   thetamax = theta
   gprime = g*cos(theta) + g_ESS*sin(theta)
   while(loglik2D(z, x, gprime, sigsq) <= logy){
      if(theta < 0)
         thetamin = theta
      else
         thetamax = theta
      theta = runif(1,thetamin,thetamax)
      gprime = g*cos(theta) + g_ESS*sin(theta)
   }
   return(gprime)
}
