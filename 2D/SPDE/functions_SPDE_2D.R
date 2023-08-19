# fast multivariate sampling from the GMRF, conditioned on the data (Using ESS)

# Q1D would not be used
Q1D = function(N, kappa, beta, tausq){
   # generating C as a sparse matrix (N=1 by N+1 matrix)
   i1 = j1 = c(1:(N+1))
   x = 1 / (2*N) * c(1, rep(2, N-1), 1)
   C = sparseMatrix(i = i1, j = j1, x = x)
   invC = sparseMatrix(i = i1, j = j1, x = 1/x)
   i2 = rep(c(2:N), each = 3)
   # generating G as a sparse matrix (N=1 by N+1 matrix)
   j2 = i2 + c(-1, 0, 1)
   x2 = N*rep(c(-1, 2, -1), times = N-1)   
   G = sparseMatrix(i = i2, j = j2, x = x2, dims = c(N+1, N+1))
   A = kappa^2 * C + G
   nu = beta - 1/2 # nu in 1D
   tau0sq = gamma(nu)/gamma(nu+1/2)/sqrt(4*pi)/kappa^{2*nu}
   if (beta == 2){
      Q = t(A) %*% invC %*% A
      return(Q * tau0sq / tausq)
   }
   else if (beta %% 2 == 0){
      Q = t(A) %*% invC %*% Q1D(N, kappa, beta - 2, tausq) %*% invC %*% A
      return(Q) # constant is already multiplied when returned from Q1D function
   }
}


Q2D = function(N, kappa, tausq, beta = 2){
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
   nu = beta - 1 # nu in 2D
   tau0sq = gamma(nu)/gamma(nu+1/2)/sqrt(4*pi)/kappa^{2*nu}
   return(Q * tau0sq/tausq)
   
}

# when w is sampled from N(0, Q) in 1D case
loglik_w1D = function(N, kappa, w, l=2){
   Q = Q1D(N, kappa)
   L = chol(Q)
   x = L %*% w
   return(sum(log(diag(L))) - (N+1)/2*log(2*pi) - 1/2* sum(x^2))
}

# sampling from N(0, Q^{-1}) using the fast Cholesky decomposition
Sampling_N_new1D = function(N, kappa, beta = 2){
   Q = Q1D(N, kappa)
   L = chol(Q)
   x = rnorm(N+1)
   vec = solve(L, x)
   return(vec)
}


# when the vector w is sampled from N(0, Q), calculate its log likelihood in 2D, input should be (N+1)^2 long vector
loglik_w2D = function(N, kappa, tausq, w){
   Q = Q2D(N, kappa, tausq)
   L = chol(Q)
   x = L %*% w
   return(sum(log(diag(L))) - (N+1)^2/2*log(2*pi) - 1/2* sum(x^2))
}


# sampling from N(0, Q^{-1}) using the fast Cholesky decomposition in 2D - returns (N+1)^2 vector
Sampling_N_new2D = function(N, kappa, tausq, beta = 2){
   Q = Q2D(N, kappa, tausq)
   L = chol(Q)
   x = rnorm((N+1)^2)
   vec = solve(L, x)
   return(vec)
}

# tent basis function

# generating Phi matrix
# X: n by 2 matrix
# N: number of grid - 1
# output: n by N matrix
Phi_2D = function(X, N){
   n = nrow(X)
   gridmat = cbind(rep(c(0:N)/N, each = N + 1),
                   rep(c(0:N)/N, N + 1))
   knot_N = c(0:N)/N
   Phi = matrix(0, nrow = n, ncol = (N+1)^2)
   # for the sparse matrix
   ilist = c()
   jlist = c()
   plist = c()
   for(index in 1:n){
      Phi_row = matrix(0, nrow = N+1, ncol = N+1)
      xx = X[index, 1]
      yy = X[index, 2]
      i = min(1 + floor(xx * N), N)
      j = min(1 + floor(yy * N), N)
      wx1 = (1 - abs(xx - knot_N[i]) * N)
      wx2 = (1 - abs(xx - knot_N[i+1]) * N)
      wy1 = (1 - abs(yy - knot_N[j]) * N)
      wy2 = (1 - abs(yy - knot_N[j+1]) * N)
      Phi_row[i, j] = wx1 * wy1
      Phi_row[i, j+1] = wx1 * wy2
      Phi_row[i+1, j] = wx2 * wy1
      Phi_row[i+1, j+1] = wx2 * wy2
      # instead of byrow = T option
      Phi[index, ] = as.vector(t(Phi_row))
   }
   Phi = as(Phi, "sparseMatrix")
   return(Phi)
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

ESS_post2D = function(z, x, g, g_ESS, sigsq, u = 0.5, Temper = 1){
   thetamin = 0; 
   thetamax = 2*pi;
   # to control the randomness within the ESS algorithm
   if(u == 0.5){
      u = runif(1)  
   }
   logy = loglik2D(z, x, g, sigsq) / Temper + log(u); 
   theta = runif(1,thetamin,thetamax); 
   thetamin = theta - 2*pi 
   thetamax = theta
   gprime = g*cos(theta) + g_ESS*sin(theta)
   while(loglik2D(z, x, gprime, sigsq) / Temper <= logy){
      if(theta < 0)
         thetamin = theta
      else
         thetamax = theta
      theta = runif(1,thetamin,thetamax)
      gprime = g*cos(theta) + g_ESS*sin(theta)
   }
   return(gprime)
}



# gridmat is a (gridsize^2) by 2 matrix!
# input: g_N value ((N+1)^2 by 2 matrix) and grid matrix at which we want to obtain f_N value
# ouput: matrices of f_N values - low1, upp1, true, lowsup, uppsup, means
glist_to_plotdf_2D = function(g, gridmat, truefun, alpha1 = 0.95, alpha2 = 0.9){
   Nsample = length(g)
   N = nrow(gridmat)
   f0 = truefun
   y.tot = matrix(nrow = Nsample, ncol = N)
   for(i in 1:Nsample){
      gi.mat = matrix(g[[i]], nrow = sqrt(length(g[[i]])), byrow = TRUE)
      y.tot[i, ] = f_N_h_2D_multi(gridmat, gi.mat) 
   }
   y.plot = matrix(nrow = nrow(gridmat), ncol = 10)
   y.plot[, c(1, 2)] = gridmat
   for(i in 1:N){
      y.plot[i, 3:7] = quantile(y.tot[, i], probs = c((1-alpha2)/2, (1-alpha1)/2, 0.5, (1+alpha1)/2, (1+alpha2)/2))
   }
   ## replace median to mean
   y.plot[, 5] = colMeans(y.tot)
   y.plot[, 8] = f0(gridmat[, 1], gridmat[, 2])
   ## add sup norm ribbon bands
   diff.tot = t(abs(t(y.tot) - y.plot[, 5])) # difference with the posterior mean at each gridpoints
   index = max.col(diff.tot)
   supnorms = vector(length = Nsample)
   for(i in 1:Nsample){
      supnorms[i] = diff.tot[i, index[i]]
   }
   supband = quantile(supnorms, probs = alpha1)
   y.plot[, 9] = y.plot[, 5] - supband
   y.plot[, 10] = y.plot[, 5] + supband
   y.plot = data.frame(y.plot)
   colnames(y.plot) <- c('x1', 'x2', 'low2', 'low1', 'mean', 'upp1', 'upp2', 'truefun', 'lowsup', 'uppsup')
   return(y.plot)
}