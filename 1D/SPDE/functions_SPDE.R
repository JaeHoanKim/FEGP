# fast multivariate sampling from the GMRF, conditioned on the data (Using ESS)
# generating the Q matrix

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

loglik_w = function(N, kappa, w, l=2){
   Q = Q1D(N, kappa)
   L = chol(Q)
   x = L %*% w
   return(sum(log(diag(L))) - (N+1)/2*log(2*pi) - 1/2* sum(x^2))
}

Sampling_N_new = function(N, kappa, beta = 2){
   Q = Q1D(N, kappa)
   L = chol(Q)
   x = rnorm(N+1)
   vec = solve(L, x)
   return(vec)
}

Phi_1D = function(X, N){
   n = length(X)
   knot_N = c(0:N)/N
   Phi = matrix(0, nrow = n, ncol = N+1)
   # for the sparse matrix
   for(index in 1:n){
      x = X[index]
      i = min(1 + floor(x*N), N)
      wx1 = (1 - abs(x - knot_N[i]) * N)
      wx2 = (1 - abs(x - knot_N[i+1]) * N)
      Phi[index, i] = wx1
      Phi[index, i+1] = wx2
   }
   Phi = as(Phi, "sparseMatrix")
   return(Phi)
}

f_N_h = function(x, g){
   N = length(g) - 1
   knot_N = c(0:N)/N
   i = min(1 + floor(x*N), N)
   value = (1 - abs(x - knot_N[i]) * N) * g[i] +
      (1 - abs(x - knot_N[i+1]) * N) * g[i+1]
   return(value)
}

f_N_h_vec = function(x, g){
   N = length(g) - 1
   knot_N = c(0:N)/N
   i = pmin(1 + floor(x*N), N)
   vec1 = (1 - abs(x - knot_N[i]) * N) * g[i]
   vec2 = (1 - abs(x - knot_N[i + 1]) * N) * g[i+1]
   return(vec1 + vec2)
}

loglik = function(y, x, g, sigsq){
   ## X: n by 1 vector, y: n by 1 vector
   n = length(y)
   f_N_x = f_N_h_vec(x, g)
   loglik_value = - sum((y - f_N_x)^2) / (2 * sigsq) - n/2 * log((2*pi) * sigsq) 
   return(loglik_value)
}

# when given two samples from the prior, return the posterior sample using ESS algorithm
ESS_post = function(y, x, g, g_ESS, sigsq){
   thetamin = 0; 
   thetamax = 2*pi;
   u = runif(1)
   logy = loglik(y, x, g, sigsq) + log(u); 
   theta = runif(1,thetamin,thetamax); 
   thetamin = theta - 2*pi 
   thetamax = theta
   gprime = g*cos(theta) + g_ESS*sin(theta)
   while(loglik(y, x, gprime, sigsq) <= logy){
      if(theta < 0)
         thetamin = theta
      else
         thetamax = theta
      theta = runif(1,thetamin,thetamax)
      gprime = g*cos(theta) + g_ESS*sin(theta)
   }
   return(gprime)    
}

## Tempered version of ESS - when g, g_ESS : sample from N(f; 0, \Sigma),
## this returns the sample from N(f; 0, \Sigma) {L(f)}^{1/Temp}
ESS_post_tempered = function(y, x, g, g_ESS, sigsq, Temp=1){
   thetamin = 0; 
   thetamax = 2*pi;
   u = runif(1)
   logy = loglik(y, x, g, sigsq)/Temp + log(u); 
   theta = runif(1,thetamin,thetamax); 
   thetamin = theta - 2*pi 
   thetamax = theta
   gprime = g*cos(theta) + g_ESS*sin(theta)
   while(loglik(y, x, gprime, sigsq)/Temp <= logy){
      if(theta < 0)
         thetamin = theta
      else
         thetamax = theta
      theta = runif(1,thetamin,thetamax)
      gprime = g*cos(theta) + g_ESS*sin(theta)
   }
   return(gprime)    
}




glist_to_plotdf = function(g, grid, true, alpha1 = 0.95, alpha2 = 0.9){
   Nsample = length(g)
   N = length(grid)
   f0 = true
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
   colnames(y.plot) <- c('x', 'low2', 'low1', 'mean', 'upp1', 'upp2', 'true', 'lowsup', 'uppsup')
   return(y.plot)
}
