### Amir Nikooienejad, Jan 2018
library(fields)
#================= Necessary Functions =================

powerset <- function(s){
   # Producing all possible subsets of c(1:n) set. It returns a list with 2^n elements.
   # Note than the i^th element of that list has a subset whose complement is in the ((2^n)-i)^th element.
   len <- length(s)
   l <- vector(mode="list",length=2^len) ; l[[1]] <- numeric()
   counter <- 1L
   for(x in 1L:length(s)){
      for(subset in 1L:counter){
         counter <- counter+1L
         l[[counter]] <- c(l[[subset]],s[x])
      }
   }
   return(l)
}
#==================

f1 <- function(X) {
   ## The first objective function in simulations
   x1 <- X[,1]; x2 <- X[,2]
   out <- 2*sin(30*(x1^4+x2^2)) + 2*cos(15*(x1-x2)) + sin(30*(x1-x2^2)) + cos(15*(x1^3+x2^2)) + 2*(1-x1^2) + 2*x2
   return(out/4)
}
#==================

f2 <- function(X) { ## The second objective function in simulations
   x1 <- X[,1]; x2 <- X[,2]
   out <- 2*sin(30*(x1^4+x2^2)) + 2*cos(15*(x1-x2)) + sin(30*(x1-x2^2)) + cos(15*(x1^3+x2^2)) + 2*(1-x1^2) + 2*x2 -
      2*(x1+x2)*(0.2 < x1 & x1 < 0.6 & 0.6 < x2 & x2 < 1) +
      2*(x1+x2)*(0.5 < x1 & x1 < 0.9 & 0.1 < x2 & x2 < 0.5)
   return(out/4)
}
#==================

y1 <- function(X,sigma2){ ## The observations based on the first objective function
   n <- dim(X)[1]
   out <- f1(X) + sqrt(sigma2) * rnorm(n)
   return(out)
}
#==================

y2 <- function(X,sigma2){ ## The observations based on the second objective function
   n <- dim(X)[1]
   out <- f2(X) + sqrt(sigma2) * rnorm(n)
   return(out)
}
#==================

thekernel <- function(X,nu,lambda){
   ## Matern kernel with parameters nu and lambda for range
   ## X could be either the matrix for point coordinates with each point in the rows or
   ## X could be a vector where it is considered as the difference at each coordinate. Thus, sqrt(sum(x^2)) is used as the input to Matern.
   d <- dim(X)
   if(!is.null(d)){
      X <- rdist(X)
   } else {X <- sqrt(sum(X^2))}
   out <- Matern(X,nu=nu,range=lambda)
   return(out)
}
#==================

grid_gen <- function(dimvecs){
   ## Generating a multi dimensional grid based on vectors on each dimension.
   ## The input is a list of vectors of points for each dimension in order and the output is a matrix
   ## where each row is a point and the coordinates are in the columns.
   grd <- expand.grid(dimvecs,KEEP.OUT.ATTRS = F)
   grd <- as.matrix(grd) # The size of the first dimension of this matrix is \bar{m}, based on Wood and Chan paper notations.
   colnames(grd) <- NULL
   grd <- grd + 0 # to make it a numeric matrix, instead of integer.
   return(grd)
}
#==================


# for 2D, generate grid
plot_grid2D = function(ndim){
   nx = ndim[1] - 1
   ny = ndim[2] - 1
   x = c(0:nx)/nx
   y = c(0:ny)/ny
   return(list(x = x, y = y))
}


filter_vec <- function(vec,mdim,ndim){ ## This function is obsolete and not used anymore. I use a faster way in my code now!
   ## Producing \tilde{h}/n for a vector h, based on equation (3.8) and (3.9) in Wood and Chan paper.
   ## The input is vector h with its values at each dimension, number of augmented points and the original number of grid points at each dimension.
   inds <- vec > mdim/2 & vec <= mdim-1
   vec[inds] <- (vec-mdim)[inds]
   return(vec/ndim)
}
#==================

myfft2d <- function(h,inv=F){
   
   ## Taking 2-d FFT of a 2 dimensional matrix h. Note that at this point this is the fastest way to do that, until I have time to write a
   ## wrapper for multidimensional FFT of FFTW3 C library. It returns a vector of all fourier transformation values in the order of writing
   ## each point of I(m) in \bar{m} rows where each column was a dimension. Another reason for returning one dimension is that it is better handled and
   ## storing is also more efficient than two dimensional case.
   h <- as.matrix(h)
   out <- t(mvfft(t(mvfft(h,inverse=inv)),inverse = inv))
   out <- as.vector(out)
   if(inv) out <- out/length(out)
   return(out)
}
#==================

# denominator of the grid: ndim!
eigvals <- function(ndim,gvec,nu,lambda){
   ## This function computes the eigenvalues of the augmented BCCB matrix based on the algorithm
   ## described in the Wood and Chan paper. The inputs are gvector where m[i]=2^gvec[i] for dimension i
   ## and ndim is a vector containing the number of grid points at each dimension.
   
   mdim <- 2^gvec;
   # Augmenting each dimension and finding \tilde{h}/n which is then used as input to the kernel function
   dimvecs <- lapply(mdim,function(x)seq(0,x-1)) # a list containing vectors for each dimension to make grids out of.
   aug_grid <- grid_gen(dimvecs) # gives the matrix containing each of the I(m) grid points in each row. Every j \in I(m) is in row j.
   
   # Filtering the matrix based on equation (3.8) and (3.9) in Wood and Chan paper and producing \tilde{h}/n:
   q <- prod(mdim)
   ndim_aug <- matrix(rep(ndim,q),q,2,byrow=T)
   mdim_aug <- matrix(rep(mdim,q),q,2,byrow=T)
   inds <- (aug_grid > mdim_aug/2) & (aug_grid <= mdim_aug-1)
   aug_grid[inds] <- (aug_grid-mdim_aug)[inds]
   tilde_h_n <- aug_grid/ndim_aug
   #tilde_h_n <- t(apply(aug_grid,1,filter_vec,mdim=mdim,ndim=ndim)) # the slower function that is not used anymore.
   
   tild_aux <- sqrt(rowSums(tilde_h_n^2))
   small_c <- Matern(tild_aux,nu=nu,range=lambda) # The c(j) used in equation (3.12) for all j\in\I(m)
   #small_c <- apply(tilde_h_n,1,thekernel,nu=nu,lambda=lambda) # Obselete part as it is very slow comparing to its equivalent above.
   
   cc <- array(small_c,dim=mdim) # This is the d-dimensional matrix which we take the d-dimensional Fourier transform. There are mdim[l] points at dimension 'l'.
   lambda <- myfft2d(cc) # Now we only support 2-d Fourier transform, (refer to my notes). The output here is a one dimensional vector corresponding to each element of small_c.
   lambda <- Re(lambda) # They are all real, but we want to check if they are all positive.
   return(lambda)
}
#==================

eigvals_exact <- function(ndim,nu,lambda_g){
   ## This function finds the optimized value for the gvector such that m[i] = 2^g[i], in the Wood and Chen paper, for
   ## obtaining the exact value of a positive definite covariance matrix.
   ## Note that our covariance structure is even in all dimesions so m[i] is a power of 2.
   ## It checks eigen values to make sure the matrix is positive definite. It then returns the gvector along
   ## with the vector of eigen values which is then used in sampling step. It uses multidimensional FFT to derive those eigen
   ## values. The input to this function is the number of grid points in each dimension in a 'd' dimensional vector of 'ndim'
   
   # Finding the smallest power of 2 for each dimension that satisfies the condition m/2 >= n-1, as the initial vlaue.
   mdim <- 2^(ceiling(log(ndim + 1, 2)))
   ginit <- log(mdim,2)
   badind <- which(mdim/2 < (ndim-1)); ginit[badind] = ginit[badind]+1
   i0 <- ginit[1]; j0 <- ginit[2]
   # initialize value
   i = i0
   j = j0
   N <- 13 # until 2^13; search m
   gout <- matrix(0,N,N)
   while (i + j <= 2 * N){
      while (j <= j0 + i - i0){
         gout[i,j] <- all(eigvals(ndim,c(i,j),nu,lambda_g)>0)
         if (gout[i, j] == TRUE){break}
         j <- j+1
      }
      if (j == N+1){
         j = j0 ## revert to the initial j
      }
      if (gout[i, j] == TRUE){
         break}
      i <- i+1
   }
   gout2 <- as.vector(gout)
   inds <- which(gout2>0) # save only True indices
   if(!length(inds)){
      stop("No combination of g vector is found to make the matrix positive definite!")
   } else {
      ind <- min(inds)
      rr <- ind%%N
      cc <- ind%/%N + (rr>0)
      gvecout <- c(rr,cc)
      eg <- eigvals(ndim,gvecout,nu,lambda_g)
      mvec <- 2^gvecout
      return(list(egvals=eg,mvec=mvec))
   }
}

calc_pos <- function(j,dimref){
   ## This function finds the position of a given grid point. The numbering is from dimension 1 up to the last point and
   ## the next dimension does not start until all the previous ones are finished. It counts from 1 to mbar. The input `j'
   ## is either the vector containing numbers 0 to dimref[i]-1 at dimension i, or a matrix containing set of such points.
   ## and `dimref' is the number of points at each dimension and a reference to the grid we want to find the number of the given vector in.
   d <- length(dimref)
   cpd <- cumprod(dimref)
   if(!length(dim(j))){
      pos <- j[1] + sum(cpd[-d]*j[-1]) +1
      return(pos)
   } else {
      pos <- j[,1] + j[,-1] %*% as.matrix(cpd[-d]) + 1
      return(pos)
   }
}
#==================

samp_from_grid <- function(ndim,mdim,egs,nu,lambda_g, seed = 0){
   ## Efficient sampling from a regular grid by constructing a BCCB matrix from a BTTB matrix,
   ## based on Wood and Chan paper, section 5.2.4. `ndim' is a vector showing number of
   ## grid points at each dimension. `mdim' is the vector extended grid points in each dimension
   ## and `egs' is the eigen values corresponding to `mdim' obtained by Wood and Chan paper algorithm.
   ## Note that these two inputs are generated by egs_gen() function.
   ## This algorithm currently works only for two dimensions.
   
   d <- length(mdim) # this is the dimension of the problem
   L <- c(1:d)
   mplus <- floor(mdim/2) + 1;
   mbar <- prod(mdim); mplusbar <- prod(mplus)
   a_vec <- numeric(length(egs)) # The final a vector
   tt <- numeric(length(egs))
   
   # Generating p(r) set of grid points, basically the I(m+) set.
   dimvecs <- lapply(mplus,function(x)seq(0,x-1))
   main_dimvecs <- lapply(mdim,function(x)seq(0,x-1))
   main_grid <- grid_gen(main_dimvecs)
   set_pr <- grid_gen(dimvecs)
   
   r <- 0
   set.seed(seed)
   while (r < mplusbar){
      r <- r+1
      j <- set_pr[r,]
      a_j <- which(j < mdim/2 & j>0)
      len <- length(a_j)
      if(len){ # if the cardinality of A(j) is positive.
         a <- c(1:len)
         pset <- powerset(a)
         s <- 0
         while (s < 2^(len-1)){
            s  <- s+1
            b_s <- a_j[pset[[s]]]
            bc_s <- a_j[pset[[2^len-s+1]]]
            j1 <- j; j2 <- j;
            j1[L%in%b_s] <- mdim[L%in%b_s]-j[L%in%b_s]
            j2[L%in%bc_s] <- mdim[L%in%bc_s]-j[L%in%bc_s]
            pos1 <- calc_pos(j1,mdim)
            pos2 <- calc_pos(j2,mdim)
            u <- rnorm(1); v <- rnorm(1)
            a_vec[pos1] <- (2*mbar)^(-0.5)*sqrt(egs[pos1])*(u+v*1i)
            a_vec[pos2] <- (2*mbar)^(-0.5)*sqrt(egs[pos1])*(u-v*1i)
         }
      } else {
         jmat <- matrix(rep(j,mbar),mbar,d,byrow = T) # to find the position of j in the main_grid
         pos <- which(rowSums(main_grid==jmat) > 0)
         u <- rnorm(1)
         a_vec[pos] <- mbar^(-0.5)*sqrt(egs[pos])*u
      }
   }
   a_mat <- array(a_vec,dim=mdim)
   samp1 <- Re(myfft2d(a_mat)) # This is the Fourier transform on all grid points of the extended grid.
   # According to equation (5.2) in the Wood and Chan paper, we only have to choose those that are in
   # in the original grid as k changes in I(n) and not I(m). The following chooses the subset of those
   # sample points that are in the original grid.
   n_dimvecs <- lapply(ndim + 1,function(x)seq(0,x-1)) ## to sample from 0/n to n/n
   n_grid <- grid_gen(n_dimvecs)
   
   npos <- calc_pos(n_grid,mdim) # finding which positions in I(m) is infact I(n).
   outsample <- samp1[npos] # This is infact equation (5.2) in Wood and Chan paper where only those points that are in I(n) are being extracted.
   return(outsample)
}
#==================

h_j = function(x, my_knot, delta_N){
   N = length(my_knot) - 1 # 0 / N to N / N
   i = pmin(1 + floor(x * N), N)
   k = rep(0, N+1) # store the value of h_0(x) to h_N(x)
   k[i] = 1 - abs(x - my_knot[i]) / delta_N
   k[i+1] = 1 - abs(x - my_knot[i+1]) / delta_N
   return(k)
}

# f_N - with tent basis
# input x has two components : (x, y) coordinates - only one point. Not used so far.
# g is N1 by N2 matrix
# output: function values at one points x = (x, y) when the function value at the grid points are given as g

f_N_h_2D = function(x, g){
   N1 = dim(g)[1] - 1
   N2 = dim(g)[2] - 1
   knot_N1 = c(0:N1)/N1
   knot_N2 = c(0:N2)/N2
   i = min(1 + floor(x[1] * N1), N1)
   j = min(1 + floor(x[2] * N2), N2)
   wx1 = (1 - abs(x[1] - knot_N1[i]) * N1)
   wx2 = (1 - abs(x[1] - knot_N1[i+1]) * N1)
   wy1 = (1 - abs(x[2] - knot_N1[j]) * N2)
   wy2 = (1 - abs(x[2] - knot_N1[j+1]) * N2)
   value =  wx1 * wy1 * g[i, j] +
      wx1 * wy2 * g[i, j+1] +
      wx2 * wy1 * g[i+1, j] +
      wx2 * wy2 * g[i+1, j+1]
   return(value)
}
# The same function as f_N_h_2D without that this function works for multiple points x. 
# input x : n by 2 matrix (x1, y1), (x2, y2), ..., (xn, yn)
# output: function values at points x when the function value at the grid points are given as g
# f_N values at x_1, ..., x_n using the values of g.
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

Phi_2D = function(X, N){
   n = nrow(X)
   gridmat = cbind(rep(c(0:N)/N, each = N + 1),
                   rep(c(0:N)/N, N + 1))
   knot_N = c(0:N)/N
   Phi = matrix(0, nrow = n, ncol = (N+1)^2)
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

## log-likelihood for z = g(x) + epsilon given g at grid points
## g : (N1 + 1) by (N2 + 1) matrix
## x : n by 2 matrix (each row are x, y coordinates)
loglik = function(z, x, g, sigsq){
   f_N_x = f_N_h_2D_multi(x, g) # n by 1 vector when x is n by 2 vector
   loglik_value = - sum((z - f_N_x)^2) / (2 * sigsq) - n / 2 *log((2*pi) * sigsq)
   return(loglik_value)
}

# gdouble : 2N + 1 by 2N + 1 matrix
# g : N + 1 by N + 1 matrix
# calculate the log likelihood to check if doubling should be preferred
loglik_double = function(gdouble, g, nu, l, result_N = NULL, result_2N = NULL){
   if (all.equal(dim(gdouble), 2 * dim(g) - 1) != TRUE){
      stop("the dimension of g double should be 2 * dim(g) - 1 !")
   }
   N1 = dim(g)[1] - 1
   N2 = dim(g)[2] - 1
   # embeded vector gdouble_embed = (u, g)
   # original vector gdouble = (u, x)
   gdouble_embed = gdouble
   gdouble_embed[2 * c(1:(N1 + 1)) - 1, 2 * c(1:(N2 + 1)) - 1] = g
   ###################################
   ## align the index order of X_N to be the same as g 
   X_N = cbind(rep(c(0:N1)/N1, N2 + 1), 
               rep(c(0:N2)/N2, each = N1 + 1)) # (N1+1)*(N2+1) by 2 matrix
   X_2N = cbind(rep(c(0:(2*N1))/(2*N1), 2*N2 + 1), 
                rep(c(0:(2*N2))/(2*N2), each = 2*N1 + 1))
   Sigma_N = thekernel(X_N, nu, lambda = l)
   Sigma_2N = thekernel(X_2N, nu, lambda = l)
   ###################################
   ## 1. apply the BTTB inversion algorithm
   if (is.null(result_N) == TRUE){
      result_N = inv_BTTB(Sigma_N, N1 + 1)
      result_2N = inv_BTTB(Sigma_2N, 2*N1 + 1)
   }
   ## 2. calculate det and log likelihood
   result = -0.5 * sum(as.vector(crossprod(as.vector(gdouble_embed), result_2N[[1]])) * gdouble_embed) +
      0.5 * sum(as.vector(crossprod(as.vector(g), result_N[[1]])) * g) +
      0.5 * sum(as.vector(crossprod(as.vector(gdouble), result_2N[[1]])) * gdouble)
   # + 0.5 * result_2N[[2]] - 0.5 * result_N[[2]] - N1 / 2 * log(2*pi)
   return(result)
}

ESS_double = function(gdouble, g, gdouble_ess, nu, l, result_N = NULL, result_2N = NULL){
   ## unify the overlapped points 
   N1 = dim(g)[1] - 1
   N2 = dim(g)[2] - 1
   ## ESS cycle
   thetamin = 0; 
   thetamax = 2*pi;
   u = runif(1)
   logy = loglik_double(gdouble, g, nu, l, result_N, result_2N) + log(u); 
   theta = runif(1,thetamin,thetamax); 
   thetamin = theta - 2*pi; 
   thetamax = theta;
   gdouble_prime = gdouble*cos(theta) + gdouble_ess*sin(theta)
   while(loglik_double(gdouble_prime, g, nu, l, result_N, result_2N) <= logy){
      if(theta < 0)
         thetamin = theta
      else
         thetamax = theta
      theta = runif(1,thetamin,thetamax)
      gdouble_prime = gdouble*cos(theta) + gdouble_ess*sin(theta)
   }
   # replace the nested part with g
   gdouble_prime[2 * c(1:(N1 + 1)) - 1, 2 * c(1:(N2 + 1)) - 1] = g
   return(gdouble_prime)
}

ESS = function(g, nu_ess, z, x, sigsq, Temper = 1, seed = 0){
   thetamin = 0; 
   thetamax = 2*pi;
   set.seed(seed)
   u = runif(1)
   logy = loglik(z, x, g, sigsq) / Temper + log(u); 
   theta = runif(1,thetamin,thetamax); 
   thetamin = theta - 2*pi; 
   thetamax = theta;
   gprime = g*cos(theta) + nu_ess*sin(theta);
   while(loglik(z, x, gprime, sigsq) / Temper <= logy){
      if(theta < 0)
         thetamin = theta
      else
         thetamax = theta
      theta = runif(1,thetamin,thetamax)
      gprime = g*cos(theta) + nu_ess*sin(theta)
   }
   return(gprime)       
}

circulant=function(x){
   n = length(x)
   mat = matrix(0, n, n)
   for (j in 1:n) {
      mat[j, ] <- c(x[-(1:(n+1-j))], x[1:(n+1-j)])
   }
   return(mat)
}

# 1D Levinson algorithm
# gamma_mat: p by p toeplitz 
# y: a p-dimensional vector
# inversion: boolean input. if TRUE, it returns (gamma)^{-1 else:
# output: a solution of (gamma)b = y (when inversion is FALSE)
Levinson = function (gamma_mat, y, inversion = FALSE){
   gamma = gamma_mat[1, ]
   ## Levinson algorithm
   M = length(gamma) - 2 # suppose M = 3
   if (M <= -1){
      stop("gamma should have length of at least 2!")
   }
   if (inversion == FALSE){
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
      return(alpha)
   } else{
      p = M + 2 # nrow of the mat
      inv_mat = matrix(0, p, p)
      for (i in 1:p){
         e_i = rep(0, p)
         e_i[i] = 1
         inv_mat[, i] = Levinson(gamma_mat, e_i, inversion = FALSE)
      }
      return(inv_mat)
   }
}


# input: pd by d matrix
# output: block transposed & reversed order d by pd matrix
BlockTrans_hattilde = function(Mat, long = TRUE){
   Mat_hat = BlockTrans_hat(Mat, long) # pd by d matrix
   if (long == TRUE){
      d = ncol(Mat)
      pd = nrow(Mat)
      p = pd / d
      Mat_out = matrix(0, d, pd)
      for (i in 1:p){
         Mat_out[1:d, (i-1)*d + 1:d] = Mat_hat[(i-1)*d + 1:d, 1:d] # reversing blocks
      }
      return(Mat_out)
   }
}

# input: pd by d matrix / d by pd matrix
# output: pd by d matrix with reversed order / d by pd matrix
BlockTrans_hat = function(Mat, long = TRUE){
   if (long == TRUE){
      d = ncol(Mat)
      pd = nrow(Mat)
      p = pd / d
      Matout = matrix(0, pd, d)
      for (i in 1:p){
         # i is a block index
         Matout[(i-1)*d + 1:d, 1:d] = 
            Mat[(p-i)*d + 1:d, 1:d]
      }
   } else {
      pd = ncol(Mat)
      d = nrow(Mat)
      p = pd / d
      Matout = matrix(0, d, pd)
      for (i in 1:p){
         # i is a block index
         Matout[1:d, (i-1)*d + 1:d] = 
            Mat[1:d, (p-i)*d + 1:d]
      }
   }
   return(Matout)
}


# input: pd by pd size matrix and the block size
# output: the inverse of the large matrix - O(p^2d^2) complexity, and the log of the det
inv_BTTB = function(Mat, blocksize){
   pd = nrow(Mat) # following the notation in Akaike
   if (pd - as.integer(pd) != 0){
      stop("The number of the row of the matrix should the square of the integer!")
   }
   d = blocksize
   if(pd %% d != 0){
      stop("nrow(Mat) should be a multiple of the blocksize!")
   }
   p = pd / d 
   R_list = list()
   for (i in 1:p){
      R_list[[i]] = Mat[1:d, (i-1)*d + 1:d]
   }
   # initial values
   inv_L2 = solve(Mat[1:(2*d), 1:(2*d)])
   q = s = inv_L2[1:d, 1:d]
   j = inv_L2[(d+1):(2*d), 1:d] 
   f = inv_L2[1:d, (d+1):(2*d)]
   M = inv_L2[(d+1):(2*d), (d+1):(2*d)]
   a_tot = r_tot = Mat[(d+1) : pd, 1:d] # (p-1)*d by d matrix
   inv_L_list = list()
   det_Linv_list = list()
   inv_L_list[[1]] = solve(R_list[[1]]) # can be replaced with the Levinson algorithm
   inv_L_list[[2]] = inv_L2
   det_Linv_list[[1]] = log(det(inv_L_list[[1]]))
   det_Linv_list[[2]] = log(det(inv_L_list[[2]]))
   # starting values of qinv, f, sinv, l
   qinv = solve(q)
   qinvf = qinv %*% f
   jsinv = j %*% qinv # since q and s are same
   for (index in 2:(p-1)){
      a_j = a_tot[1:((index-1)*d), 1:d] # d by (j-1)d
      r_j = r_tot[1:((index-1)*d), 1:d] # d by (j-1)d
      ## when qinv (d*d), f(d * (j-1)d), sinv, l is given
      qinvf = cbind(matrix(0, d, d), qinvf) - 
         (R_list[[index + 1]] + qinvf %*% a_j) %*% q %*% cbind(diag(d), BlockTrans_hat(qinvf, long = FALSE)) # matrix size d * (p+1)d
      jsinv = rbind(jsinv, matrix(0, d, d)) - rbind(BlockTrans_hat(jsinv, long = TRUE), diag(d)) %*% 
         q %*% (R_list[[index + 1]] +  BlockTrans_hattilde(a_j)%*% jsinv)
      qinv = qinv - qinvf[1:d, 1:d] %*% qinvf[1:d, 1:d] %*% qinv
      # restoring f and j
      q = solve(qinv)
      f = q %*% qinvf # fat matrix
      j = jsinv %*% q # long matrix
      M = inv_L_list[[index]] + jsinv %*% BlockTrans_hat(f, long = FALSE)
      ## augment elements
      inv_L_list[[index+1]] = rbind(cbind(q, BlockTrans_hat(f, long = FALSE)), cbind(j, M))
      det_Linv_list[[index+1]] = det_Linv_list[[index]] + log(det(q))
   }
   # does not return det since it is not needed for ESS
   return(list(inv_L_list[[p]], det_Linv_list[[p]]))
}

# gridmat is a (gridsize^2) by 2 matrix!
# input: g_N value ((N+1)^2 by 2 matrix) and grid matrix at which we want to obtain f_N value
# ouput: matrices of f_N values - low1, upp1, true, lowsup, uppsup, means
glist_to_plotdf_2D = function(g, gridmat, truefun = f0, alpha1 = 0.95, alpha2 = 0.9){
   # Nsample = number of samples (em)
   Nsample = length(g)
   N = nrow(gridmat)
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