sample_fullGP <- function(X, Z, sigma, nu, 
                          kappa_pr, kappa_sampler, 
                          tausq_pr, tausq_sampler, grid_size = 40, target) {
   # Covariance matrix function for Matern kernel
   CovMat_XX <- function(X1, X2, nu, lambda, tausq) {
      if (ncol(X1) != ncol(X2)) stop("X1, X2 dimensions not aligned")
      X <- rdist(X1, X2)
      out <- Matern(X, nu = nu, range = lambda) * tausq
      return(out)
   }
   # Create grid for predictions
   gridmat <- cbind(
      rep(seq(0, 1, length.out = grid_size + 1), each = grid_size + 1),
      rep(seq(0, 1, length.out = grid_size + 1), grid_size + 1)
   )
   N_list = kappa_list = tausq_list = vector(length = target)
   PostSamples = matrix(0, target, nrow(gridmat))
   # initial kappa
   kappa = kappa_sampler()
   tausq = tausq_sampler()
   K_XX <- CovMat_XX(X, X, nu, lambda = 1 / kappa, tausq = tausq)
   K_XnewXnew <- CovMat_XX(gridmat, gridmat, nu, lambda = 1 / kappa, tausq = tausq)
   K_XnewX <- CovMat_XX(gridmat, X, nu, lambda = 1 / kappa, tausq = tausq)
   inv_K_XX_noise <- solve(K_XX + sigma^2 * diag(nrow(K_XX)))
   chol_noise = chol(inv_K_XX_noise)
   log_prob_now = log(kappa_pr(kappa)) + log(tausq_pr(tausq)) - sum(log(diag(chol_noise))) -
      1/2 * Z %*% crossprod(inv_K_XX_noise, Z) 
   Xnew_mean <- (K_XnewX %*% inv_K_XX_noise %*% Z)[,]
   Xnew_var <- K_XnewXnew - K_XnewX %*% inv_K_XX_noise %*% t(K_XnewX)
   # PostSamples <- matrix(0, nrow = target, ncol = length(gridmat))
   for (i in 1:target){
      # kappa selection using MH algorithm
      kappa_cand = kappa_sampler()
      tausq_cand = tausq_sampler()
      # Gibbs sampler for kappa
      K_XX_cand <- CovMat_XX(X, X, nu, lambda = 1 / kappa_cand, tausq = tausq_cand)
      inv_K_XX_noise_cand <- solve(K_XX_cand + sigma^2 * diag(nrow(K_XX_cand)))
      chol_noise_cand = chol(inv_K_XX_noise_cand)
      log_prob_cand = log(kappa_pr(kappa_cand)) + log(tausq_pr(tausq_cand)) - sum(log(diag(chol_noise_cand))) - 
         1/2 * Z %*% crossprod(inv_K_XX_noise_cand, Z)
      u = runif(1)
      if (u < exp(log_prob_cand - log_prob_now)){
         kappa = kappa_cand
         tausq = tausq_cand
         log_prob_now = log_prob_cand
         K_XX <- K_XX_cand
         K_XnewXnew <- CovMat_XX(gridmat, gridmat, nu, lambda = 1 / kappa, tausq = tausq)
         K_XnewX <- CovMat_XX(gridmat, X, nu, lambda = 1 / kappa, tausq = tausq)
         Xnew_mean <- (K_XnewX %*% inv_K_XX_noise %*% Z)[,]
         inv_K_XX_noise <- inv_K_XX_noise_cand
         Xnew_var <- K_XnewXnew - K_XnewX %*% inv_K_XX_noise %*% t(K_XnewX)
         print(kappa)
         print(tausq)
      }
      PostSamples[i, ] = rmvnorm(1, mean = Xnew_mean, sigma = Xnew_var)[,] + rnorm(n = length(Xnew_mean), sd = sigma)
      kappa_list[i] = kappa
      tausq_list[i] = tausq
      print(i)
      print(kappa)
      print(tausq)
   }
   PostMean <- colMeans(PostSamples)
   # Return posterior mean 
   return(list(PostMean = PostMean, kappa_list = kappa_list, tausq_list = tausq_list))
}
