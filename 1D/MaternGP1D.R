get_r <- function(x1, x2) {
   outer(x1, x2, FUN = "-")
}

matern_kernel <- function(r, l = 1, v = 1) {
   r <- abs(r)
   r[r == 0] <- 1e-8
   part1 = 2^(1-v) / gamma(v)
   part2 = (sqrt(2*v) * r / l)^v
   part3 = besselK(sqrt(2*v) * r / l, nu = v)
   return(part1 * part2 * part3)
}


sample_posterior <- function(
      x, y, x_star, l, v,
      sigma_n = 0.1, random_seed = -1, n_samples = 5
) {
   k_xx = matern_kernel(get_r(x, x), l, v)
   k_xxs = matern_kernel(get_r(x, x_star), l, v)
   k_xsx = matern_kernel(get_r(x_star, x), l, v)
   k_xsxs = matern_kernel(get_r(x_star, x_star), l, v)
   
   I = diag(1, dim(k_xx)[1])
   k_xx_noise = solve(k_xx + sigma_n ^ 2 * I)
   kxsx_kxxNoise = k_xsx %*% k_xx_noise
   # Eq.2.23, 24
   fsb = kxsx_kxxNoise %*% y
   cov_fs = k_xsxs - kxsx_kxxNoise %*% k_xxs
   random_seed <- as.integer(random_seed)
   if (random_seed > 0) set.seed(random_seed)
   return(MASS::mvrnorm(n_samples, fsb, cov_fs))
}

# Example
# grid = c(1:50)/50
# l = 1
# v = 1
# tmp = sample_posterior(X, Z, grid, l, v, sigma_n = 0.1, n_samples = 20)
