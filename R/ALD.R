#' The Asymmetric Laplacian Distribution
#'
#' Density, distribution function, quantile function and random generation for the Asymmetric Laplacian Distribution (ALD). Also contains random number generation for a truncated ALD. Our code here is heavily derived from the "ald" package by Galarza and Lachos (2021), adapted to vectorize the mu parameter.
#'
#' @usage
#' dALD(x, mu = 0, sigma = 1, p = 0.5)
#' pALD(q, mu = 0, sigma = 1, p = 0.5)
#' qALD(prob, mu = 0, sigma = 1, p = 0.5)
#' rALD(n, mu = 0, sigma = 1, p = 0.5)
#' rTALD(n, rtrunc, mu = 0, sigma = 1, p = 0.5)
#' @param x,q vector of quantiles.
#' @param prob vector of probabilities.
#' @param n number of observations.
#' @param mu vector of location parameters.
#' @param sigma vector of scale parameters.
#' @param p ALD skew parameter.
#' @param rtrunc Logical denoting what portion of the distribution is discarded. If TRUE, only positive values are sampled, and negative if FALSE.
#' @details These functions are based on the three parameter ALD:
#' \deqn{f(x|\mu, \sigma, p) = \frac{p(1-p)}{\sigma}exp\left(-\rho_p(\frac{x-\mu}{\sigma})\right)}
#' where
#' \deqn{\rho_p(z) = z(p - I_{z<0})}
#' These functions differ than the ones provided in the "ald" package by vectorising the mu parameter. There's also the addition of rTALD which generates random numbers from the truncated ALD, which is useful for Bayesian quantile regression. Note that truncation is fixed at zero; the user only decides whether the negative or positive axis is discarded.
#' @examples
#' set.seed( 100 )
#'
#' ### Vectorised input
#' random_ALD <- rALD( 1000, mu = runif(1000, -100, 100) )
#' plot( random_ALD )
#'
#' ### Truncated version
#' trunc_random_ALD <- rTALD( 1000, TRUE, mu = 2, sigma = 3, p = 0.75 )
#' plot( dALD(sort(trunc_random_ALD), mu = 2, sigma = 3, p = 0.75) )
#'
#' @author Original code written by Christian E. Galarza and Victor H. Lachos. Edits made by Daniel Dempsey.
#' @aliases pALD qALD rALD dALD rTALD
#' @import stats
#' @export
pALD <- function(q, mu = 0, sigma = 1, p = 0.5) {
  ifelse( test = q < mu,
          yes = p * exp((1 - p) * (q - mu)/sigma),
          no = 1 - (1 - p) * exp(-p * (q - mu)/sigma) )
}

#' @rdname pALD
#' @export
qALD <- function(prob, mu = 0, sigma = 1, p = 0.5)  {
  ifelse( test = prob < p,
          yes = mu + (sigma*log(prob/p))/(1 - p),
          no = mu - sigma*log((1 - prob)/(1 - p))/p )
}

#' @rdname pALD
#' @export
rALD <- function(n, mu = 0, sigma = 1, p = 0.5) {
  u <- runif(n)
  mapply(qALD, prob = u, mu = mu, sigma = sigma, p = p)
}

#' @rdname pALD
#' @export
dALD <- function(x, mu = 0, sigma = 1, p = 0.5) {
  ifelse(test = x < mu,
         yes = (p * (1 - p)/sigma) * exp((1 - p) * (x - mu)/sigma),
         no = (p * (1 - p)/sigma) * exp(-p * (x - mu)/sigma))
}

#' @rdname pALD
#' @export
rTALD <- function(n, rtrunc, mu = 0, sigma = 1, p = 0.5) {
  bound <- pALD( 0, mu = mu, sigma = sigma, p = p )
  u <- runif( n, min = ifelse( rtrunc, bound, 0 ), max = ifelse( rtrunc, 1, bound ) )
  qALD( prob = u, mu = mu, sigma = sigma, p = p )
}

