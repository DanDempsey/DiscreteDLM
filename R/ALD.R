pALD <- function(q, mu = 0, sigma = 1, p = 0.5) {
  ifelse( test = q < mu,
          yes = p * exp((1 - p) * (q - mu)/sigma),
          no = 1 - (1 - p) * exp(-p * (q - mu)/sigma) )
}

qALD <- function(prob, mu = 0, sigma = 1, p = 0.5)  {
  ifelse( test = prob < p,
          yes = mu + (sigma*log(prob/p))/(1 - p),
          no = mu - sigma*log((1 - prob)/(1 - p))/p )
}

rALD <- function(n, mu = 0, sigma = 1, p = 0.5) {
  u <- runif(n)
  mapply(qALD, prob = u, mu = mu, sigma = sigma, p = p)
}

dALD <- function(x, mu = 0, sigma = 1, p = 0.5) {
  ifelse(test = x < mu,
         yes = (p * (1 - p)/sigma) * exp((1 - p) * (x - mu)/sigma),
         no = (p * (1 - p)/sigma) * exp(-p * (x - mu)/sigma))
}

rTALD <- function(n, upper_tail, mu = 0, sigma = 1, p = 0.5) {
  bound <- pALD( 0, mu = mu, sigma = sigma, p = p )
  u <- runif( n, min = ifelse(upper_tail, bound, 0), max = ifelse(upper_tail, 1, bound) )
  qALD( prob = u, mu = mu, sigma = sigma, p = p )
}

