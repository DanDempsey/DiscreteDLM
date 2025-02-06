\name{ALD}
\alias{pALD}
\alias{qALD}
\alias{rALD}
\alias{dALD}
\alias{rTALD}
\title{The Asymmetric Laplacian Distribution}
\usage{
dALD(x, mu = 0, sigma = 1, p = 0.5)
pALD(q, mu = 0, sigma = 1, p = 0.5)
qALD(prob, mu = 0, sigma = 1, p = 0.5)
rALD(n, mu = 0, sigma = 1, p = 0.5)
rTALD(n, upper_tail, mu = 0, sigma = 1, p = 0.5)
}
\arguments{
\item{mu}{vector of location parameters.}

\item{sigma}{vector of scale parameters.}

\item{p}{ALD skew parameter.}

\item{prob}{vector of probabilities.}

\item{n}{number of observations.}

\item{x, q}{vector of quantiles.}

\item{upper_tail}{Logical denoting what portion of the distribution is retained. If TRUE, only positive values are sampled, and negative if FALSE.}
}
\description{
Density, distribution function, quantile function and random generation for the Asymmetric Laplacian Distribution (ALD). Also contains random number generation for a truncated ALD. Our code here is heavily derived from the "ald" package by Galarza and Lachos (2021), adapted to vectorize the mu parameter.
}
\details{
These functions are based on the three parameter ALD:
\deqn{f(x|\mu, \sigma, p) = \frac{p(1-p)}{\sigma}exp\left(-\rho_p(\frac{x-\mu}{\sigma})\right)}
where
\deqn{\rho_p(z) = z(p - I_{z<0})}
These functions differ than the ones provided in the \emph{ald} package by vectorising the mu parameter. In addition we have implemented a function to sample from the truncated ALD (rTALD), which is needed for Bayesian quantile regression. Note that truncation is fixed at zero; the user only decides whether the negative or positive axis is discarded.
}
\value{
Each function returns a vector. dALD, pALD and qALD returns the PDF, CDF and inverse CDF of the ALD respectively. rALD returns a random sample from the ALD distribution, and rTALD returns a random sample from the ALD truncated at 0.
}
\examples{
set.seed( 100 )

### Vectorised input
random_ALD <- rALD( 1000, mu = runif(1000, -100, 100) )
plot( random_ALD )

### Truncated version
trunc_random_ALD <- rTALD( 1000, TRUE, mu = 2, sigma = 3, p = 0.75 )
plot( dALD(sort(trunc_random_ALD), mu = 2, sigma = 3, p = 0.75) )

}
\author{
Daniel Dempsey (\href{mailto:daniel.dempsey0@gmail.com}{daniel.dempsey0@gmail.com}). Original code from the \emph{ald} package were written by Christian E. Galarza and Victor H. Lachos.
}
