\name{NB_MCMC}
\alias{NB_MCMC}
\title{Negative Binomial Regression via MCMC}
\usage{
NB_MCMC(
  formula,
  data = NULL,
  nsamp = 1000,
  nburn = 1000,
  thin = 1,
  standardize = TRUE,
  prior_beta_mu = 0,
  prior_beta_sigma = 100,
  prior_gamma_p = 0.5,
  prior_xi_shape = 2,
  prior_xi_scale = 1/50,
  init_beta = 0,
  init_gamma = FALSE,
  init_xi = 1
)
}
\arguments{
\item{formula}{Formula object to set the symbolic description of the model to be fitted.}

\item{data}{Optional dataframe. Ideally this should be a \code{\link[DiscreteDLM:dataframe_DLM]{dataframe_DLM}} object if doing distributed lag modelling.}

\item{nsamp}{The desired sample size from the posterior. Set to 5000 by default.}

\item{nburn}{The number of iterations of the MCMC to be discarded as burn-in. Set to 5000 by default.}

\item{thin}{Thinning factor of the MCMC chain after burn-in. Set to 1 by default.}

\item{standardize}{Logical indicating if the data should be standardised prior to fitting the model to zero mean and unit variance. If there is no intercept, then only scaling is applied (with a warning). The posterior sample of beta is transformed back to the original scale upon completion of the algorithm.}

\item{prior_beta_mu}{Mean of the Gaussian prior for the regression coefficients, beta. Either a vector of length equal to the number of predictors or a single numeric to represent a constant vector.}

\item{prior_beta_sigma}{Covariance matrix of the Gaussian prior for the regression coefficients, beta. Either a square matrix of dimension equal to the number of predictors or a single numeric to represent an isotropic covariance matrix.}

\item{prior_gamma_p}{Probability parameter of the Bernoulli prior for the predictor inclusion parameter, gamma. Either a vector of length equal to the number of predictors or a single numeric, to represent that all predictors have the same prior probability of inclusion.}

\item{prior_xi_shape}{Shape parameter of the Gamma distribution prior for the negative binomial stopping parameter, xi.}

\item{prior_xi_scale}{Scale parameter of the Gamma distribution prior for the negative binomial stopping parameter, xi.}

\item{init_beta}{Initial MCMC values for the beta parameters. Either a vector of length equal to the number of predictors or a single numeric representing the same starting value for each beta component.}

\item{init_gamma}{Initial MCMC values for the gamma parameters. Either a vector of length equal to the number of predictors or a single numeric representing the same starting value for each gamma component.}

\item{init_xi}{Initial MCMC value for xi.}
}
\value{
An "MCMC_DLM" object; a list containing the MCMC-derived posterior sample, as well as the data that were used.
}
\description{
Fits a negative binomial regression model using Gibbs sampling and performs Bayesian variable selection.
}
\details{
This function fits a negative binomial regression model via MCMC. Latent variable representation allows for Gibbs sampling of the parameters; see Pillow and Scott (2012) and Zhou et. al. (2012). The algorithm also includes predictor inclusion uncertainty inference (inferred via a Metroplis step) adapted from Holmes and Held (2006).

The parameters of interest for this model are the regression slopes (beta), the binary predictor inclusion indicator (gamma) and the negative binomial stopping parameter (xi). For further details, see Dempsey and Wyse (2024).

Note that when setting initial values and priors for gamma, the intercept must be included and therefore its initial value and prior are forced to be equal to 1. If an intercept is not used, then the first variable is forced to be included instead.
}
\seealso{
\link{QB_MCMC}, \link{dataframe_DLM}, \link{plot.MCMC_DLM}
}
\examples{
### Set up data
X <- dplyr::select( dlnm::chicagoNMMAPS, c('cvd', 'dow', 'temp', 'dptp', 'o3') )
X <- na.omit( X )
arglag <- list( fun = 'bs', df = 4 )
DLM_dat <- dataframe_DLM( X, lag = 40, dynamic_vars =  c('temp', 'dptp', 'o3'), arglag = arglag )

### Fit model
# NOTE: Only using 100 total iterations for illustration purposes!
myfit <- NB_MCMC( cvd ~ ., data = DLM_dat, nsamp = 50, nburn = 50 )
summary( myfit )
}
\references{
Chris C. Holmes and Leonhard Held. 2006. "Bayesian auxiliary variable models for binary and multinomial regression." Bayesian Analysis 1(1): 145-168.

Daniel Dempsey and Jason Wyse. 2025. "Bayesian Variable Selection in Distributed Lag Models: A Focus on Binary Quantile and Count Data Regressions." https://doi.org/10.48550/arXiv.2403.03646

Jonathan Pillow and James Scott. 2012. "Fully Bayesian inference for neural models with negative-binomial spiking." Advances in neural information processing systems 25.

Mingyuan Zhou, Lingbo Li, David Dunson and Lawrence Carin. 2012. "Lognormal and gamma mixed negative binomial regression." Proceedings of the... International Conference on Machine Learning. International Conference on Machine Learning. Vol. 2012. NIH Public Access.
}
\author{
Daniel Dempsey (\href{mailto:daniel.dempsey0@gmail.com}{daniel.dempsey0@gmail.com})
}

