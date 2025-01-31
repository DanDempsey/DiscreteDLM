\name{DiscreteDLM-overview}
\alias{DiscreteDLM-overview}
\alias{DiscreteDLM}
\docType{package}
\title{
Fitting Discrete Distributed Lag Models with Variable Selection via MCMC
}
\description{
DiscreteDLM contains functionality for fitting and visualising Bayesian Distributed Lag Models (DLMs) when the response variable is either count or binary. More specifically, this package contains an implementation of Bayesian quantile binary regression (for binary response data) and negative-binomial regression (for count response data). Additionally, these functions implement variable uncertainty inference for each predictor. Multinomial regression is not implemented.

While these functions can be applied to any standard regression problem, the package was set up with DLMs in mind; the functions therein allow for an easy-to-use pipeline to go from defining the DLM, fitting the DLM, and visualising the results.
}
\details{
\tabular{ll}{
Package: \tab DiscreteDLM\cr
Type: \tab Package\cr
Version: \tab 0.9.8\cr
Date: \tab 2024-11-06\cr
License: \tab GPL (>=3)\cr
}

}

\author{
Daniel Dempsey <\email{daniel.dempsey0@gmail.com}>, with contributions and supervision by Jason Wyse <\email{wyseja@tcd.ie}>
}

\references{
Daniel Dempsey and Jason Wyse. 2025. "Bayesian Variable Selection in Distributed Lag Models: A Focus on Binary Quantile and Count Data Regressions." https://doi.org/10.48550/arXiv.2403.03646.
}

\keyword{package}

\seealso{
\code{\link{dataframe_DLM}},\code{\link{NB_MCMC}},\code{\link{QB_MCMC}}
}

