\name{plot.MCMC_DLM}
\alias{plot.MCMC_DLM}
\title{Visualises Results from MCMC_DLM Fits}
\usage{
\method{plot}{MCMC_DLM}(x, type = "beta", include_intercept = FALSE, print_output = TRUE, ...)
}
\arguments{
\item{x}{An MCMC_DLM object.}

\item{type}{One of "beta", "gamma", "xi", or "lags". Partial matching is supported. See Details below.}

\item{include_intercept}{Logical indicating if the intercept should be included in the graphs.}

\item{print_output}{Logical indicating whether or not the ggplot objects should be printed to the screen.}

\item{...}{Extra arguments; currently unused.}
}
\value{
A ggplot object, except when type = 'lags', in which case it will be a list of ggplot objects. Will display the ggplot visual if print_output = TRUE.
}
\description{
Provides some ggplot visualisations of the posterior from MCMC_DLM objects.
}
\details{
The type of plot generated is determined by the \emph{type} argument. The options are:
\itemize{
\item \emph{beta}: A ridgeplot of the regression slopes,
\item \emph{gamma}: A bar chart of the mean probability of inclusion,
\item \emph{xi}: A kernel density estimate plot of the negative binomial stopping parameter (of course, this only works for negative binomial regression fits),
\item \emph{lags}: A collection of plots (supplied in a list), each containing a line chart of the lag-response for each dynamic variable.
}
}
\author{
Daniel Dempsey (\href{mailto:daniel.dempsey0@gmail.com}{daniel.dempsey0@gmail.com})
}
