\name{dataframe_DLM}
\alias{dataframe_DLM}
\title{Creating a DLM-Ready Dataframe}
\usage{
dataframe_DLM(X, lag, dynamic_vars = NULL, arglag = list(fun = "bs"), ...)
}
\arguments{
\item{X}{A dataframe, or something that can be coerced into one.}

\item{lag}{Lag length of the dynamic variable. See \code{\link[dlnm:crossbasis]{crossbasis}} for more details.}

\item{dynamic_vars}{The column name or indices that correspond to the longitudinal variables. If left missing, the dataset is not altered in any way and a warning is supplied.}

\item{arglag}{A list that is passed into \code{\link[dlnm:onebasis]{onebasis}} for generating a basis matrix for the lag space. See \code{\link[dlnm:crossbasis]{crossbasis}}.}

\item{...}{Further arguments to be passed into the \code{\link[dlnm:crossbasis]{crossbasis}} function.}
}
\description{
Transforms a dataframe to prepare it for distributed lag modelling.
}
\details{
The purpose of this function is to streamline the preperation of the data for distributed lag modelling in the case where every dynamic variable is handled the same way.

If no dynamic variables are given, the input dataframe is returned unaltered with a warning. Otherwise, the function returns a \emph{dataframe_DLM} object, which is essentially treated like a standard dataframe but contains extra information to simplify the process of distributed lag modelling when using the \emph{NB_MCMC} or \emph{QB_MCMC} functions.
}
\value{
The input is returned unchanged if no dynamic variables are given. Otherwise a \emph{dataframe_DLM} object is returned. See details.
}
\examples{
X <- dplyr::select( dlnm::chicagoNMMAPS, c('cvd', 'dow', 'temp', 'dptp', 'o3') )
X <- na.omit( X )
arglag <- list( fun = 'bs', df = 4 )
DLM_dat <- dataframe_DLM( X, lag = 40, dynamic_vars =  c('temp', 'dptp', 'o3'), arglag = arglag )
}
\author{
Daniel Dempsey (\href{mailto:daniel.dempsey0@gmail.com}{daniel.dempsey0@gmail.com})
}
