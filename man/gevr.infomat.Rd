% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univdist.R
\name{gevr.infomat}
\alias{gevr.infomat}
\title{Observed information matrix for GEV distribution (return levels)}
\usage{
gevr.infomat(par, dat, method = c("obs", "exp"), p, nobs = length(dat))
}
\arguments{
\item{par}{vector of \code{retlev}, \code{scale} and \code{shape}}

\item{dat}{sample vector}

\item{method}{string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} - the default) information matrix.}

\item{p}{tail probability, corresponding to \eqn{(1-p)}th quantile for \eqn{z}}

\item{nobs}{integer number of observations}
}
\description{
The information matrix is parametrized in terms of return level (\eqn{(1-p)}th quantile), scale and shape.
}
\seealso{
\code{\link{gevr}}
}
\keyword{internal}
