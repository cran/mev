% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univdist.R
\name{gpdr.infomat}
\alias{gpdr.infomat}
\title{Observed information matrix for GP distribution (return levels)}
\usage{
gpdr.infomat(par, dat, m, method = c("obs", "exp"), nobs = length(dat))
}
\arguments{
\item{par}{vector of length 2 containing \eqn{y_m} and \eqn{\xi}, respectively the \eqn{m}-year return level and the shape parameter.}

\item{dat}{sample vector}

\item{m}{number of observations of interest for return levels. See \strong{Details}}

\item{method}{string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} - the default) information matrix.}

\item{nobs}{number of observations}
}
\description{
The information matrix is parametrized in terms of rate of \code{m}-year return level and shape
}
\seealso{
\code{\link{gpdr}}
}
\keyword{internal}
