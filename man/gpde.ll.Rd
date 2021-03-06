% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univdist.R
\name{gpde.ll}
\alias{gpde.ll}
\alias{gpde.ll.optim}
\title{Negative log likelihood of the generalized Pareto distribution (expected shortfall)}
\usage{
gpde.ll(par, dat, m)

gpde.ll.optim(par, dat, m)
}
\arguments{
\item{par}{vector of length 2 containing \eqn{e_m} and \eqn{\xi}, respectively the expected shortfall at probability 1/(1-\eqn{\alpha}) and the shape parameter.}

\item{dat}{sample vector}

\item{m}{number of observations of interest for return levels. See \strong{Details}}
}
\description{
Negative log likelihood of the generalized Pareto distribution (expected shortfall)

Negative log likelihood of the generalized Pareto distribution (expected shortfall) - optimization
The negative log likelihood is parametrized in terms of log expected shortfall and shape in order to perform unconstrained optimization
}
\seealso{
\code{\link{gpde}}
}
\keyword{internal}
