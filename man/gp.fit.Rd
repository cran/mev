% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gp.R
\name{gp.fit}
\alias{gp.fit}
\title{Maximum likelihood estimate of generalized Pareto applied to threshold exceedances}
\usage{
gp.fit(
  xdat,
  threshold,
  method = c("Grimshaw", "auglag", "nlm", "optim", "ismev", "zs", "zhang"),
  show = FALSE,
  MCMC = NULL,
  fpar = NULL,
  warnSE = TRUE
)
}
\arguments{
\item{xdat}{a numeric vector of data to be fitted.}

\item{threshold}{the chosen threshold.}

\item{method}{the method to be used. See \bold{Details}. Can be abbreviated.}

\item{show}{logical; if \code{TRUE} (the default), print details of the fit.}

\item{MCMC}{\code{NULL} for frequentist estimates, otherwise a boolean or a list with parameters passed. If \code{TRUE}, runs a Metropolis-Hastings sampler to get posterior mean estimates. Can be used to pass arguments \code{niter}, \code{burnin} and \code{thin} to the sampler as a list.}

\item{fpar}{a named list with fixed parameters, either \code{scale} or \code{shape}}

\item{warnSE}{logical; if \code{TRUE}, a warning is printed if the standard errors cannot be returned from the observed information matrix when the shape is less than -0.5.}
}
\description{
The function \code{\link[mev]{fit.gpd}} is a wrapper around \code{gp.fit}
}
\keyword{internal}
