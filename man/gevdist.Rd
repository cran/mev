% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{gevdist}
\alias{gevdist}
\alias{qgev}
\alias{rgev}
\alias{dgev}
\alias{pgev}
\title{Generalized extreme value distribution}
\usage{
qgev(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)

rgev(n, loc = 0, scale = 1, shape = 0)

dgev(x, loc = 0, scale = 1, shape = 0, log = FALSE)

pgev(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
\item{p}{vector of probabilities}

\item{loc}{scalar or vector of location parameters whose length matches that of the input}

\item{scale}{scalar or vector of positive scale parameters whose length matches that of the input}

\item{shape}{scalar shape parameter}

\item{lower.tail}{logical; if \code{TRUE} (default), returns the distribution function, otherwise the survival function}

\item{n}{scalar number of observations}

\item{x, q}{vector of quantiles}

\item{log, log.p}{logical; if \code{TRUE}, probabilities \eqn{p} are given as
\eqn{\log(p)}.}
}
\description{
Density function, distribution function, quantile function and
random number generation for the generalized extreme value
distribution.
}
\details{
The distribution function of a GEV distribution with parameters
 \code{loc} = \eqn{\mu}, \code{scale} = \eqn{\sigma} and
 \code{shape} = \eqn{\xi} is
 \deqn{F(x) = \exp\{-[1 + \xi (x - \mu) / \sigma] ^ {-1/\xi} \}}{%
       F(x) = exp{ -[1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)} }
 for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  If \eqn{\xi = 0} the
 distribution function is defined as the limit as \eqn{\xi} tends to zero.

 The quantile function, when evaluated at zero or one,
 returns the lower and upper endpoint, whether the latter is finite or not.
}
\references{
Jenkinson, A. F. (1955) The frequency distribution of the
  annual maximum (or minimum) of meteorological elements.
  \emph{Quart. J. R. Met. Soc.}, \strong{81}, 158-171.
  Chapter 3: \doi{10.1002/qj.49708134804}

Coles, S. G. (2001) \emph{An Introduction to Statistical
  Modeling of Extreme Values}, Springer-Verlag, London.
  \doi{10.1007/978-1-4471-3675-0_3}
}
\author{
Leo Belzile, with code adapted from Paul Northrop
}
