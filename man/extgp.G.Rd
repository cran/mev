% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/egp.R
\name{extgp.G}
\alias{extgp.G}
\title{Carrier distribution for the extended GP distributions of Naveau et al.}
\arguments{
\item{u}{vector of observations (\code{dextgp.G}), probabilities (\code{qextgp.G}) or quantiles (\code{pextgp.G}), in \eqn{[0,1]}}

\item{prob}{mixture probability for model \code{type} \code{4}}

\item{kappa}{shape parameter for \code{type} \code{1}, \code{3} and \code{4}}

\item{delta}{additional parameter for \code{type} \code{2}, \code{3} and \code{4}}

\item{type}{integer between 0 to 5 giving the model choice}

\item{log}{logical; should the log-density be returned (default to \code{FALSE})?}

\item{n}{sample size}

\item{unifsamp}{sample of uniform; if provided, the data will be used in place of new uniform random variates}

\item{censoring}{numeric vector of length 2 containing the lower and upper bound for censoring}

\item{direct}{logical; which method to use for sampling in model of \code{type} \code{4}?}
}
\description{
Density, distribution function, quantile function and random number
generation for the carrier distributions of the extended Generalized Pareto distributions.
}
\section{Usage}{
 \code{pextgp.G(u, type=1, prob, kappa, delta)}

 \code{dextgp.G(u, type=1, prob=NA, kappa=NA, delta=NA, log=FALSE)}

 \code{qextgp.G(u, type=1, prob=NA, kappa=NA, delta=NA)}

 \code{rextgp.G(n, prob=NA, kappa=NA, delta=NA,
type=1, unifsamp=NULL, direct=FALSE, censoring=c(0,1))}
}

\seealso{
\code{\link{extgp}}
}
\author{
Raphael Huser and Philippe Naveau
}
