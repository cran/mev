% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rparpcs.R
\name{rparpcs}
\alias{rparpcs}
\title{Simulation from Pareto processes (max) using composition sampling}
\usage{
rparpcs(n, Lambda = NULL, Sigma = NULL, df = NULL, model = c("br",
  "xstud"), riskf = c("max", "min"), shape = 1)
}
\arguments{
\item{n}{sample size.}

\item{Lambda}{parameter matrix for the Brown--Resnick model. See \bold{Details}.}

\item{Sigma}{correlation matrix if \code{model = 'xstud'}, otherwise
the covariance matrix formed from the stationary Brown-Resnick process.}

\item{df}{degrees of freedom for extremal Student process.}

\item{model}{string indicating the model family.}

\item{riskf}{string indicating the risk functional. Only \code{max} and \code{min} are currently supported.}

\item{shape}{tail index of the Pareto variates (reciprocal shape parameter). Must be strictly positive.}
}
\value{
an \code{n} by \code{d} matrix of samples, where \code{d = ncol(Sigma)}, with \code{attributes} \code{mixt.weights}.
}
\description{
The algorithm performs forward sampling by simulating first from a
mixture, then sample angles conditional on them being less than one.
The resulting sample from the angular distribution is then multiplied by
Pareto variates with tail index \code{shape}.
}
\details{
Only extreme value models based on elliptical processes are handled. The \code{Lambda} matrix
is formed by evaluating the semivariogram \eqn{\gamma} at sites \eqn{s_i, s_j}, meaning that
 \eqn{\Lambda_{i,j} = \gamma(s_i, s_j)/2}.

The argument \code{Sigma} is ignored for the Brown-Resnick model
if \code{Lambda} is provided by the user.
}
\examples{
\dontrun{
#Brown-Resnick, Wadsworth and Tawn (2014) parametrization
D <- 20L
coord <- cbind(runif(D), runif(D))
semivario <- function(d, alpha = 1.5, lambda = 1){0.5 * (d/lambda)^alpha}
Lambda <- semivario(as.matrix(dist(coord))) / 2
rparpcs(n = 10, Lambda = Lambda, model = 'br', shape = 0.1)
#Extremal Student
Sigma <- stats::rWishart(n = 1, df = 20, Sigma = diag(10))[,,1]
rparpcs(n = 10, Sigma = cov2cor(Sigma), df = 3, model = 'xstud')
}
}
\seealso{
\code{\link{rparp}} for general simulation of Pareto processes based on an accept-reject algorithm.
}
\author{
Leo Belzile
}