% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Wadsworth_Technometrics_2016.R
\name{W.diag}
\alias{W.diag}
\title{Wadsworth's univariate and bivariate exponential threshold diagnostic}
\usage{
W.diag(xdat, model = c("nhpp", "exp", "invexp"), u = NULL, k, q1 = 0,
  q2 = 1, par = NULL, M = NULL, nbs = 1000, alpha = 0.05,
  plots = c("LRT", "WN", "PS"), UseQuantiles = TRUE, pmar = c(5.5, 7, 3,
  3), tikz = FALSE, ...)
}
\arguments{
\item{xdat}{A numeric vector of data to be fitted.}

\item{model}{String specifying whether the univariate or bivariate diagnostic should be used. Either \code{nhpp}
for the univariate model, \code{exp} (\code{invexp}) for the bivariate exponential model with rate (inverse rate) parametrization. See details.}

\item{u}{Optional; vector of candidate thresholds.}

\item{k}{Number of thresholds to consider (if \code{u} unspecified).}

\item{q1}{Lowest quantile for the threshold sequence.}

\item{q2}{Upper quantile limit for the threshold sequence (\code{q2} itself is not used as a threshold,
but rather the uppermost threshold will be at the \eqn{(q_2-1/k)}{q2-1/k} quantile).}

\item{par}{Parameters of the NHPP likelihood. If \code{missing}, the \code{\link[evd]{fpot}} routine will be run to obtain values}

\item{M}{Number of superpositions or "blocks" / "years" the process corresponds to (can affect the optimization)}

\item{nbs}{Number of simulations used to assess the null distribution of the LRT, and produce the p-value}

\item{alpha}{Significance level of the LRT}

\item{plots}{Which plots to produce; \code{LRT}= likelihood ratio test, \code{WN} = white noise, \code{PS} = parameter stability}

\item{UseQuantiles}{Logical; use quantiles as the thresholds in the plot?}

\item{pmar}{Vector of length 4 giving the arguments for the plot margins in \code{par(mar=c(*,*,*,*))}.}

\item{tikz}{Logical; if \code{TRUE}, axis labels are replaced with \code{LaTeX} code}

\item{...}{Additional parameters passed to \code{plot}.}
}
\value{
plots of the requested diagnostics and a list with components
\itemize{
\item MLE  Maximum likelihood estimates from all thresholds
\item Cov  Joint asymptotic covariance matrix for \eqn{\xi}{xi}, \eqn{\eta}{eta} or \eqn{\eta^{-1}}{1/eta}.
\item WN  Values of the white noise process.
\item LRT Values of the likelihood ratio test statistic vs threshold.
\item pval P-value of the likelihood ratio test.
\item k  Final number of thresholds used.
\item thresh Threshold selected by the likelihood ratio procedure.
\item mle.u Maximum likelihood estimates from selected threshold.
}
}
\description{
Function to produce diagnostic plots and test statistics for the
threshold diagnostics exploiting structure of maximum likelihood estimators
based on the non-homogeneous Poisson process likelihood
}
\details{
The function is a wrapper for the univariate (non-homogeneous Poisson process model) and bivariate exponential dependence model.
For the latter, the user can select either the rate or inverse rate parameter  (the inverse rate parametrization  works better for uniformity
of the p-value distribution under the \code{LR} test.

There are two options for the bivariate diagnostic: either provide pairwise minimum of marginally
exponentially distributed margins or provide a \code{n} times 2 matrix with the original data, which
is transformed to exponential margins using the empirical distribution function.
}
\examples{
\dontrun{
set.seed(123)
W.diag(rexp(1000), model="nhpp", k=30, q1=0)
# Parameter Stability only
W.diag(abs(rnorm(5000)), model="nhpp", k=30, q1=0, plots=c("PS"))
library(mvtnorm)
xbvn<-rmvnorm(6000, sigma=matrix(c(1,0.7,0.7,1),2,2))
# Transform margins to exponential manually
xbvn.exp<- -log(1-pnorm(xbvn))
W.diag(apply(xbvn.exp,1,min), model="exp", k=30, q1=0) #rate parametrization
W.diag(xbvn, model="exp", k=30, q1=0)
W.diag(apply(xbvn.exp,1,min), model="invexp", k=30, q1=0) #inverse rate parametrization
}
\dontrun{
library(ismev)
data(rain)
u <- quantile(rain, seq(0.85,0.99,by=0.01))
W.diag(xdat=rain, u=u, plots="PS")
}
}
\author{
Jennifer L. Wadsworth
}
