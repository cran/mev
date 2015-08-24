#' Exact simulations of multivariate extreme value distributions
#'
#' Implementation of the random number generators for multivariate extreme-value distributions
#' and max-stable processes based on the two algorithms described in
#' Dombry, Engelke and Oesting (2015).
#'
#'@param n number of observations
#'@param d dimension of sample
#'@param param parameter vector for the logistic, bilogistic, negative bilogistic and extremal Dirichlet (Coles and Tawn) model.
#' Parameter matrix for the Dirichlet mixture. Degree of freedoms for extremal student model.
#'@param sigma covariance matrix for Husler-Reiss and extremal Student-t distributions
#'@param alg algorithm, either simulation via extremal function or via the spectral measure. The extremal Dirichlet model is only implemented with \code{sm}.
#'@param model choice between 1-parameter logistic and negative logistic, bilogistic and the extremal Dirichlet model of Coles and Tawn,
#' the Brown-Resnick and extremal Student max-stable process (which generate the Husler-Reiss MEV distribution), or the Dirichlet mixture.
#'@param vario function specifying the variogram. Used only if provided in conjonction with \code{loc} and if \code{sigma} is missing
#'@param loc \code{d} by \code{k} matrix of location, used as input in the variogram \code{vario}.
#'@param weights vector of length \code{m} for the \code{m} mixture components. Must sum to one
#'
#'@author Leo Belzile
#'@details The vector param differs depending on the model
#' \itemize{
#'  \item \code{log}: one dimensional parameter greater than 1
#'  \item \code{neglog}: one dimensional positive parameter
#'  \item \code{bilog}: \code{d}-dimensional vector of parameters in \eqn{[0,1]}
#'  \item \code{negbilog}: \code{d}-dimensional vector of negative parameters
#'  \item \code{ct}: \code{d}-dimensional vector of positive (a)symmetry parameters. Alternatively, a \eqn{d+1} 
#'  vector consisting of the \code{d} Dirichlet parameters and the last entry is an index of regular variation in \code{(0, 1]} treated as scale
#'  \item \code{xstud}: one dimensional parameter corresponding to degrees of freedom \code{alpha}
#'  \item \code{dirmix}: \code{d} by \code{m}-dimensional matrix of positive (a)symmetry parameters
#' }
#'@return an \code{n} by \code{d} exact sample from the corresponding multivariate extreme value model
#'@references
#'Dombry, Engelke and Oesting (2015). Exact simulation of max-stable processes, \emph{arXiv:1506.04430v1}, 1--24.
#'@examples
#'set.seed(1)
#'rmev(n=100, d=3, param=2.5, model="log", alg="ef")
#'rmev(n=100, d=4, param=c(0.2,0.1,0.9,0.5), model="bilog", alg="sm")
#'## Spatial example using variogram, from Clement Dombry
#'#Variogram gamma(h) = scale*||h||^alpha
#'scale <- 0.5; alpha <- 1
#'vario <- function(x) scale*sqrt(sum(x^2))^alpha
#'#grid specification
#'grid.loc <- as.matrix(expand.grid(runif(4), runif(4)))
#'rmev(n=100, vario=vario,loc=grid.loc, model="hr")
#'## Example with Dirichlet mixture
#'alpha.mat <- cbind(c(2,1,1),c(1,2,1),c(1,1,2))
#'rmev(n=100, param=alpha.mat, weights=rep(1/3,3), model="dirmix")
rmev <-function(n, d, param, sigma, model=c("log","neglog","bilog","negbilog","hr","xstud","ct","dirmix"),
                alg=c("ef","sm"), weights, vario, loc){
	if(!missing(param) && mode(param) != "numeric") stop("Invalid parameter")
  alg <- match.arg(alg)
  model <- match.arg(model)
	m1 <- c("log","neglog")
  m2 <- c("bilog","negbilog","ct")
  m3 <- c("hr","xstud")
  if(model %in% m1){
    d <- as.integer(d)
    sigma = cbind(0)
    if(missing(param) || param < 0 || d < 1){
      stop("Invalid parameter value")
    }
    if(length(param)!=1){
      warning("Only first entry of param vector considered")
      param <- param[1]
    }
    if(model=="log"){
      if(param < 1.0){
      param <- 1.0/param
      }
      mod <- 1
    } else{
      mod <- 2
    }
  } else if(model %in% m2){
    d <- as.integer(d)
    sigma = cbind(0)
  	if(model %in% c("bilog","negbilog")){
    	if(missing(param) || length(param)!=d) stop("Invalid parameter value")
  		#Check whether arguments are valid
  		if(model=="bilog" && all(param>=1)){ param <- 1/param}
  		if(model=="negbilog" && all(param>=0)){ param <- -param}
      if(any(param>1.0)) stop("Invalid param vector for bilogistic or negative bilogistic")
    	if(any(param<0) && model=="bilog") warning("Negative parameter values in bilogistic");
    	if(any(param>0) && model=="negbilog") warning("Positive parameter values in negative bilogistic");
      mod <- 4
    } else{
    	if(missing(param) || (length(param)!=d && length(param)!=d+1)) stop("Invalid parameter value")
    	if(alg=="ef") warning("Not implemented using extremal functions")
    		alg = "sm"
        mod <- 7
    }
  } else if(model %in% m3){
    if(missing(sigma) && !missing(vario) && !missing(loc)){
      if(!is.matrix(loc)) loc <- matrix(loc, ncol=1)
      stopifnot(is.function(vario))
      sigma <- sapply(1:nrow(loc), function(i) sapply(1:nrow(loc), function(j)
        vario(loc[i,]) + vario(loc[j,]) - vario(loc[i,]-loc[j,])))
    }
    d <- ncol(sigma)
    if(missing(sigma) || ncol(sigma)!=nrow(sigma)) stop("Invalid covariance matrix")
    if(any(diag(sigma)<=0)) stop("Degenerate covariance matrix; negative or zero entries found")
    if(model=="xstud" && any(diag(sigma)!=1)) warning("Extremal student requires correlation matrix")
    if(model=="xstud" && (missing(param) || length(param)!=1)) stop("Degrees of freedom argument missing or invalid")
    if(model=="xstud"){
      mod <- 5
    } else{
      mod <- 6; param = 0
    }
  } else if(model=="dirmix"){
    if(any(missing(param),
           length(weights)!=ncol(param) && ncol(param)!=1,
           any(param<0))){
      stop("Invalid arguments for the Dirichlet mixture")
    }
    if(!missing(weights)){
     if(any(weights<0))stop("Negative weights provided")
      if(sum(weights)!=1) warning("weights do not sum to one")
      weights <- weights/sum(weights)
    }
    if(missing(d)){ d <- nrow(param)
    } else if(d != nrow(param)){
      stop("Dimension of d and provided param do not match")
    }
    #Checking for the mean constraints
    mar_mean <- colSums(t(param)/ colSums(param)*weights)-1/d
		if(any(mar_mean!=0)) stop("Invalid mixture components")
  	#Switching parameters around to pass them to Rcpp function
  	sigma <- param
  	param <- weights
  	mod <- 3
  }
  samp <- switch(alg,
                 ef=.rmevA2(n=n, d=d, param=param, model=mod, Sigma=sigma),
                 sm=.rmevA1(n=n, d=d, param=param, model=mod, Sigma=sigma)
                )
  return(samp)
}




#' Random samples from spectral distributions of multivariate extreme value models.
#'
#' Generate from \eqn{Q_i}{Qi}, the spectral measure of a given multivariate extreme value model
#'
#'@param n number of observations
#'@param d dimension of sample
#'@param param parameter vector for the logistic, bilogistic, negative bilogistic and Dirichlet (Coles and Tawn) model.
#' Parameter matrix for the Dirichlet mixture. Degree of freedoms for extremal student model.
#'@param sigma covariance matrix for Husler-Reiss and extremal Student-t distributions
#'@param model choice between 1-parameter logistic and negative logistic, bilogistic, negative bilogistic and extremal Dirichlet,
#' the Brown-Resnick and extremal Student max-stable process (which generate the Husler-Reiss MEV distribution), or the Dirichlet mixture.
#'@param vario function specifying the variogram. Used only if provided in conjonction with \code{loc} and if \code{sigma} is missing
#'@param loc \code{d} by \code{k} matrix of location, used as input in the variogram \code{vario}.
#'@param weights vector of length \code{m} for the \code{m} mixture components. Must sum to one
#'
#'@author Leo Belzile
#'@details The vector param differs depending on the model
#' \itemize{
#'  \item \code{log}: one dimensional parameter greater than 1
#'  \item \code{neglog}: one dimensional positive parameter
#'  \item \code{bilog}: \code{d}-dimensional vector of parameters in \eqn{[0,1]}
#'  \item \code{negbilog}: \code{d}-dimensional vector of negative parameters
#'  \item \code{ct}: \code{d}-dimensional vector of positive (a)symmetry parameters. Alternatively, a \eqn{d+1} 
#'  vector consisting of the \code{d} Dirichlet parameters and the last entry is an index of regular variation in \code{(0, 1]} treated as scale
#'  \item \code{xstud}: one dimensional parameter corresponding to degrees of freedom \code{alpha}
#'  \item \code{dirmix}: \code{d} by \code{m}-dimensional matrix of positive (a)symmetry parameters
#' }
#'@return an \code{n} by \code{d} exact sample from the corresponding multivariate extreme value model
#'
#'@references Dombry, Engelke and Oesting (2015). Exact simulation of max-stable processes,
#' \emph{arXiv:1506.04430v1}, 1--24.
#'@references Boldi (2009). A note on the representation of parametric models for multivariate extremes.
#' \emph{Extremes} \bold{12}, 211--218.
#'
#' @examples
#'set.seed(1)
#'rmevspec(n=100, d=3, param=2.5, model="log")
#'rmevspec(n=100, d=3, param=2.5, model="neglog")
#'rmevspec(n=100, d=4, param=c(0.2,0.1,0.9,0.5), model="bilog") 
#'rmevspec(n=100, d=2, param=c(0.8,1.2), model="ct") #Dirichlet model
#'rmevspec(n=100, d=2, param=c(0.8,1.2,0.5), model="ct") #with additional scale parameter
#'#Variogram gamma(h) = scale*||h||^alpha
#'scale <- 0.5; alpha <- 1
#'vario <- function(x) scale*sqrt(sum(x^2))^alpha
#'#grid specification
#'grid.loc <- as.matrix(expand.grid(runif(4), runif(4)))
#'rmevspec(n=100, vario=vario,loc=grid.loc, model="hr")
#'## Example with Dirichlet mixture
#'alpha.mat <- cbind(c(2,1,1),c(1,2,1),c(1,1,2))
#'rmevspec(n=100, param=alpha.mat, weights=rep(1/3,3), model="dirmix")
rmevspec <-function(n, d, param, sigma, model=c("log","neglog","bilog","negbilog","hr","xstud","ct","dirmix"),
                     weights, vario, loc){
  if(!missing(param) && mode(param) != "numeric") stop("Invalid parameter")
  model <- match.arg(model)
  m1 <- c("log","neglog")
  m2 <- c("bilog","negbilog","ct")
  m3 <- c("hr","xstud")
  if(model %in% m1){
    d <- as.integer(d)
    sigma = cbind(0)
    if(missing(param) || param < 0 || d < 1){
      stop("Invalid parameter value")
    }
    if(length(param)!=1){
      warning("Only first entry of param vector considered")
      param <- param[1]
    }
    if(model=="log"){
      if(param < 1.0){
        param <- 1.0/param
      }
      mod <- 1
    } else{
      mod <- 2
    }
  } else if(model %in% m2){
    d <- as.integer(d)
    sigma = cbind(0)
    if(model %in% c("bilog","negbilog")){
    	if(missing(param) || length(param)!=d) stop("Invalid parameter value")
      if(model=="bilog" && all(param>=1)){ param <- 1/param}
    	if(model=="negbilog" && all(param>=0)){param <- -param}
    	if(any(param>1.0)) stop("Invalid param vector for bilogistic or negative bilogistic")
    	if(any(param<0) && model=="bilog") warning("Negative parameter values in bilogistic");
    	if(any(param>0) && model=="negbilog") warning("Positive parameter values in negative bilogistic");
      mod <- 4
    } else{
    	if(missing(param) || any(param<0) || (length(param)!=d && length(param)!=d+1)) stop("Invalid parameter value for extremal Dirichlet")
    	if(length(param)==d+1 && param[d+1]>1)  stop("Invalid parameter value")
        mod <- 7
    }
  } else if(model %in% m3){
    if(missing(sigma) && !missing(vario) && !missing(loc)){
      if(!is.matrix(loc)) loc <- matrix(loc, ncol=1)
      stopifnot(is.function(vario))
      sigma <- sapply(1:nrow(loc), function(i) sapply(1:nrow(loc), function(j)
        vario(loc[i,]) + vario(loc[j,]) - vario(loc[i,]-loc[j,])))
    }
    d <- ncol(sigma)
    if(missing(sigma) || ncol(sigma)!=nrow(sigma)) stop("Invalid covariance matrix")
    if(any(diag(sigma)<=0)) stop("Degenerate covariance matrix; negative or zero entries found")
    if(model=="xstud" && any(diag(sigma)!=1)) warning("Extremal student requires correlation matrix")
    if(model=="xstud" && (missing(param) || length(param)!=1)) stop("Degrees of freedom argument missing or invalid")
    if(model=="xstud"){
      mod <- 5
    } else{
      mod <- 6; param = 0
    }
  } else if(model=="dirmix"){
    if(any(missing(param),
           length(weights)!=ncol(param) && ncol(param)!=1,
           any(param<0))){
      stop("Invalid arguments for the Dirichlet mixture")
    }
    if(!missing(weights)){
      if(any(weights<0))stop("Negative weights provided")
      if(sum(weights)!=1) warning("weights do not sum to one")
      weights <- weights/sum(weights)
    }
    if(missing(d)){ d <- nrow(param)
    } else if(d != nrow(param)){
      stop("Dimension of d and provided param do not match")
    }
    #Checking for the mean constraints
    mar_mean <- colSums(t(param)/ colSums(param)*weights)-1/d
    if(any(mar_mean!=0)) stop("Invalid mixture components")
    #Switching parameters around to pass them to Rcpp function
    sigma <- param
    param <- weights
    mod <- 3
  }
  .rmevspec_cpp(n=n, d=d, param=param, model=mod, Sigma=sigma)
}
