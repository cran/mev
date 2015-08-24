// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// Sample an index. Simpler version of sample with equal weights
// @param d positive integer
//
// @return an integer between 0 and \eqn{d-1}
int sampleone(int d){
  NumericVector index(1);
  index[0] = (double)d *runif(1)[0];
  index[0] = floor(index)[0];
  return (int)index[0];
}


//' Random variate generation for Dirichlet distribution on \eqn{S_{d}}{Sd}
//'
//' A function to sample Dirichlet random variables, based on the representation as ratios of Gamma.
//' Note that the RNG will generate on the full simplex and the sum to one constraint is respected
//' here
//'
//' @param n sample size
//' @param alpha vector of parameter
//' @param normalize boolean. If \code{FALSE}, the function returns Gamma variates with parameter \code{alpha}.
//'
//' @return sample of dimension \code{d} (size of alpha) from the Dirichlet distribution.
//' @examples rdir(n=100, alpha=c(0.5,0.5,2),TRUE)
//' rdir(n=100, alpha=c(3,1,2),FALSE)
// [[Rcpp::export]]
NumericMatrix rdir(int n, NumericVector alpha, bool normalize = true){
  NumericMatrix sample(n, alpha.size());
  for(int j=0; j<alpha.size(); j++){
    sample(_,j) = rgamma(n, alpha[j], 1.0);
  }
  if(normalize == true){
    for(int i=0; i<n; i++){
      sample(i,_) = sample(i,_)/sum(sample(i,_));
    }
  }
  return sample;
}


//' Multivariate Normal distribution sampler
//'
//' Sampler derived using the eigendecomposition of the covariance
//' matrix \code{Sigma}. The function uses the Armadillo random normal generator
//'
//' @param n sample size
//' @param mu mean vector. Will set the dimension
//' @param Sigma a square covariance matrix, of same dimension as \code{mu}.
//' No sanity check is performed to validate that the matrix is p.s.d., so use at own risk
//'
//' @return an \code{n} sample from a multivariate Normal distribution
//' @examples
//' mvrnorm(n=10, mu=c(0,2), Sigma=diag(2))
// [[Rcpp::export]]
NumericMatrix mvrnorm(int n, NumericVector mu, NumericMatrix Sigma){
  if (Sigma.nrow()!=Sigma.ncol() || mu.size()!=Sigma.ncol()){
    Rcpp::stop("Incompatible arguments - mvrnorm");
  }
  int length=Sigma.nrow();
  arma::colvec Mu(mu.begin(),length,false);
  arma::mat Xmat(Sigma.begin(), length, length, false);
  arma::mat q = arma::randn(length,arma::as_scalar(n));
  arma::colvec eigval;
  arma::mat eigvec;
  //Covariance matrix must be symmetric - otherwise eig_sym throws error
  arma::eig_sym(eigval, eigvec, Xmat);
  arma::mat sample(n,length);
  arma::colvec intermed(length);
  for(int t=0;t<n;t++){
    intermed= eigvec*arma::diagmat(arma::sqrt(eigval))*q.col(t)+Mu;
    for(int i=0;i<length;i++){
      sample(t,i)=intermed(i);
    }
  }
  return Rcpp::as<Rcpp::NumericMatrix>(wrap(sample));
}
//' Multivariate Normal distribution sampler (Rcpp version), derived using the eigendecomposition
//' of the covariance matrix Sigma. The function utilizes the arma random normal generator
//'
//' @param n sample size
//' @param mu mean vector. Will set the dimension
//' @param Sigma covariance matrix, of same dimension as mu (and square matrix).
//' No sanity check is performed to validate that the matrix is p.s.d., so use at own risk
//'
//' @return an n sample from a multivariate Normal distribution
//'
// [[Rcpp::export(.mvrnorm_arma)]]
arma::mat mvrnorm_arma(int n, arma::colvec Mu, arma::mat Xmat){
  int length = Xmat.n_rows;
  arma::mat q = arma::randn(length,arma::as_scalar(n));
  arma::colvec eigval;
  arma::mat eigvec;
  //Covariance matrix must be symmetric - otherwise eig_sym throws error
  arma::eig_sym(eigval, eigvec, Xmat);
  arma::mat sample(n,length);
  arma::colvec intermed(length);
  for(int t=0;t<n;t++){
    intermed= eigvec*arma::diagmat(arma::sqrt(eigval))*q.col(t)+Mu;
    for(int i=0;i<length;i++){
      sample(t,i)=intermed(i);
    }
  }
  return sample;
}

// DISTRIBUTIONS OF EXTREMAL FUNCTION

//' Generate from logistic \eqn{Y \sim {P_x}}, where
//' \eqn{P_{x}} is probability of extremal function scaled by a Frechet variate
//'
//' @param d dimension of the 1-sample
//' @param index index of the location. An integer in {0, ..., \eqn{d-1}}
//' @param theta a one-dimensional parameter for the logistic model, strictly greater than 1.
//'
//' @return a \code{d}-vector from \eqn{P_x}
//[[Rcpp::export(.rPlog)]]
NumericVector rPlog (int d, int index, NumericVector theta){
  if(theta[0] < 1){
    Rcpp::stop("Invalid value for the logistic model");
  }
  double shape = theta[0];
  double scale = 1/tgamma(1.0-1.0/theta[0]);
  NumericVector F(d);
  NumericVector F0(1);
  F0[0] = exp(-log(rgamma(1,1.0-1.0/theta[0],1.0)[0])/theta[0]);
  F = exp(-log(Rcpp::rexp(d,1.0))/shape)/F0[0];//*scale cancels out //scale * pow(rexp(1, 1)[0],-shape);
  F[index]=1.0;
  //*scale, but cancel out
  return F;
}


//' Generate from negative logistic \eqn{Y \sim {P_x}}, where
//' \eqn{P_{x}} is probability of extremal function scaled by a Frechet variate
//'
//' @param d dimension of the 1-sample
//' @param index index of the location. An integer in {0, ..., \eqn{d-1}}
//' @param theta a one-dimensional parameter for the negative logistic model
//'
//' @return a \code{d}-vector from \eqn{P_x}
//[[Rcpp::export(.rPneglog)]]
NumericVector rPneglog (int d, int index, NumericVector theta){
  if(theta[0] <= 0){
    Rcpp::stop("Invalid value for the negative logistic model");
  }
  NumericVector W = Rcpp::rweibull(d, theta[0], 1.0/tgamma(1.0+1.0/theta[0]));
  NumericVector Wj0 = NumericVector::create(exp(log(rgamma(1,1.0+1.0/theta[0])[0])/theta[0])/tgamma(1.0+1.0/theta[0]));
  W = W/Wj0[0];
  W[index] = 1.0;
  return W;
}

//' Generate from extremal Dirichlet \eqn{Y \sim {P_x}}, where
//' \eqn{P_{x}} is probability of extremal functions from a Dirichlet mixture
//'
//' @param d dimension of the 1-sample
//' @param index index of the location. An integer in {0, ..., \eqn{d-1}}
//' @param alpha a \eqn{d \times n} dimensional vector of positive parameter values for the Dirichlet vector
//' @param weight a \code{m} vector of mixture weights, which sum to 1
//' @return a \code{d}-vector from \eqn{P_x}
//[[Rcpp::export(.rPdirmix)]]
NumericVector rPdirmix (int d, int index, NumericMatrix alpha, NumericVector weight){
  IntegerVector int_seq = seq_len(d) - 1;
  //Probability weights
  NumericVector w (weight.size());
  for(int k=0; k<weight.size(); k++){
    w[k] = weight.size()*weight[k]*alpha(index,k)/sum(alpha(_,k));
  }
  //Rf_PrintValue(w);
  //Sample an index in 0, ..., m-1 with probabilities given by above
  IntegerVector m = RcppArmadillo::sample(int_seq, 1, false, w);
  //Define containers for the random variables for Pj
  NumericVector G(d);
  //G_{j0} variable, for location j_0, from mth mixture component
  NumericVector G0 = rgamma(1, alpha(index,m[0]) + 1.0, 1.0);
  for(int j = 0; j < d; j++){
    G[j] = rgamma(1, alpha(j,m[0]), 1.0)[0]/G0[0];
  }
  G[index] = 1.0; //Resetting the j0 value to 1.0
  return G;
}

//' Generate from bilogistic \eqn{Y \sim {P_x}}, where
//' \eqn{P_{x}} is probability of extremal functions
//'
//' @param d dimension of the 1-sample
//' @param index index of the location. An integer in {0, ..., \eqn{d-1}}
//' @param alpha a \eqn{d} dimensional vector of positive parameter values for the Dirichlet vector
//' @return a \code{d}-vector from \eqn{P_x}
//[[Rcpp::export(.rPbilog)]]
NumericVector rPbilog(int d, int index, NumericVector alpha){
  NumericVector alpha_star = rep(1.0, d);
  NumericVector sample(d);
  alpha_star[index] = 1.0-alpha[index];
  sample = rdir(1, alpha_star, true)(0,_);
  for(int i=0; i<d; i++){
    sample[i] = exp(-alpha[i]*log(sample[i])+lgamma(d-alpha[i])-lgamma(1-alpha[i]));
  }
  sample = sample/sample[index];
  return sample;
}

//' Generate from extremal Student-t \eqn{Y \sim {P_x}}, where
//' \eqn{P_{x}} is probability of extremal function
//'
//' @param index index of the location. An integer in {0, ..., \eqn{d-1}}
//' @param Sigma a positive semi-definite covariance matrix with unit variance
//' @param al the alpha parameter in Proposition 7. Corresponds to degrees of freedom - 1
//'
//' @return a \code{d}-vector from \eqn{P_x}
// [[Rcpp::export(.rPexstud)]]
NumericVector rPexstud (int index, arma::mat sigma, NumericVector al){
  if(al[0]<0 || index<0 || index >= sigma.n_cols) Rcpp::stop("Invalid argument in rPexstud");
  arma::vec zeromean = arma::vec(sigma.n_cols-1);// b/c need constructor, then setter
  zeromean.zeros(); // set elements of vector to zero
  arma::mat Covar = (sigma - sigma.col(index) * sigma.row(index))/(al[0]+1.0);
  //Covar matrix is not positive definite; shred it
  Covar.shed_row(index); Covar.shed_col(index);
  //Sample from d-1 dimensional normal
  arma::vec normalsamp = mvrnorm_arma(1, zeromean, Covar).row(0).t();
  //Add the missing zero entry back
  arma::vec indexentry = arma::vec(1);
  indexentry.zeros();
  normalsamp.insert_rows(index, indexentry);
  double nu = Rcpp::rchisq(1,al[0])[0];
  arma::vec studsamp = exp(0.5*(log(al[0])-log(nu)))*normalsamp+sigma.col(index);
  //Note: this is the shifted Student as gamma mixture,
  // i.e. adding the noncentrality parameter after multiplication by sqrt(dof)
  NumericVector samp = Rcpp::as<Rcpp::NumericVector>(wrap(studsamp));
  samp = pow(pmax(samp,0),al[0]);
  samp[index] = 1.0; //Sometimes off due to rounding
  return samp;
}

//' Generate from extremal Husler-Reiss distribution (Brown-Resnick) \eqn{Y \sim {P_x}}, where
//' \eqn{P_{x}} is probability of extremal function
//'
//' @param index index of the location. An integer in {0, ..., \eqn{d-1}}
//' @param Sigma a positive semi-definite covariance matrix
//'
//' @return a \code{d}-vector from \eqn{P_x}
//[[Rcpp::export(.rPHuslerReiss)]]
NumericVector rPHuslerReiss (int index, NumericMatrix Sigma){
  if(index<0 || index >= Sigma.ncol()) Rcpp::stop("Invalid argument in rPHuslerReiss");
  NumericVector mu(Sigma.ncol());
  NumericMatrix mvnormsamp = mvrnorm(1, mu, Sigma);
  NumericVector samp(Sigma.ncol());
  for(int i=0; i < Sigma.ncol(); i++){
    samp[i] = exp(mvnormsamp(0,i)-mvnormsamp(0,index)-0.5*(Sigma(i,i)+
      Sigma(index,index)-2*Sigma(i,index)));
  }
  return samp;
}


// SPECTRAL DISTRIBUTIONS

//' Generates from \eqn{Q_i}{Qi}, the spectral measure of the logistic model
//'
//' Simulation algorithm of Dombry et al. (2015)
//'
//' @param n sample size
//' @param theta a one-dimensional parameter
//'
//' @references Dombry, Engelke and Oesting (2015). Exact simulation of max-stable
//' processes, \emph{arXiv:1506.04430v1}, 1--24.
//'
//' @return an \code{n} by \code{d} sample from the spectral distribution
//[[Rcpp::export(.rlogspec)]]
NumericMatrix rlogspec (int n, int d, NumericVector theta){
  double shape = theta[0];
  double scale = 1/tgamma(1.0-1.0/theta[0]);
  //Define containers
  NumericMatrix samp(n,d);
  int j;
  for(int r=0; r<n; r++){
    j = sampleone(d);
    samp(r,_) = exp(-log(Rcpp::rexp(d,1.0))/shape);
    samp(r,j) = exp(-log(rgamma(1,1.0-1.0/theta[0],1.0)[0])/theta[0]);
    samp(r,_) = samp(r,_)/sum(samp(r,_));
  }
  return samp;
}




//' Generates from \eqn{Q_i}{Qi}, the spectral measure of the negative logistic model
//'
//' Simulation algorithm of Dombry et al. (2015)
//'
//' @param n sample size
//' @param theta a one-dimensional parameter
//'
//' @references Dombry, Engelke and Oesting (2015). Exact simulation of max-stable
//' processes, \emph{arXiv:1506.04430v1}, 1--24.
//'
//' @return an \code{n} by \code{d} sample from the spectral distribution
//[[Rcpp::export(.rneglogspec)]]
NumericMatrix rneglogspec (int n, int d, NumericVector theta){
  NumericMatrix samp(n,d);
  int j;
  for(int r=0; r<n; r++){
    j = sampleone(d);
    samp(r,_) = Rcpp::rweibull(d, theta[0], 1.0/tgamma(1.0+1.0/theta[0]));
    samp(r,j) = exp(log(rgamma(1,1.0+1.0/theta[0])[0])/theta[0])/tgamma(1.0+1.0/theta[0]);
    samp(r,_) = samp(r,_)/sum(samp(r,_));
  }
  return samp;
}


//' Generates from \eqn{Q_i}{Qi}, the spectral measure of the Dirichlet mixture model
//'
//' Simulation algorithm of Dombry et al. (2015)
//'
//' @param n sample size
//' @param d dimension of the 1-sample
//' @param alpha a \eqn{d \times n} dimensional vector of positive parameter values for the Dirichlet vector
//' @param weight a \code{m} vector of mixture weights, which sum to 1
//'
//' @return an \code{n} by \code{d} sample from the spectral distribution
//[[Rcpp::export(.rdirmixspec)]]
NumericMatrix rdirmixspec (int n, int d, NumericMatrix alpha, NumericVector weight){
  NumericMatrix samp(n, d);
  int j;
  IntegerVector int_seq = seq_len(d) - 1;
  //Probability weights
  NumericVector w (weight.size());
  for(int r=0; r<n; r++){
    j = sampleone(d);
    for(int k=0; k<weight.size(); k++){
      w[k] = weight.size()*weight[k]*alpha(j,k)/sum(alpha(_,k));
    }
    //Rf_PrintValue(w);
    //Sample an index in 1, ..., m with probabilities given by above
    IntegerVector m = RcppArmadillo::sample(int_seq, 1, false, w);
    //Modify the mth entry to match the extremal Dirichlet model
    NumericVector alpha_sub = alpha(_,m[0]);
    alpha_sub[j] = alpha_sub[j] + 1.0;
    //Define containers for the random variables for Pj
    //G_{j0} variable, for location j_0, from mth mixture component
    for(int k = 0; k < d; k++){
      samp(r,k) = rgamma(1, alpha(k,m[0]), 1.0)[0];
    }
    samp(r,j) = rgamma(1, alpha(j, m[0]) + 1.0, 1.0)[0];
    samp(r,_) = samp(r,_)/sum(samp(r,_));
  }
  return samp;
}


//' Generates from \eqn{Q_i}{Qi}, the spectral measure of the bilogistic model
//'
//' Simulation algorithm of Boldi (2009) for the bilogistic model
//'
//' @param n sample size
//' @param alpha vector of parameter of dimension \code{d}
//'
//' @references Boldi (2009). A note on the representation of parametric models
//' for multivariate extremes. \emph{Extremes} \bold{12}, 211--218.
//'
//' @return an \code{n} by \code{d} sample from the spectral distribution
// [[Rcpp::export(.rbilogspec)]]
NumericMatrix rbilogspec(int n, NumericVector alpha){
  NumericMatrix sample(n,alpha.size());
  NumericVector alpha_star = rep(1.0, alpha.size());
  int j;
  for(int r=0; r<n; r++){
    j = sampleone(alpha.size());
    alpha_star[j] = 1.0-alpha[j];
    sample(r,_) = rdir(1, alpha_star, true)(0,_);
    for(int i=0; i<alpha.size(); i++){
      sample(r,i) = exp(-alpha[i]*log(sample(r,i))+lgamma(alpha.size()-alpha[i])-lgamma(1-alpha[i]));
    }
    sample(r,_) = sample(r,_)/sum(sample(r,_));
    //Reset value for multiple samples
    alpha_star[j] = 1.0;
  }
  return sample;
}

//' Generates from \eqn{Q_i}{Qi}, the spectral measure of the extremal Student model
//'
//' @param index index of the location. An integer in {0, ..., \eqn{d-1}}
//' @param Sigma a positive semi-definite covariance matrix with unit variance
//' @param al the alpha parameter in Proposition 7. Corresponds to degrees of freedom - 1
//'
//' @return an \code{n} by \code{d} sample from the spectral distribution
// [[Rcpp::export(.rexstudspec)]]
NumericMatrix rexstudspec (int n, arma::mat sigma, NumericVector al){
  if(al[0]<0) Rcpp::stop("Invalid dof argument in rexstudspec");
  //Define containers and auxiliary variables
  arma::vec zeromean = arma::vec(sigma.n_cols-1);// b/c need constructor, then setter
  zeromean.zeros(); // set elements of vector to zero
  int d = sigma.n_cols;
  NumericMatrix samp(n,d);
  int j;
  arma::mat Covar = arma::mat(sigma.n_rows,sigma.n_cols);
  //Need to adjust the size of Covar because it is shreded
  arma::vec normalsamp = arma::vec(d-1);
  arma::vec indexentry = arma::vec(1);
  arma::vec studsamp = arma::vec(d);
  indexentry.zeros();
  double nu;
  for(int r=0; r<0; r++){
    j = sampleone(d);
    //Redefine values
    Covar = arma::mat(sigma.n_rows,sigma.n_cols);
    normalsamp = arma::vec(d-1);
    Covar = (sigma - sigma.col(j) * sigma.row(j))/(al[0]+1.0);
    //Covar matrix is not positive definite; shred it
    Covar.shed_row(j); Covar.shed_col(j);
    //Sample from d-1 dimensional normal
    normalsamp = mvrnorm_arma(1, zeromean, Covar).row(0).t();
    normalsamp.insert_rows(j, indexentry);
    nu = Rcpp::rchisq(1,al[0])[0];
    studsamp = exp(0.5*(log(al[0])-log(nu)))*normalsamp+sigma.col(j);
    //Note: this is the shifted Student as gamma mixture,
    // i.e. adding the noncentrality parameter after multiplication by sqrt(dof)
    samp(r,_) = Rcpp::as<Rcpp::NumericVector>(wrap(studsamp));
    samp(r,_) = pow(pmax(samp(r,_),0),al[0]);
    samp(r,j) = 1.0; //Sometimes off due to rounding
    samp(r,_) = samp(r,_)/sum(samp(r,_));
  }
  return samp;
}

//' Generates from \eqn{Q_i}{Qi}, the spectral measure of the Brown-Resnick (or Husler-Reiss) model
//'
//' Simulation algorithm of Dombry et al. (2015)
//'
//' @param n sample size
//' @param Sigma \code{d}-dimensional covariance matrix
//'
//' @references Dombry, Engelke and Oesting (2015). Exact simulation of max-stable
//' processes, \emph{arXiv:1506.04430v1}, 1--24.
//'
//' @return an \code{n} by \code{d} sample from the spectral distribution
// [[Rcpp::export(.rhrspec)]]
NumericMatrix rhrspec (int n, NumericMatrix Sigma){
  int d = Sigma.ncol();
  NumericVector mu(d);
  NumericMatrix mvnormsamp = mvrnorm(n, mu, Sigma);
  NumericMatrix samp(n, d);
  int j;
  for(int r=0; r<n; r++){
    j = sampleone(d);
    for(int i=0; i < d; i++){
      samp(r,i) = exp(mvnormsamp(0,i)-mvnormsamp(0,j)-0.5*(Sigma(i,i)+
        Sigma(j,j)-2*Sigma(i,j)));
    }
    samp(r,_) = samp(r,_)/sum(samp(r,_));
  }
  return samp;
}

//' Generates from \eqn{Q_i}{Qi}, the spectral measure of the extremal Dirichlet
//' model
//'
//' This model was introduced in Coles and Tawn (1991); the
//' present method uses the simulation algorithm of Boldi (2009) for the extremal Dirichlet model
//'
//' @param n sample size
//' @param d dimension of sample
//' @param alpha vector of Dirichlet parameters of dimension \code{d}, or \eqn{d+1} vector with the \code{d} Dirichlet parameters and an index of regular variation in \eqn{[0, 1]}
//' @param rho index of regular variation
//' @param irv should the usual model (\code{FALSE}) or the general scaled version (\code{TRUE}) be used
//'
//' @references Boldi (2009). A note on the representation of parametric models
//' for multivariate extremes. \emph{Extremes} \bold{12}, 211--218.
//'
//' @return an \code{n} by \code{d} sample from the spectral distribution
// [[Rcpp::export(.rdirspec)]]
NumericMatrix rdirspec(int n, int d, NumericVector alpha, bool irv = false){
	NumericVector alpha_star(d);
  int j;
  NumericMatrix sample(n, d);
  NumericVector m(d);
  if(irv==true){
	   for(int i=0; i<d; i++){
	   m[i] = -lgamma(alpha[i])+lgamma(alpha[i]+alpha[d]);
	  	alpha_star[i]=alpha[i];
	  }
  } else{
  alpha_star = clone(alpha);
  }
  for(int r=0; r<n; r++){
    j = sampleone(d);
    if(irv==false){
	    alpha_star[j] = alpha_star[j]+1.0;
	    sample(r,_) = rdir(1, alpha_star, false)(0,_);
	    for(int i=0; i<d; i++){
	      sample(r,i) = sample(r,i)/alpha[i];
	    }
    } else{
    	alpha_star[j] += alpha[d];
	    sample(r,_) = rdir(1, alpha_star, false)(0,_);
	    for(int i=0; i<d; i++){
	      sample(r,i) = exp(alpha[d]*log(sample(r,i))-m[i]);
	    }
    }
    alpha_star[j] = alpha[j];
    sample(r,_) = sample(r,_)/sum(sample(r,_));
    
  }
  return sample;
}


//' Multivariate extreme value distribution sampling algorithm
//' corresponding to Algorithm 1 in Dombry, Engelke and Oesting (2015),
//' using the formulation of the Dirichlet mixture of Coles and Tawn (1991)
//' as described and derived in Boldi (2009) for the bilogistic and extremal
//' Dirichlet model. Models currently implemented include logistic, negative
//' logistic, extremal Dirichlet and bilogistic MEV.
//'
//' @param n sample size
//' @param d dimension of the multivariate distribution
//' @param param a vector of parameters
//' @param model integer, currently ranging from 1 to 4, corresponding respectively to
//' (1) \code{log}, (2) \code{neglog}, (3) \code{dirmix}, (4) \code{bilog},
//' (5) \code{extstud}, (6) \code{hr} and (7) \code{ct}
//' @param Sigma covariance matrix for Husler-Reiss and extremal student. Default for compatibility
//'
//' @return a \code{n} by \code{d} matrix containing the sample
// [[Rcpp::export(.rmevA1)]]
NumericMatrix rmevA1(int n, int d, NumericVector param, int model, NumericMatrix Sigma) {
  //Sanity checks
  arma::mat sigma(Sigma.begin(), Sigma.nrow(), Sigma.ncol(), false);
	bool irv = false;
	//Model 1: logistic
  if(model == 1 && param.size() != 1){
    Rcpp::warning("Logistic model currently only implemented for one argument");
  }
  //Model 2: negative logistic
  if(model == 2 && param.size() != 1){
    Rcpp::warning("Negative logistic model currently only implemented for one argument");
  }
  //Model 3: Dirichlet mixture
  //Model 4: bilogistic
  if(model == 4){
    if(param.size() != d || is_true(any(param > 1.0))){
         Rcpp::stop("Invalid input for the bilogistic or the negative bilogistic model");
    }
  }
  //Model 5: Extremal Student t
  if(model == 5){
    if(Sigma.ncol()!=Sigma.nrow()) Rcpp::stop("Provided covariance matrix is not square");
    if(param[0]<0 || param.size()!=1){Rcpp::stop("Invalid degree of freedom");}
    //Standardize the covariance to correlation matrix (do only once)
    arma::vec stdev = exp(0.5*log(sigma.diag()));
    arma::mat stdevmat = inv(diagmat(stdev));
    sigma = stdevmat * sigma * stdevmat;
  }
  //Model 6: Brown-Resnick (Husler-Reiss)
  if(model == 6){
    if(Sigma.ncol()!=Sigma.nrow()) Rcpp::stop("Provided covariance matrix is not square");
  }
  //Model 7: Coles and Tawn (Beta-Dirichlet extremal distribution)
  if(model == 7){
    if(is_true(any(param < 0.0))){
      Rcpp::stop("Invalid input for Dirichlet model");
    }
    if(param.size() == d+1){
    	if(param[d]>1.0) Rcpp::stop("Invalid index of regular variation");
    	irv = true;	
    }
  }

  NumericMatrix samp = NumericMatrix(n, d); //Initialized to zero
  NumericVector zeta_I(1);
  NumericVector Y(d);
  for(int i = 0; i < n; i ++){
    //For each sample of the max-stable distribution
    zeta_I[0] = rexp(1, d)[0];
    while(1.0/zeta_I[0] > min(samp( i, _ ))){
      //Simulate from T and Y, or spectral density D
      int T = sampleone(d);
      if(model == 1){
        Y = rPlog(d, T, param);
      } else if(model == 2){
        Y = rPneglog(d, T, param);
      } else if(model == 3){
        Y = rPdirmix(d, T, Sigma, param);
      } else if(model == 4){
        Y = rbilogspec(1, param)(0,_);
      } else if(model == 5){
        Y = rPexstud(T, sigma, param);
      } else if(model == 6){
        Y = rPHuslerReiss(T, Sigma);
      } else if(model ==7){
        Y = rdirspec(1, d, param, irv)(0,_);
      } else {
        Rcpp::stop("Model not yet implemented");
      }
      // Generating from poisson point process and
      // computing pointwise maxima
      if(model == 7 || model == 4){
        samp(i,_) = pmax(samp(i, _ ), Y/zeta_I[0]);
      } else{
        samp(i,_) = pmax(samp(i, _ ), Y/(sum(Y)*zeta_I[0]));
      }
      zeta_I[0] = zeta_I[0] + rexp(1,d)[0];//Rcpp uses scale rather than rate
    }

  }
  return samp;
}


//' Multivariate extreme value distribution sampling algorithm
//' corresponding to Algorithm 2 in Dombry, Engelke and Oesting (2015)
//'
//' Models currently implemented include logistic and negative logistic, sampling
//' from the extremal functions. This requires derivation of \eqn{P_x}
//'
//' @param n sample size
//' @param d dimension of the multivariate distribution
//' @param param a vector of parameters
//' @param model integer, currently ranging from 1 to 7, corresponding respectively to
//' (1) \code{log}, (2) \code{neglog}, (3) \code{dirmix}, (4) \code{bilog},
//' (5) \code{extstud} and (6) \code{hr}.
//' @param Sigma covariance matrix for Husler-Reiss and extremal student. Default for compatibility
//'
//' @return a \code{n} by \code{d} matrix containing the sample
// [[Rcpp::export(.rmevA2)]]
NumericMatrix rmevA2(int n, int d, NumericVector param, int model, NumericMatrix Sigma) {
  arma::mat sigma(Sigma.begin(), Sigma.nrow(), Sigma.ncol(), false);
	bool irv = false;
  //Sanity checks
  //Model 1: logistic
  if(model==1 && param.size()!=1){
    Rcpp::warning("Logistic model currently only implemented for one argument");
  }
  //Model 2: negative logistic
  if(model==2 && param.size()!=1){
    Rcpp::warning("Negative logistic model currently only implemented for one argument");
  }
  //Model 3: Dirichlet mixture
  //Checks are performed in rmev wrapper function
  //Model 4: bilogistic
  if(model == 4){
    if(param.size() != d || is_true(any(param > 1.0))){
      Rcpp::stop("Invalid input for the bilogistic or the negative bilogistic model");
    }
  }
  //Model 5: Extremal Student t
  if(model == 5){
    if(Sigma.ncol()!=Sigma.nrow()) Rcpp::stop("Provided covariance matrix is not square");
    if(param[0]<0 || param.size()!=1){Rcpp::stop("Invalid degree of freedom");}
    //Standardize the covariance to correlation matrix (do only once)
    arma::vec stdev = exp(0.5*log(sigma.diag()));
    arma::mat stdevmat = inv(diagmat(stdev));
    sigma = stdevmat * sigma * stdevmat;
  }
  //Model 6: Brown-Resnick (Husler-Reiss)
  if(model == 6){
    if(Sigma.ncol()!=Sigma.nrow()) Rcpp::stop("Provided covariance matrix is not square");
  }
  if(model > 6){
    Rcpp::stop("Model not implemented yet");
  }

  //Define the containers
  NumericMatrix samp = NumericMatrix(n, d);
  NumericVector zeta_I(1);
  NumericVector Y(d);
  for(int i = 0; i < n; i ++){
    //For each sample of the max-stable distribution
    zeta_I[0] = rexp(1, 1)[0];   //(1) Initial sample
    if(model == 1){
      Y = rPlog(d, 0, param);
    } else if(model == 2){
      Y = rPneglog(d, 0, param);
    } else if(model == 3){
      Y = rPdirmix(d, 0, Sigma, param);
    } else if(model == 4){
      Y = rPbilog(d, 0, param);
    } else if(model == 5){
      Y = rPexstud(0, sigma, param);
    } else if(model == 6){
      Y = rPHuslerReiss(0, Sigma);
    }  else{
      Rcpp::stop("Sampler not yet implemented with extremal functions");
    }

    samp(i,_) = Y/zeta_I[0]; //(2) Set starting values
    for(int j = 1; j < d; j++){//(3) per coordinate
      zeta_I[0] = rexp(1, 1.0)[0]; //(4) Poisson process generation
      //Rcpp::Rcout << "Extremal path until " << j << "has been generated" << std::endl;
      while(1.0/zeta_I[0] > samp( i, j )){ //(5) Check stopping rule
        //(6)  Simulate from Pxn
        if(model == 1){
          Y = rPlog(d, j, param);
        } else if(model == 2){
          Y = rPneglog(d, j, param);
        } else if(model == 3){
          Y = rPdirmix(d, j, Sigma, param);
        } else if(model == 4){
          Y = rPbilog(d, j, param);
        } else if(model == 5){
          Y = rPexstud(j, sigma, param);
        } else if(model == 6){
          Y = rPHuslerReiss(j, Sigma);
        }

        bool res = true;
        for(int k = 0; k < j; k++){ //(7) Check previous extremal functions
          //Rcpp::Rcout << "Extremal path until " << j << "has been generated, with " << k<< std::endl;
          if(Y[k]/zeta_I[0] >= samp(i,k)){
            res = false;
            break;
          }
        }
        if(res){ //(8) Update if true
          samp(i, _ ) = pmax(samp(i, _ ), Y/zeta_I[0]);
        }
        zeta_I[0] = zeta_I[0] + rexp(1, 1.0)[0]; //(9) Increment of Poisson process
      }
    }
  }
  return samp;
}


//' Random number generator from spectral distribution
//'
//' Generate from \eqn{Q_i}{Qi}, the spectral measure of a given multivariate extreme value model
//'
//' @param n sample size
//' @param d dimension of the multivariate distribution
//' @param param a vector of parameters
//' @param model integer, currently ranging from 1 to 7, corresponding respectively to
//' (1) \code{log}, (2) \code{neglog}, (3) \code{dirmix}, (4) \code{bilog},
//' (5) \code{extstud}, (6) \code{hr} and (7) \code{ct}.
//' @param Sigma covariance matrix for Husler-Reiss and extremal student. Default for compatibility
//'
//' @references Dombry, Engelke and Oesting (2015). Exact simulation of max-stable processes, \emph{arXiv:1506.04430v1}, 1--24.
//' @references Boldi (2009). A note on the representation of parametric models for multivariate extremes. \emph{Extremes} \bold{12}, 211--218.
//'
//' @return a \code{n} by \code{d} matrix containing the sample
// [[Rcpp::export(.rmevspec_cpp)]]
NumericMatrix rmevspec_cpp(int n, int d, NumericVector param, int model, NumericMatrix Sigma) {
  bool irv = false;
  arma::mat sigma(Sigma.begin(), Sigma.nrow(), Sigma.ncol(), false);
	//Sanity checks
  //Model 1: logistic
  if(model == 1 && param.size() != 1){
    Rcpp::warning("Logistic model currently only implemented for one argument");
  }
  //Model 2: negative logistic
  if(model == 2 && param.size() != 1){
    Rcpp::warning("Negative logistic model currently only implemented for one argument");
  }
  //Model 3: Dirichlet mixture
  //Model 4: bilogistic
  if(model == 4){
    if(param.size() != d || is_true(any(param > 1.0))){
     Rcpp::stop("Invalid input for the bilogistic or the negative bilogistic model");
    }
  }
  //Model 5: Extremal Student t
  if(model == 5){
    if(Sigma.ncol()!=Sigma.nrow()) Rcpp::stop("Provided covariance matrix is not square");
    if(param[0]<0 || param.size()!=1){Rcpp::stop("Invalid degree of freedom");}
    //Standardize the covariance to correlation matrix (do only once)
    arma::vec stdev = exp(0.5*log(sigma.diag()));
    arma::mat stdevmat = inv(diagmat(stdev));
    sigma = stdevmat * sigma * stdevmat;
  }
  //Model 6: Brown-Resnick (Husler-Reiss)
  if(model == 6){
    if(Sigma.ncol()!=Sigma.nrow()) Rcpp::stop("Provided covariance matrix is not square");
  }
  //Model 7: Coles and Tawn (Beta-Dirichlet extremal distribution)
  if(model == 7){
    if(is_true(any(param < 0.0))){
      Rcpp::stop("Invalid input for Dirichlet model");
    }
    if(param.size()==d+1){
    	irv = true;
  	}
  } else if (model>7){
  	Rcpp::stop("Model not currently implemented");
  }
	//Sampling
  
  NumericMatrix samp = NumericMatrix(n, d); //Initialized to zero
  if(model == 1){
    samp = rlogspec(n, d, param);
  } else if(model == 2){
    samp = rneglogspec(n, d, param);
  } else if(model == 3){
    samp = rdirmixspec(n, d, Sigma, param);
  } else if(model == 4){
    samp = rbilogspec(n, param);
  } else if(model == 5){
    samp = rexstudspec(n, sigma, param);
  } else if(model == 6){
    samp = rhrspec(n, Sigma);
  } else if(model == 7){
    samp = rdirspec(n, d, param, irv);
  }  else{
    Rcpp::stop("Invalid model");
  }
  return samp;
}
