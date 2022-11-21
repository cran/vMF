// [[Rcpp::depends(RcppArmadillo, BH)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>

using namespace Rcpp;
using namespace arma;
using namespace std;

// FOR SAMPLING von Mises Fisher Distribution 

// function computes the W following :
// Andrew T.A Wood (1994) Simulation of the von mises fisher distribution, 
// Communications in Statistics - Simulation and Computation

void rw(const int& size, const double& lambda, const int& d, arma::vec& W){
  // Step 0
  // Algebraically equivalent to
  // (-2. * l + sqrt(4. * l * l + d * d)) / d
  // but numerically more stable. See 
  // Hornik, K., & Grün, B. (2014). movMF: An R package for fitting mixtures
  // of von Mises-Fisher distributions. Journal of Statistical Software, 
  // 58(10), 1-31.
  double b = d/ (sqrt(4. * lambda * lambda + d*d) + 2. * lambda);
  double x = (1. - b) / (1. + b);
  double c = lambda * x + d * log(1. - x * x);
  
  // Step 1
  // Let's declare the variables we will use
  double w, Z, U;   // distinguish w from W. W is a vector of w
  
  // Start the loop
  for(int i(0);i<size;++i){
    Step1: Z = (rbeta(1,d/2.,d/2.))(0);
    w = (1.-(1.+b)*Z)/(1.-(1.-b)*Z);
    U = (runif(1,0.,1.))(0);
    
    //step 2ra
    if(lambda*w+d*log(1-x*w)-c < log(U)){goto Step1;}
    W(i)=w;
  }
}


// The following function complete the algorithm by performing the step 4 and
// the rotation toward the mean directional.
// It needs the sample size and theta as intensity parameter x mean directional

// [[Rcpp::export]]
SEXP cpprvMF(const int& size,const arma::vec& theta){
  int p=theta.n_rows;            //hypersphere dimension
  double lambda=norm(theta);     //intensity parameter
  arma::mat X;                         //Output matrix
  
  // if lambda=0 sample uniform; that is normal/norm
  if(lambda==0){
    NumericMatrix Xtemp(size,p,rnorm(size*p).begin());
    X=arma::normalise(as<arma::mat>(Xtemp),2,1);  //nomalize rows by their norm   
  }
  else{
    double d=p-1;
    // Compute W
    arma::vec W(size);        //Void W
    rw(size, lambda, d, W);   //Fill W using rw
    arma::mat Wplus=repmat(W,1,d);  //Reshape to [W W W ... W] of dimension (n,d)
    //mean direction parameter
    arma::vec mu=theta/lambda;           
    // Necessary variables declaration
    NumericMatrix Vtemp(size,d,rnorm(size*d).begin());
    arma::mat V=arma::normalise(as<arma::mat>(Vtemp),2,1);
    arma::mat X1=sqrt(1-Wplus % Wplus) % V;
    X=arma::join_rows(X1,W);
    //Rotation
    // To get samples from a vMF distribution with arbitrary mean direction
    // parameter µ, X is multiplied from the right with a matrix where the
    // first (m − 1) columns consist of unitary basis vectors of
    // the subspace orthogonal to µ and the last column is equal to µ. See
    // Hornik, K., & Grün, B. (2014). movMF: An R package for fitting mixtures
    // of von Mises-Fisher distributions. Journal of Statistical Software, 
    // 58(10), 1-31.
    arma::mat Q,R;
    arma::qr(Q, R, mu);    //QR decomposition to get subsâce orthogonal to µ
    IntegerVector seqcol=seq_len(d);
    Q=Q.cols(as<arma::uvec>(seqcol));
    Q=join_rows(Q,mu);
    X=X*Q.t();
  }
  return wrap(X);
}


// The normalization constant which depends on modified Bessel function
// in von Mises-Fisher distribution

// [[Rcpp::export]]

double cppCpvMF(const int& p, const double& k){
  if(k==0){ /*If k=0 return 1*/  return 1;}
  NumericVector temp=gamma(NumericVector::create(p/2.));
  return pow(k/2.,(p/2.)-1.)/(temp(0)*boost::math::cyl_bessel_i((p/2.)-1.,k));
}

// Compute the von Mises-Fisher density
// [[Rcpp::export]]

SEXP cppdvMF(arma::mat& z, arma::vec& theta){
  int Nrow=theta.n_rows, Ncol=theta.n_cols;
  if(Ncol>Nrow){
    theta=theta.t();                //to be sure that theta is a column
  }
  return wrap(cppCpvMF(z.n_cols,arma::norm(theta))*exp(z*theta));
}
