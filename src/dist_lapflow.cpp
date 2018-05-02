#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double lfdistance(arma::mat& L1, arma::mat& L2, double inct){
  // 1. parameters
  double distvalue = 0.0;
  double incvalue  = 100.0;
  const double thr = 0.01;
  const int N = L1.n_cols;

  // 2. eigendecomposition
  arma::vec val1(N,fill::zeros);
  arma::vec val2(N,fill::zeros);
  arma::mat mat1(N,N,fill::zeros);
  arma::mat mat2(N,N,fill::zeros);
  arma::eig_sym(val1,mat1,L1);
  arma::eig_sym(val2,mat2,L2);

  arma::mat c1old(N,N,fill::eye);
  arma::mat c2old(N,N,fill::eye);
  arma::mat c1new(N,N,fill::zeros);
  arma::mat c2new(N,N,fill::zeros);
  arma::mat dmatold(N,N,fill::zeros); // difference between c1old and c2old
  arma::mat dmatnew(N,N,fill::zeros); // difference between c1new and c2new
  arma::mat dmatprocess(N,N,fill::zeros); // abs(dmatold-dmatnew) with zeroed diagonals

  // 3. use while loop
  double timevalue = inct;
  while (incvalue > thr){
    // 3-1. update c1new, c2new, and dmatnew
    c1new = mat1*arma::diagmat(exp(-timevalue*val1))*mat1.t();
    c2new = mat2*arma::diagmat(exp(-timevalue*val2))*mat2.t();
    dmatnew = c1new-c2new;

    // 3-2. dmatprocess
    dmatprocess = arma::abs(dmatold-dmatnew);
    dmatprocess.diag().zeros(); // fill with zeros

    // 3-3. update the distvalue
    distvalue = distvalue + arma::accu(dmatprocess);

    // 3-4. update variables
    incvalue = dmatprocess.max(); // maximum discrepancy
    c1old = c1new;
    c2old = c2new;
    dmatold = dmatnew;
    timevalue += inct;

    // 3-5. conditional stopping : timevalue
    if (timevalue > 100){
      break;
    }
  }
  return(distvalue);
}
