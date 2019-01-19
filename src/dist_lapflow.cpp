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



//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double lfdistance_new(arma::mat& L1, arma::vec D1, arma::mat& L2, arma::vec D2, arma::vec timestamp){
  // 1. parameters
  const int nT = timestamp.n_elem;
  const int N  = L1.n_cols;

  // 2. prepare
  // 2-1. c1old : depending on the initial value of timestamp
  double t0 = timestamp(0);
  arma::mat c1old(N,N,fill::zeros);
  arma::mat c2old(N,N,fill::zeros);
  c1old = L1*arma::diagmat(exp(-t0*D1))*L1.t();
  c2old = L2*arma::diagmat(exp(-t0*D2))*L2.t();
  arma::mat c1new(N,N,fill::zeros);
  arma::mat c2new(N,N,fill::zeros);
  arma::mat dmatold = c1old-c2old;; // difference between c1old and c2old
  arma::mat dmatnew(N,N,fill::zeros); // difference between c1new and c2new
  arma::mat dmatprocess(N,N,fill::zeros); // abs(dmatold-dmatnew) with zeroed diagonals

  // 3. iterate
  double distvalue = 0.0;
  for (int i=1;i<nT;i++){
    // 3-1. current timestamp and new "c" matrices
    double tnow = timestamp(i);
    c1new = L1*arma::diagmat(exp(-tnow*D1))*L1.t();
    c2new = L2*arma::diagmat(exp(-tnow*D2))*L2.t();

    // 3-2. incremental change
    dmatnew = c1new-c2new;

    // 3-3. final update for incremental change
    dmatprocess = arma::abs(dmatold-dmatnew);
    dmatprocess.diag().zeros(); // fill with zeros

    // 3-4. update distance values
    distvalue += arma::accu(dmatprocess);

    // 3-5. update old&new ones
    c1old = c1new;
    c2old = c2new;
    dmatold = dmatnew;
  }
  return(distvalue);
}



//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat lfdistance_new_faster(arma::cube& vecs, arma::mat& vals, arma::vec timestamps){
  // 1. parameters
  const int N = vecs.n_slices; // number of networks
  const int p = vecs.n_cols;   // dimensionality

  arma::mat L1(p,p,fill::zeros);
  arma::mat L2(p,p,fill::zeros);
  arma::vec D1(p,fill::zeros);
  arma::vec D2(p,fill::zeros);

  // 2. compute with lfdistance_new
  arma::mat distmat(N,N,fill::zeros);
  for (int i=0;i<(N-1);i++){
    L1 = vecs.slice(i);
    D1 = vals.col(i);
    for (int j=(i+1);j<N;j++){
      L2 = vecs.slice(j);
      D2 = vals.col(j);

      double dval = lfdistance_new(L1, D1, L2, D2, timestamps);
      distmat(i,j) = dval;
      distmat(j,i) = dval;
    }
  }
  return(distmat);
}
