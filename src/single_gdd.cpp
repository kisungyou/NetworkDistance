#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List cpp_gdd(arma::cube& vecs, arma::mat& vals, arma::vec timestamps){
  // 1. parameters
  const int p  = vals.n_rows;
  const int N  = vals.n_cols;
  const int nT = timestamps.n_elem;

  // 2. prepare for recording
  arma::mat dist_mat(N,N,fill::zeros);
  arma::mat time_mat(N,N,fill::zeros);

  // 3. nesting the loops
  arma::mat L1(p,p,fill::zeros);
  arma::mat L2(p,p,fill::zeros);
  arma::vec D1(p,fill::zeros);
  arma::vec D2(p,fill::zeros);
  arma::vec tvalues(nT,fill::zeros);
  for (int i=0;i<(N-1);i++){
    L1 = vecs.slice(i);
    D1 = vals.col(i);
    for (int j=(i+1);j<N;j++){
      L2 = vecs.slice(j);
      D2 = vals.col(j);

      for (int k=0;k<nT;k++){
        double ctvalue = timestamps(k);
        tvalues(k) = arma::norm((L1*arma::diagmat(exp(-ctvalue*D1))*L1.t())-(L2*arma::diagmat(exp(-ctvalue*D2))*L2.t()), "fro");
      }

      double distmax   = tvalues.max();
      double timemax   = timestamps(tvalues.index_max());
      dist_mat(i,j) = distmax;
      dist_mat(j,i) = distmax;
      time_mat(i,j) = timemax;
      time_mat(j,i) = timemax;
    }
  }

  // 4. return output
  List output;
  output["distmat"] = dist_mat;
  output["timemat"] = time_mat;
  return(output);
}
