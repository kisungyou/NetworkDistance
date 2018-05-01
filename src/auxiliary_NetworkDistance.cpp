/*
 * 01. aux_FrobeniusDiff
 *
 *
 */
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double aux_FrobeniusDiff(arma::mat& A, arma::mat& B){
  double output = arma::norm(A-B,"fro");
  return(output);
}
