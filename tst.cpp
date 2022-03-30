#include <Rmath.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

// [[Rcpp::export]]
arma::Row<int> sample_vec(const arma::Row<int>& x, const arma::Mat<double>& psi) {
  arma::Row<int> out(x.size(), arma::fill::zeros);
  arma::Row<int> temp(x.size(), arma::fill::zeros);
  arma::Mat<double> psit = psi.t();
  for (auto i = 0u; i < x.size(); ++i) {
    R::rmultinom(x(i), psit.colptr(i), psi.n_cols, temp.memptr());
    out += temp;
  }
  return out;
};