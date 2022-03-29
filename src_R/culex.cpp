#include "culex.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

using culex_stochastic = culex<int>;


// [[Rcpp::export]]
Rcpp::XPtr<culex_stochastic> create_culex_stochastic(const int p, const std::vector<int>& tau_E, const std::vector<int>& tau_L, const std::vector<int>& tau_P, const double dt) {
  return Rcpp::XPtr<culex_stochastic>(
    new culex<int>(p, tau_E, tau_L, tau_P, dt),
    true
  );
};

// [[Rcpp::export]]
void step_culex_stochastic(Rcpp::XPtr<culex_stochastic> mod, const Rcpp::List& parameters) {
  mod->update(parameters);
}

// [[Rcpp::export]]
void set_A_stochastic(Rcpp::XPtr<culex_stochastic> mod, arma::Row<int> A) {
  mod->A = A;
};

// [[Rcpp::export]]
arma::Row<int> get_A_stochastic(Rcpp::XPtr<culex_stochastic> mod) {
  return mod->A;
};

// [[Rcpp::export]]
arma::Row<int> get_E_stochastic(Rcpp::XPtr<culex_stochastic> mod) {
  return arma::sum(mod->E, 0);
};

// [[Rcpp::export]]
arma::Row<int> get_L_stochastic(Rcpp::XPtr<culex_stochastic> mod) {
  return arma::sum(mod->L, 0);
};

// [[Rcpp::export]]
arma::Row<int> get_P_stochastic(Rcpp::XPtr<culex_stochastic> mod) {
  return arma::sum(mod->P, 0);
};