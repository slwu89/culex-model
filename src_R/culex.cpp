#include "culex.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]


// ---------- stochastic model interface ----------
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


// ---------- deterministic model interface ----------
using culex_deterministic = culex<double>;

// [[Rcpp::export]]
Rcpp::XPtr<culex_deterministic> create_culex_deterministic(const int p, const std::vector<int>& tau_E, const std::vector<int>& tau_L, const std::vector<int>& tau_P, const double dt) {
  return Rcpp::XPtr<culex_deterministic>(
    new culex<double>(p, tau_E, tau_L, tau_P, dt),
    true
  );
};

// [[Rcpp::export]]
void step_culex_deterministic(Rcpp::XPtr<culex_deterministic> mod, const Rcpp::List& parameters) {
  mod->update(parameters);
}

// [[Rcpp::export]]
void set_A_deterministic(Rcpp::XPtr<culex_deterministic> mod, arma::Row<double> A) {
  mod->A = A;
};

// [[Rcpp::export]]
arma::Row<double> get_A_deterministic(Rcpp::XPtr<culex_deterministic> mod) {
  return mod->A;
};

// [[Rcpp::export]]
arma::Row<double> get_E_deterministic(Rcpp::XPtr<culex_deterministic> mod) {
  return arma::sum(mod->E, 0);
};

// [[Rcpp::export]]
arma::Row<double> get_L_deterministic(Rcpp::XPtr<culex_deterministic> mod) {
  return arma::sum(mod->L, 0);
};

// [[Rcpp::export]]
arma::Row<double> get_P_deterministic(Rcpp::XPtr<culex_deterministic> mod) {
  return arma::sum(mod->P, 0);
};

