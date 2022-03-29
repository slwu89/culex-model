/*
 * culex.h
 * struct to define the model and temperature/photoperiod functions (from https://github.com/davewi13/Temperate-Mosquito-DDE/blob/master/Chapter%202%20DDE%20code.f90)
 * mathematical model by David Ewing
 * stochastic adaption by Sean L. Wu
 * March 2022
 */

#ifndef CULEX_HPP
#define CULEX_HPP

#include <Rmath.h>
#include <RcppArmadillo.h>

#include <algorithm> // for max_element
#include <vector>


// --------------------------------------------------------------------------------
//  external forcing (temp, photoperiod)
// --------------------------------------------------------------------------------

// temperature as modified cosine function
inline double temperature(const double t, const Rcpp::List& pars) {
  
  double phi = pars["phi"]; // PHASE
  double lambda = pars["lambda"]; // A
  double mu = pars["mu"]; // M
  double gamma = pars["gamma"]; // POWER
  
  double temp = 0.0;
  
  if (t < 0.0) {
    temp = (mu - lambda) + lambda * 2.0 * pow(0.5 * (1.0 + cos(2.0 * M_PI * (0.0 - phi) / 365.0)), gamma);
  } else {
    temp = (mu - lambda) + lambda * 2.0 * pow(0.5 * (1.0 + cos(2.0 * M_PI * (t - phi) / 365.0)), gamma);
  }
  
  return temp;
}

// photoperiod
inline double photoperiod(const double t, const Rcpp::List& pars) {
  
  double L = pars["L"]; // latitude (51 in thesis)
  
  // define photoperiod values
  double EPS = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (t - 3.5)))));
  double NUM = sin(0.8333 * M_PI/ 180.0) + (sin(L * M_PI / 180.0) * sin(EPS));
  double DEN = cos(L * M_PI / 180.0) * cos(EPS);
  double DAYLIGHT = 24.0 - (24.0 / M_PI) * acos(NUM / DEN);
  
  return DAYLIGHT;
}


// --------------------------------------------------------------------------------
//   diapause
// --------------------------------------------------------------------------------

// pp: photoperiod
inline double diapause_spring(const double pp){
  return 1.0 / (1.0 + exp(5.0 * (14.0 - pp)));
}

inline double diapause_autumn(const double pp) {
  return 1.0 / (1.0 + exp(5.0 * (13.0 - pp)));
}


// --------------------------------------------------------------------------------
//   mortality
// --------------------------------------------------------------------------------

// egg mortality
inline double death_egg_rate(const double temp, const Rcpp::List& pars) {
  
  double nu_0E = pars["nu_0E"]; // U3
  double nu_1E = pars["nu_1E"]; // U4
  double nu_2E = pars["nu_2E"]; // U5
  double death_max = pars["death_max"];
  
  // calculate egg death rate
  double egg_d = nu_0E * exp(pow((temp - nu_1E) / nu_2E, 2.0));
  
  if (egg_d > death_max) {
    egg_d = death_max;
  }
  
  return egg_d;
}

// larvae mortality
inline double death_larvae_rate(const double temp, const Rcpp::List& pars) {
  
  double nu_0L = pars["nu_0L"]; // U3
  double nu_1L = pars["nu_1L"]; // U4
  double nu_2L = pars["nu_2L"]; // U5
  double death_max = pars["death_max"];
  
  // calculate egg death rate
  double larvae_d = nu_0L * exp(pow((temp - nu_1L) / nu_2L, 2.0));
  
  if (larvae_d > death_max) {
    larvae_d = death_max;
  }
  
  return larvae_d;
}

// pupal mortality
inline double death_pupae_rate(const double temp, const Rcpp::List& pars) {
  
  double nu_0P = pars["nu_0P"]; // U3
  double nu_1P = pars["nu_1P"]; // U4
  double nu_2P = pars["nu_2P"]; // U5
  double death_max = pars["death_max"];
  
  // calculate egg death rate
  double pupal_d = nu_0P * exp(pow((temp - nu_1P)/nu_2P, 2.0));
  
  if (pupal_d > death_max) {
    pupal_d = death_max;
  }
  
  return pupal_d;
}

// adult mortality
inline double death_adult_rate(const double temp, const Rcpp::List& pars) {
  
  double alpha_A = pars["alpha_A"]; // ALPHA
  double beta_A = pars["beta_A"]; // # BETA
  double death_min_a = pars["death_min_a"];
  
  // calculate adult death rate
  double adult_d = alpha_A * pow(temp, beta_A);
  
  if(adult_d < death_min_a){
    adult_d = death_min_a;
  }
  
  return adult_d;
}


// --------------------------------------------------------------------------------
//   lifecycle stage progression rates
// --------------------------------------------------------------------------------

// G: duration of gonotrophic cycle
inline double gonotrophic(const double temp, const Rcpp::List& pars) {
  
  double q1 = pars["q1"]; // KG
  double q2 = pars["q2"]; // QG
  double q3 = pars["q3"]; // BG
  double gon_min = pars["gon_min"];
  
  // calculate gonotrophic cycle length
  double grate;
  if (temp < 0.0) {
    grate = 0.0333;
  } else {
    grate = q1 / (1.0 + q2*exp(-q3*temp));
  }
  
  if(grate < gon_min){
    grate = gon_min;
  }
  
  return 1.0 / grate;
}


// --------------------------------------------------------------------------------
//   model struct
// --------------------------------------------------------------------------------

// fields:
// E: eggs (steps X patches)
// L: larvae (steps X patches)
// P: pupae (steps X patches)
// A: adults (patches)
// shiftE: sparse matrix when multiplies E on the left to advance by 1 step
// shiftL: sparse matrix when multiplies L on the left to advance by 1 step
// shiftP: sparse matrix when multiplies P on the left to advance by 1 step
// step: current step
// p: number of patches
// dt: size of time step
// tau_E: tau_E[t] gives the number of steps required for a new E at step t to mature
// tau_L: tau_L[t] gives the number of steps required for a new L at step t to mature
// tau_P: tau_P[t] gives the number of steps required for a new P at step t to mature
template <typename T>
struct culex {
  
  arma::Mat<T> E;
  arma::Mat<T> L;
  arma::Mat<T> P;
  arma::Row<T> A;
  
  arma::SpMat<int> shiftE;
  arma::SpMat<int> shiftL;
  arma::SpMat<int> shiftP;
  
  int step;
  int p;
  double dt;
  
  std::vector<int> tau_E;
  std::vector<int> tau_L;
  std::vector<int> tau_P;
  
  culex(const int p_, const std::vector<int>& tau_E_, const std::vector<int>& tau_L_, const std::vector<int>& tau_P_, const double dt_);
  ~culex() = default;
  
};

#endif