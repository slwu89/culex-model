/*
 * culex.hpp
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
  
  if (grate < gon_min) {
    grate = gon_min;
  }
  
  return 1.0 / grate;
}

// per-capita oviposition rate
// d: diapause
// G: duration of gonotrophic cycle
// pars: parameters
inline double oviposition(const double d, const double G, const Rcpp::List& pars) {
  
  double max_egg = pars["max_egg"];
  
  double egg_raft = d * max_egg * 0.5;
  double ovi = egg_raft / G;
  
  return ovi;
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
  
  void update(const Rcpp::List& parameters);
  
};

// constructor
template <typename T>
inline culex<T>::culex(const int p_, const std::vector<int>& tau_E_, const std::vector<int>& tau_L_, const std::vector<int>& tau_P_, const double dt_) :
  step(0), p(p_), dt(dt_), tau_E(tau_E_), tau_L(tau_L_), tau_P(tau_P_)
{
  int maxE = *max_element(tau_E.begin(), tau_E.end());
  int maxL = *max_element(tau_L.begin(), tau_L.end());
  int maxP = *max_element(tau_P.begin(), tau_P.end());
  
  // state
  E = arma::Mat<T>(maxE, p, arma::fill::zeros);
  L = arma::Mat<T>(maxL, p, arma::fill::zeros);
  P = arma::Mat<T>(maxP, p, arma::fill::zeros);
  A = arma::Row<T>(p, arma::fill::zeros);
  
  // shift matrices (multiply on right)
  arma::umat fillE(2, maxE);
  for (auto i = 0u; i < (maxE - 1); ++i) {
    fillE(0, i) = i;
    fillE(1, i) = i+1;
  }
  arma::Col<int> fillEvals(maxE, arma::fill::ones);
  shiftE = arma::SpMat<int>(fillE, fillEvals, maxE, maxE);
  
  arma::umat fillL(2, maxL);
  for (auto i = 0u; i < (maxL - 1); ++i) {
    fillL(0, i) = i;
    fillL(1, i) = i+1;
  }
  arma::Col<int> fillLvals(maxL, arma::fill::ones);
  shiftL = arma::SpMat<int>(fillL, fillLvals, maxL, maxL);
  
  arma::umat fillP(2, maxP);
  for (auto i = 0u; i < (maxP - 1); ++i) {
    fillP(0, i) = i;
    fillP(1, i) = i+1;
  }
  arma::Col<int> fillPvals(maxP, arma::fill::ones);
  shiftP = arma::SpMat<int>(fillP, fillPvals, maxP, maxP);
  
};

// stochastic update
template <>
inline void culex<int>::update(const Rcpp::List& parameters) {
  
  double tnow = this->step * this->dt;
  int tau_E = this->tau_E[this->step];
  int tau_L = this->tau_L[this->step];
  int tau_P = this->tau_P[this->step];
  
  double p0 = parameters["p0"];
  double p1 = parameters["p1"];
  
  // temperature
  double temp = temperature(tnow, parameters);
  
  // photoperiod
  double pp = photoperiod(tnow, parameters);
  double pp_1 = photoperiod(tnow - 1.0, parameters);
  
  // gonotrophic cycle
  double gon = gonotrophic(temp, parameters);
  
  // mortality
  double death_egg = death_egg_rate(temp, parameters);
  double death_larvae = death_larvae_rate(temp, parameters);
  double death_pupae = death_pupae_rate(temp, parameters);
  double death_adult = death_adult_rate(temp, parameters);
  
  arma::Row<int> larvae_tot = arma::sum(this->L, 0);
  std::vector<double> death_larvae_tot(this->p, death_larvae);
  for (auto i = 0u; i < this->p; ++i) {
    death_larvae_tot[i] += p0 * larvae_tot(i) / (p1 + larvae_tot(i));
  }
  
  // diapause and egg laying
  double dia;
  if (pp > pp_1) {
    dia = diapause_spring(pp);
  } else {
    dia = diapause_autumn(pp);
  }
  
  // egg laying
  arma::Row<int> lambda(this->p, arma::fill::zeros);
  int i{0};
  lambda.for_each([&i, this, dia, gon, &parameters](arma::Row<int>::elem_type& val) {
    double lambda_mean = oviposition(dia, gon, parameters) * this->A(i) * this->dt;
    val = R::rpois(lambda_mean);
    i++;
  });
  
  // survival 
  this->E.for_each([death_egg = death_egg, dt = this->dt](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, R::pexp(death_egg * dt, 1.0, 0, 0));  
    }
  });
  
  i = 0;
  this->L.for_each([&death_larvae_tot, &i, dt = this->dt](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, R::pexp(death_larvae_tot[i] * dt, 1.0, 0, 0));  
    }
    i++;
  });
  
  this->P.for_each([death_pupae = death_pupae, dt = this->dt](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, R::pexp(death_pupae * dt, 1.0, 0, 0));  
    }
  });
  
  this->A.for_each([death_adult = death_adult, dt = this->dt](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, R::pexp(death_adult * dt, 1.0, 0, 0));  
    }
  });
  
  // advancement
  arma::Row<int> E2L = this->E.row(0);
  this->E.row(0).zeros();
  this->E = this->shiftE * this->E;
  this->E.row(tau_E-1) = lambda;
  
  arma::Row<int> L2P = this->L.row(0);
  this->L.row(0).zeros();
  this->L = this->shiftL * this->L;
  this->L.row(tau_L-1) = E2L;
  
  arma::Row<int> P2A = this->P.row(0);
  this->P.row(0).zeros();
  this->P = this->shiftP * this->P;
  this->P.row(tau_P-1) = L2P;;
  
  this->A += P2A;
  
  this->step += 1;
  
};

// deterministic update

#endif