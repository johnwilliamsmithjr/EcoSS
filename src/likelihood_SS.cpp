#include <Rcpp.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include "lmvnd_ll.h"
// load header files and libraries
//#include "lgd.h"
using namespace std;

// [[Rcpp::export]]
RcppExport SEXP SSLL(SEXP Cobs_, SEXP C_, SEXP Sigma_, SEXP Cpred_, SEXP sd_, SEXP init_mean_, SEXP init_sd_, int N, SEXP G_, SEXP G_e_, SEXP g_sig_,
                     SEXP g_tau_, SEXP Gpred_, std::vector<int> ind) {
  // Inputs: Matrix Cobs_ of carbon stock observations
  // Matrix C_ of latent carbon stocks
  // Vector Sigma_ of variances for observational error
  // Vector sd_ of process variances
  double log_lik = 0;
  // initialize log likelihood
  Rcpp::NumericMatrix Cobs = Rcpp::as<Rcpp::NumericMatrix >(Cobs_);
  Rcpp::NumericMatrix C = Rcpp::as<Rcpp::NumericMatrix >(C_);
  Rcpp::NumericVector Sigma = Rcpp::as<Rcpp::NumericVector >(Sigma_);
  Rcpp::NumericVector sd = Rcpp::as<Rcpp::NumericVector >(sd_);
  Rcpp::NumericVector init_mean = Rcpp::as<Rcpp::NumericVector>(init_mean_);
  Rcpp::NumericVector init_sd = Rcpp::as<Rcpp::NumericVector>(init_sd_);
  Rcpp::NumericMatrix Cpred = Rcpp::as<Rcpp::NumericMatrix >(Cpred_);
  Rcpp::NumericVector Cobsrow = Cobs(0, Rcpp::_ );
  Rcpp::NumericVector Crow = C(0,Rcpp::_);
  Rcpp::NumericVector Cprow = Cpred(0,Rcpp::_);
  Rcpp::NumericVector G = Rcpp::as<Rcpp::NumericVector>(G_);
  Rcpp::NumericVector G_e = Rcpp::as<Rcpp::NumericVector>(G_e_);
  Rcpp::NumericVector g_sig = Rcpp::as<Rcpp::NumericVector>(g_sig_);
  Rcpp::NumericVector g_tau = Rcpp::as<Rcpp::NumericVector>(g_tau_);
  Rcpp::NumericVector G_pred = Rcpp::as<Rcpp::NumericVector>(Gpred_);
  // coerce vector and matrix valued inputs to Rcpp objects
  log_lik = LMVND_LL::lmvnd(Rcpp::as<vector<double> >(Sigma), Rcpp::as<vector<double> > (Cobsrow), Rcpp::as<vector<double> > (Crow), 5);
  log_lik = log_lik + LMVND_LL::lmvnd(Rcpp::as<vector<double> >(init_sd), Rcpp::as<vector<double> >(Crow), Rcpp::as<vector<double> >(init_mean), 5);
  for (int i=1; i<N; i++) {
    Cobsrow = Cobs(i,Rcpp::_);
    Crow = C(i,Rcpp::_);
    Cprow = Cpred(i, Rcpp::_);
    log_lik = log_lik + LMVND_LL::lmvnd(Rcpp::as<vector<double> >(sd), Rcpp::as<vector<double> > (Crow), Rcpp::as<vector<double> > (Cprow), 5);
    if (std::count(ind.begin(), ind.end(), i)){
      Cobsrow = Cobs(i,Rcpp::_);
      log_lik = log_lik + LMVND_LL::lmvnd(Rcpp::as<vector<double> >(Sigma), Rcpp::as<vector<double> > (Cobsrow), Rcpp::as<vector<double> > (Crow), 5);
    }
  }
  log_lik = log_lik + LMVND_LL::lmvnd(Rcpp::as<vector<double> >(g_sig), Rcpp::as<vector<double> > (G_pred), Rcpp::as<vector<double> > (G_e), N);
  //log_lik = log_lik + LMVND_LL::lmvnd(Rcpp::as<vector<double> >(g_tau), Rcpp::as<vector<double> > (G), Rcpp::as<vector<double> > (G_pred), N);
  //double a = Rcpp::as<double>(a_);
  //double b = Rcpp::as<double>(b_);
  //vector<double> X = Rcpp::as<vector<double> >(sd_);
  
  //for (int j=0; j<5; j++) {
  //  double lgx = X[j];
  //  log_lik = log_lik + LMVND_LL::lgammad(lgx, a, b);
  //}
  return Rcpp::wrap(log_lik);
}