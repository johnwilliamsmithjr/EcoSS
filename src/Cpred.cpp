#include <Rcpp.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
//#include <cstdio.h>
#include <math.h>
using namespace std;

// [[Rcpp::export]]
RcppExport SEXP Cpred(SEXP p_, SEXP C_, SEXP LAI_, SEXP max_t_, SEXP min_t_, SEXP Ca_, SEXP lat_e_,
                      SEXP yearday_, SEXP Nit_e_, SEXP rad_, int N, SEXP G_) {
  //vector<double> mean = Rcpp::as<vector<double> >(mean_);
  //vector<double> Sigma = Rcpp::as<vector<double> >(Sigma_);
  //vector<double> X = Rcpp::as<vector<double> >(X_);
  Rcpp::NumericMatrix C = Rcpp::as<Rcpp::NumericMatrix >(C_);
  vector<double> p = Rcpp::as<vector<double> >(p_);
  vector<double> G = Rcpp::as<vector<double> >(G_);
  vector<double> LAI = Rcpp::as<vector<double> >(LAI_);
  vector<double> max_t = Rcpp::as<vector<double> >(max_t_);
  vector<double> min_t = Rcpp::as<vector<double> >(min_t_);
  vector<double> Ca = Rcpp::as<vector<double> >(Ca_);
  vector<double> yearday = Rcpp::as<vector<double> >(yearday_);
  double Nit_e = Rcpp::as<double>(Nit_e_);
  double lat_e = Rcpp::as<double>(lat_e_);
  vector<double> rad = Rcpp::as<vector<double> >(rad_);
  Rcpp::NumericVector Crow_ = C(0, Rcpp::_);
  vector<double> Crow = Rcpp::as<vector<double> > (Crow_);
  vector<double> Cf_pred(N);
  Cf_pred[0] = Crow[0];
  vector<double> A_cf(N);
  A_cf[0] = (1-p[4]);
  vector<double> b_cf(N);
  b_cf[0] = G[0]*(1-p[1])*p[2];
  vector<double> Cw_pred(N);
  Cw_pred[0] = Crow[1];
  vector<double> A_cw(N);
  A_cw[0] = (1-p[5]);
  vector<double> b_cw(N);
  b_cw[0] = G[0]*(1-p[1])*(1-p[2])*(1-p[3]);
  vector<double> Cr_pred(N);
  Cr_pred[0] = Crow[2];
  vector<double> A_cr(N);
  A_cr[0] = (1-p[6]);
  vector<double> b_cr(N);
  b_cr[0] = G[0]*(1-p[1])*(1-p[2])*p[3];
  vector<double> Clit_pred(N);
  Clit_pred[0] = Crow[3];
  vector<double> A_clit(N);
  A_clit[0] = (1 - (.5*exp(p[9]*.5*(min_t[0] + max_t[0])))*(p[7] + p[0]));
  vector<double> b_clit(N);
  b_clit[0] = Crow[0]*p[4] + Crow[2]*p[6];
  vector<double> Csom_pred(N);
  Csom_pred[0] = Crow[4];
  vector<double> A_csom(N);
  A_csom[0] = (1 - p[8]*.5*exp(p[9]*.5*(min_t[0] + max_t[0])));
  vector<double> b_csom(N);
  b_csom[0] = p[5]*Crow[1] + .5*exp(p[9]*.5*(min_t[0] + max_t[0]))*p[0]*Crow[3];
  //double Cpred[730][5] = {
  //  {0,0,0,0,0}
  //};
  //for (int j=0; j<5; j++){
  //  std::printf ("%f", Cpred(0,j));
 // }
  //for (int i=0; i<N; i++){
  //  std::printf ("%f", rad[i]);
  //}
  double psid = -2;
  double rtot = 1;
  double trange;
  double gs;
  double pp;
  double qq = -204.6453;
  double ci;
  double e0;
  double dec;
  double mult;
  double dayl;
  double cps;
  //vector<double> a[10];
  double a0 = p[10];
  double a1 = 0.0156935; 
  double a2 = 4.22273;
  double a3 = 208.868;
  double a4 = 0.0453194; 
  double a5 = 0.37836;
  double a6 = 7.19298;
  double a7 = 0.011136;
  double a8 = 2.1001; 
  double a9 = 0.789798;
  double pi = std::atan(1)*4;
  vector<double> Gpred(N);
  Gpred[1] = 0;
  //G[0] = 0;
  for (int i = 1; i<N; i++){
  Crow_ = C(i-1, Rcpp::_);
  Crow = Rcpp::as<vector<double> > (Crow_);
  //std::printf("%f \n", Crow[0]);
  trange = .5*(max_t[i] - min_t[i]);
  gs = pow(fabs(psid),(0.789798)) / (0.37836*rtot + trange);
  pp = LAI[i]*Nit_e/gs*a0*exp(a7*(max_t[i]));
  ci = .5*(Ca[i] + qq - pp + pow(pow(Ca[i]+qq-pp,2) - 4*(Ca[i]*qq - pp*a2), .5 ));
  //std::printf("%f", ci);
  e0 = (a6*pow(LAI[i],2)) / (pow(LAI[i],2) + a8);
  dec = -23.4*cos((360*(yearday[i] + 10)/365)*pi/180)*pi/180;
  mult = tan(lat_e)*tan(dec);
  if (mult >= 1){
    dayl = 24;
  } else if(mult <= -1){
    dayl = 0;
  } else{
    dayl = 24*acos(-mult)/pi;
  }
  cps = e0*rad[i]*gs*(Ca[i] - ci) / (e0*rad[i] + gs*(Ca[i] - ci)) ;
  Gpred[i] = cps*(a1*dayl + a4);
  A_cf[i] = (1-p[4]);
  b_cf[i] = G[i]*(1-p[1])*p[2];
  Cf_pred[i] = A_cf[i]*Crow[0] + b_cf[i];
  
  Cw_pred[i] = (1-p[5])*Crow[1] + G[i]*(1-p[1])*(1-p[2])*(1-p[3]);
  A_cw[i] = (1-p[5]);
  b_cw[i] = G[i]*(1-p[1])*(1-p[2])*(1-p[3]);
  
  Cr_pred[i] = (1-p[6])*Crow[2] + G[i]*(1-p[1])*(1-p[2])*p[3];
  A_cr[i] = (1-p[6]);
  b_cr[i] = G[i]*(1-p[1])*(1-p[2])*p[3];
  
  Clit_pred[i] = Crow[0]*p[4] + Crow[2]*p[6] + (1 - (.5*exp(p[9]*.5*(min_t[i] + max_t[i])))*(p[7] + p[0]))*Crow[3];
  A_clit[i] = (1 - (.5*exp(p[9]*.5*(min_t[i] + max_t[i])))*(p[7] + p[0]));
  b_clit[i] = Crow[0]*p[4] + Crow[2]*p[6];
  
  Csom_pred[i] = p[5]*Crow[1] + .5*exp(p[9]*.5*(min_t[i] + max_t[i]))*p[0]*Crow[3] + (1 - p[8]*.5*exp(p[9]*.5*(min_t[i] + max_t[i])))*Crow[4];
  A_csom[i] = (1 - p[8]*.5*exp(p[9]*.5*(min_t[i] + max_t[i])));
  b_csom[i] = p[5]*Crow[1] + .5*exp(p[9]*.5*(min_t[i] + max_t[i]))*p[0]*Crow[3];
  // Fix the p stuff!!
  //std::printf("%f \n", Csom_pred[i]);
  }
  Rcpp::List ret;
  ret["Cfpred"] = Cf_pred;
  ret["Cwpred"] = Cw_pred;
  ret["Crpred"] = Cr_pred;
  ret["Clitpred"] = Clit_pred;
  ret["Csompred"] = Csom_pred;
  ret["G"] = Gpred;
  ret["A_cf"] = A_cf;
  ret["A_cw"] = A_cw;
  ret["A_cr"] = A_cr;
  ret["A_clit"] = A_clit;
  ret["A_csom"] = A_csom;
  ret["b_cf"] = b_cf;
  ret["b_cw"] = b_cw;
  ret["b_cr"] = b_cr;
  ret["b_clit"] = b_clit;
  ret["b_csom"] = b_csom;
  return Rcpp::wrap(ret);
}
