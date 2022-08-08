#include <Rcpp.h>
using namespace Rcpp;

#include "../phase_gen_estep_grouppoi.h"
#include "../phase_gen_mstep.h"
#include "../emfit.h"

// #include "gperftools/profiler.h"
// 
// // [[Rcpp::export]]
// SEXP start_profiler(SEXP str) {
//   ProfilerStart(as<const char*>(str));
//   return R_NilValue;
// }
// 
// // [[Rcpp::export]]
// SEXP stop_profiler() {
//   ProfilerStop();
//   return R_NilValue;
// }

// [[Rcpp::export]]
void test_estep_group_poi(
    double omega,
    NumericVector alpha,
                      NumericMatrix Q,
                      NumericVector xi,
                      NumericMatrix P,
                      NumericMatrix H,
                      List data) {
  int n = alpha.length();
  int m = 5;
  
  auto tdat = as<NumericVector>(data["time"]);
  auto gdat = as<IntegerVector>(data["counts"]);
  auto idat = as<IntegerVector>(data["indicators"]);
  int glast = as<int>(data["last"]);
  double maxtime = as<double>(data["maxtime"]);
  auto eres = GPHEres<NumericVector, NumericMatrix>(
    NumericVector(n),
    NumericVector(n),
    NumericVector(n),
    H);
  auto options = EMOptions();
  options.poisson_eps = 1.0e-8;
  options.ufactor = 1.01;
  double llfv;
  auto dat = PHGroupSample<NumericVector,IntegerVector,IntegerVector>(tdat, gdat, idat, maxtime, glast);
  
  IntegerVector di(n);
  diag(Q, di);
  copy(Q, P);
  double qv = unif(P, di, 1.01);
  auto gph = GPH<NumericVector, NumericMatrix, IntegerVector>(alpha, Q, P, xi, qv, di);
  auto model = GPHPoi<GPH<NumericVector, NumericMatrix, IntegerVector>>(gph, omega);
  
  auto work = GPHWorkSpaceGroupPoi(m, n);
  llfv = estep(model, dat, eres, options, work);
  Rcout << llfv << std::endl;
  Rcout << eres.eb << std::endl;
  Rcout << sum(eres.eb) << std::endl;
  Rcout << eres.ey << std::endl;
  Rcout << sum(eres.ey) << std::endl;
  Rcout << eres.ez << std::endl;
  Rcout << sum(eres.ez) << std::endl;
  Rcout << eres.en << std::endl;
  mstep(eres, model, options);
  llfv = estep(model, dat, eres, options, work);
  Rcout << llfv << std::endl;
  Rcout << "LLF=" << llf(model, dat, options.poisson_eps) << std::endl;
}

/*** R
alpha <- c(0.2, 0.6, 0.2)
Q <- rbind(
  c(-2.0, 1.0, 0.0),
  c(2.0, -5.0, 1.0),
  c(3.0, 2.0, -8.0))
omega <- 10
P <- matrix(0, 3, 3)
H <- matrix(0, 3, 3)
xi <- c(1.0, 2.0, 3.0)
dat <- list(time=c(1,2,1,3,4), counts=c(1,3,-1,2,4), indicators=c(0,0,0,1,0), last=10, maxtime=4)
test_estep_group_poi(omega, alpha, Q, xi, P, H, dat)
print(alpha)
print(Q)
print(xi)
*/
