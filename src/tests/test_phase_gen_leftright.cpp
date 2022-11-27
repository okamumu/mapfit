#include <Rcpp.h>
using namespace Rcpp;

#include "../phase_gen_estep_leftright.h"
#include "../phase_gen_mstep.h"
#include "../emfit.h"

// [[Rcpp::export]]
void test_estep_leftright(NumericVector alpha,
                      NumericMatrix Q,
                      NumericVector xi,
                      NumericMatrix P,
                      NumericMatrix H,
                      List data) {
  int n = alpha.length();
  int m = as<NumericVector>(data["time"]).length();
  
  auto tdat = as<NumericVector>(data["time"]);
  auto nu = as<IntegerVector>(data["nu"]);
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
  auto dat = PHLeftRightSample<NumericVector,IntegerVector>(tdat, nu, maxtime);
  
  IntegerVector di(n);
  diag(Q, di);
  copy(Q, P);
  double qv = unif(P, di, 1.01);
  auto model = GPH<NumericVector, NumericMatrix, IntegerVector>(alpha, Q, P, xi, qv, di);
  
  auto work = GPHWorkSpaceLeftRight(m, n);
  llfv = estep(model, dat, eres, options, work);
  Rcout << "llf=" << llfv << std::endl;
  Rcout << eres.eb << std::endl;
  Rcout << sum(eres.eb) << std::endl;
  Rcout << eres.ey << std::endl;
  Rcout << sum(eres.ey) << std::endl;
  Rcout << eres.ez << std::endl;
  Rcout << sum(eres.ez) << std::endl;
  Rcout << eres.en << std::endl;
  Rcout << "eres.etotal=" << eres.etotal << std::endl;
  mstep(eres, model, options);
  Rcout << "alpha=" << model.alpha << std::endl;
  Rcout << "Q=" << model.Q << std::endl;
  Rcout << "P=" << model.P << std::endl;
  Rcout << "xi=" << model.xi << std::endl;
  Rcout << "qv=" << model.qv << std::endl;
  llfv = estep(model, dat, eres, options, work);
  Rcout << "llf=" << llfv << std::endl;
  mstep(eres, model, options);
  llfv = estep(model, dat, eres, options, work);
  Rcout << "llfv=" << llfv << std::endl;
  mstep(eres, model, options);
  llfv = estep(model, dat, eres, options, work);
  Rcout << "llfv=" << llfv << std::endl;
  Rcout << "LLF=" << llf(model, dat, options.poisson_eps) << std::endl;
}

// [[Rcpp::export]]
void test_estep_leftright2(NumericVector alpha,
                          S4 Q0,
                          NumericVector xi,
                          S4 P0,
                          S4 H0,
                          List data) {
  int n = alpha.length();
  int m = as<NumericVector>(data["time"]).length();
  
  using MatrixT = S4matrix<CSCMatrixT>;
  auto Q = MatrixT(Q0);
  auto P = MatrixT(P0);
  auto H = MatrixT(H0);

  auto tdat = as<NumericVector>(data["time"]);
  auto nu = as<IntegerVector>(data["nu"]);
  double maxtime = as<double>(data["maxtime"]);
  auto eres = GPHEres<NumericVector, MatrixT>(
    NumericVector(n),
    NumericVector(n),
    NumericVector(n),
    H);
  auto options = EMOptions();
  options.poisson_eps = 1.0e-8;
  options.ufactor = 1.01;
  double llfv;
  auto dat = PHLeftRightSample<NumericVector,IntegerVector>(tdat, nu, maxtime);
  
  IntegerVector di(n);
  diag(Q, di);
  copy(Q, P);
  double qv = unif(P, di, 1.01);
  auto model = GPH<NumericVector, MatrixT, IntegerVector>(alpha, Q, P, xi, qv, di);
  
  auto work = GPHWorkSpaceLeftRight(m, n);
  for (int k=0; k<20; k++) {
    llfv = estep(model, dat, eres, options, work);
    mstep(eres, model, options);
    Rcout << "llf=" << llfv << std::endl;
  }
  // Rcout << eres.eb << std::endl;
  // Rcout << sum(eres.eb) << std::endl;
  // Rcout << eres.ey << std::endl;
  // Rcout << sum(eres.ey) << std::endl;
  // Rcout << eres.ez << std::endl;
  // Rcout << sum(eres.ez) << std::endl;
  //Rcout << eres.en << std::endl;
  // Rcout << "eres.etotal=" << eres.etotal << std::endl;
  // mstep(eres, model, options);
  // Rcout << "alpha=" << model.alpha << std::endl;
  // //Rcout << "Q=" << model.Q << std::endl;
  // //Rcout << "P=" << model.P << std::endl;
  // Rcout << "xi=" << model.xi << std::endl;
  // Rcout << "qv=" << model.qv << std::endl;
  // llfv = estep(model, dat, eres, options, work);
  // Rcout << "llf=" << llfv << std::endl;
  // mstep(eres, model, options);
  // llfv = estep(model, dat, eres, options, work);
  // Rcout << "llfv=" << llfv << std::endl;
  // mstep(eres, model, options);
  // llfv = estep(model, dat, eres, options, work);
  // Rcout << "llfv=" << llfv << std::endl;
  // Rcout << "LLF=" << llf(model, dat, options.poisson_eps) << std::endl;
}

/*** R

alpha <- c(0.2, 0.6, 0.2)
Q <- rbind(
  c(-0.2, 0.1, 0.0),
  c(0.2, -0.5, 0.1),
  c(0.3, 0.2, -0.8))
P <- matrix(0, 3, 3)
H <- matrix(0, 3, 3)
xi <- c(0.1, 0.2, 0.3)
dat <- list(time=c(0.0, 0.0, 1.0, 0.2, 0.8, 2.5), nu=c(3, 3, 0, 1, 3, 1), maxtime=2.5)
test_estep_leftright(alpha, Q, xi, P, H, dat)
print(alpha)
print(Q)
print(xi)
*/
