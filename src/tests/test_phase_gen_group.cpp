#include <Rcpp.h>
using namespace Rcpp;

#include "../phase_gen_estep_group.h"
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
void test_estep_group(NumericVector alpha,
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
  auto model = GPH<NumericVector, NumericMatrix, IntegerVector>(alpha, Q, P, xi, qv, di);
  
  auto work = GPHWorkSpaceGroup(m, n);
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
  llfv = llf(model, dat, options.poisson_eps);
  Rcout << "LLF=" << llfv << std::endl;
}

// [[Rcpp::export]]
void test_estep_group_csc(NumericVector alpha,
                      S4 Q0,
                      NumericVector xi,
                      S4 P0,
                      S4 H0,
                      List data) {
  int n = alpha.length();
  int m = 5;
  
  using SparseT = S4matrix<CSCMatrixT>;
  auto Q = SparseT(Q0);
  auto P = SparseT(P0);
  auto H = SparseT(H0);
  
  auto tdat = as<NumericVector>(data["time"]);
  auto gdat = as<IntegerVector>(data["counts"]);
  auto idat = as<IntegerVector>(data["indicators"]);
  int glast = as<int>(data["last"]);
  double maxtime = as<double>(data["maxtime"]);
  auto eres = GPHEres<NumericVector, SparseT>(
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
  auto model = GPH<NumericVector, SparseT, IntegerVector>(alpha, Q, P, xi, qv, di);
  
  auto work = GPHWorkSpaceGroup(m, n);
  llfv = estep(model, dat, eres, options, work);
  Rcout << llfv << std::endl;
  Rcout << eres.eb << std::endl;
  Rcout << sum(eres.eb) << std::endl;
  Rcout << eres.ey << std::endl;
  Rcout << sum(eres.ey) << std::endl;
  Rcout << eres.ez << std::endl;
  Rcout << sum(eres.ez) << std::endl;
  // Rcout << eres.en << std::endl;
  mstep(eres, model, options);
  llfv = estep(model, dat, eres, options, work);
  Rcout << llfv << std::endl;
}

// [[Rcpp::export]]
void test_estep_group_csr(NumericVector alpha,
                          S4 Q0,
                          NumericVector xi,
                          S4 P0,
                          S4 H0,
                          List data) {
  int n = alpha.length();
  int m = 5;
  
  using SparseT = S4matrix<CSRMatrixT>;
  auto Q = SparseT(Q0);
  auto P = SparseT(P0);
  auto H = SparseT(H0);
  
  auto tdat = as<NumericVector>(data["time"]);
  auto gdat = as<IntegerVector>(data["counts"]);
  auto idat = as<IntegerVector>(data["indicators"]);
  int glast = as<int>(data["last"]);
  double maxtime = as<double>(data["maxtime"]);
  auto eres = GPHEres<NumericVector, SparseT>(
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
  auto model = GPH<NumericVector, SparseT, IntegerVector>(alpha, Q, P, xi, qv, di);
  
  auto work = GPHWorkSpaceGroup(m, n);
  llfv = estep(model, dat, eres, options, work);
  Rcout << llfv << std::endl;
  Rcout << eres.eb << std::endl;
  Rcout << sum(eres.eb) << std::endl;
  Rcout << eres.ey << std::endl;
  Rcout << sum(eres.ey) << std::endl;
  Rcout << eres.ez << std::endl;
  Rcout << sum(eres.ez) << std::endl;
  // Rcout << eres.en << std::endl;
  mstep(eres, model, options);
  llfv = estep(model, dat, eres, options, work);
  Rcout << llfv << std::endl;
}

// [[Rcpp::export]]
void test_estep_group_coo(NumericVector alpha,
                          S4 Q0,
                          NumericVector xi,
                          S4 P0,
                          S4 H0,
                          List data) {
  int n = alpha.length();
  int m = 5;
  
  using SparseT = S4matrix<COOMatrixT>;
  auto Q = SparseT(Q0);
  auto P = SparseT(P0);
  auto H = SparseT(H0);
  
  auto tdat = as<NumericVector>(data["time"]);
  auto gdat = as<IntegerVector>(data["counts"]);
  auto idat = as<IntegerVector>(data["indicators"]);
  int glast = as<int>(data["last"]);
  double maxtime = as<double>(data["maxtime"]);
  auto eres = GPHEres<NumericVector, SparseT>(
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
  auto model = GPH<NumericVector, SparseT, IntegerVector>(alpha, Q, P, xi, qv, di);
  
  auto work = GPHWorkSpaceGroup(m, n);
  llfv = estep(model, dat, eres, options, work);
  Rcout << llfv << std::endl;
  Rcout << eres.eb << std::endl;
  Rcout << sum(eres.eb) << std::endl;
  Rcout << eres.ey << std::endl;
  Rcout << sum(eres.ey) << std::endl;
  Rcout << eres.ez << std::endl;
  Rcout << sum(eres.ez) << std::endl;
  // Rcout << eres.en << std::endl;
  mstep(eres, model, options);
  llfv = estep(model, dat, eres, options, work);
  Rcout << llfv << std::endl;
}

/*** R
alpha <- c(0.2, 0.6, 0.2)
Q <- rbind(
  c(-2.0, 1.0, 0.0),
  c(2.0, -5.0, 1.0),
  c(3.0, 2.0, -8.0))
P <- matrix(0, 3, 3)
H <- matrix(0, 3, 3)
xi <- c(1.0, 2.0, 3.0)
dat <- list(time=c(1,2,1,3,4), counts=c(1,3,-1,2,4), indicators=c(0,0,0,1,0), last=10, maxtime=4)
test_estep_group(alpha, Q, xi, P, H, dat)
print(alpha)
print(Q)
print(xi)

alpha <- c(0.2, 0.6, 0.2)
Q0 <- rbind(
  c(-2.0, 1.0, 0.0),
  c(2.0, -5.0, 1.0),
  c(3.0, 2.0, -8.0))
Q <- as(Q0, "dgCMatrix")
P <- as(Q0, "dgCMatrix")
H <- as(Q0, "dgCMatrix")
xi <- c(1.0, 2.0, 3.0)
dat <- list(time=c(1,2,1,3,4), counts=c(1,3,-1,2,4), indicators=c(0,0,0,1,0), last=10, maxtime=4)
test_estep_group_csc(alpha, Q, xi, P, H, dat)
print(alpha)
print(Q)
print(xi)

alpha <- c(0.2, 0.6, 0.2)
Q0 <- rbind(
  c(-2.0, 1.0, 0.0),
  c(2.0, -5.0, 1.0),
  c(3.0, 2.0, -8.0))
Q <- as(Q0, "dgRMatrix")
P <- as(Q0, "dgRMatrix")
H <- as(Q0, "dgRMatrix")
xi <- c(1.0, 2.0, 3.0)
dat <- list(time=c(1,2,1,3,4), counts=c(1,3,-1,2,4), indicators=c(0,0,0,1,0), last=10, maxtime=4)
test_estep_group_csr(alpha, Q, xi, P, H, dat)
print(alpha)
print(Q)
print(xi)
  
alpha <- c(0.2, 0.6, 0.2)
Q0 <- rbind(
  c(-2.0, 1.0, 0.0),
  c(2.0, -5.0, 1.0),
  c(3.0, 2.0, -8.0))
Q <- as(Q0, "dgTMatrix")
P <- as(Q0, "dgTMatrix")
H <- as(Q0, "dgTMatrix")
xi <- c(1.0, 2.0, 3.0)
dat <- list(time=c(1,2,1,3,4), counts=c(1,3,-1,2,4), indicators=c(0,0,0,1,0), last=10, maxtime=4)
test_estep_group_coo(alpha, Q, xi, P, H, dat)
print(alpha)
print(Q)
print(xi)
*/
