#include <Rcpp.h>
using namespace Rcpp;

#include "phase_cf1.h"
#include "emfit.h"

// [[Rcpp::export]]
List emfit_cf1_wtime(NumericVector alpha,
                     NumericVector rate,
                     List data,
                     List options,
                     S4 Q0,
                     S4 P0,
                     S4 H0) {
  using MatrixT = S4matrix<CSCMatrixT>;
  auto Q = MatrixT(Q0);
  auto P = MatrixT(P0);
  auto H = MatrixT(H0);
  
  auto maxiter = as<int>(options["maxiter"]);
  auto atol = as<double>(options["abstol"]);
  auto rtol = as<double>(options["reltol"]);
  auto verbose = as<bool>(options["em.verbose"]);
  auto steps = as<int>(options["steps"]);
  auto ufactor = as<double>(options["uniform.factor"]);
  auto eps = as<double>(options["poisson.eps"]);
  
  int n = alpha.length();
  IntegerVector di(n);
  diag(Q, di);
  copy(Q, P);
  double qv = unif(P, di, ufactor);
  NumericVector al(n);
  copy(alpha, al);
  NumericVector xi(n);
  xi[n-1] = rate[n-1];
  auto gph = GPH<NumericVector, MatrixT, IntegerVector>(al, Q, P, xi, qv, di);
  auto model = CF1<NumericVector, GPH<NumericVector, MatrixT, IntegerVector>>(alpha, rate, gph);
  
  auto tdat = as<NumericVector>(data["time"]);
  auto wdat = as<NumericVector>(data["weights"]);
  double maxtime = as<double>(data["maxtime"]);
  auto m = tdat.length();
  auto dat = PHWeightSample<NumericVector,NumericVector>(tdat, wdat, maxtime);
  
  auto eres = GPHEres<std::vector<double>, MatrixT>(
    std::vector<double>(n),
    std::vector<double>(n),
    std::vector<double>(n),
    H);
  auto work = GPHWorkSpace1<std::vector<double>>(m);
  
  auto opts = EMOptions();
  opts.maxiter = maxiter;
  opts.atol = atol;
  opts.rtol = rtol;
  opts.steps = steps;
  opts.verbose = verbose;
  opts.ufactor = ufactor;
  opts.poisson_eps = eps;
  
  emfit(model, dat, opts, eres, work);
  
  return List::create(
    Named("alpha") = alpha,
    Named("rate") = rate,
    Named("iter") = opts.iter,
    Named("aerror") = opts.aerror,
    Named("rerror") = opts.rerror,
    Named("llf") = opts.llf,
    Named("convergence") = opts.status == Convergence);
}

// [[Rcpp::export]]
List emfit_cf1_group(NumericVector alpha,
                     NumericVector rate,
                     List data,
                     List options,
                     S4 Q0,
                     S4 P0,
                     S4 H0) {
  using MatrixT = S4matrix<CSCMatrixT>;
  auto Q = MatrixT(Q0);
  auto P = MatrixT(P0);
  auto H = MatrixT(H0);
  
  auto maxiter = as<int>(options["maxiter"]);
  auto atol = as<double>(options["abstol"]);
  auto rtol = as<double>(options["reltol"]);
  auto verbose = as<bool>(options["em.verbose"]);
  auto steps = as<int>(options["steps"]);
  auto ufactor = as<double>(options["uniform.factor"]);
  auto eps = as<double>(options["poisson.eps"]);
  
  int n = alpha.length();
  IntegerVector di(n);
  diag(Q, di);
  copy(Q, P);
  double qv = unif(P, di, ufactor);
  NumericVector al(n);
  copy(alpha, al);
  NumericVector xi(n);
  xi[n-1] = rate[n-1];
  auto gph = GPH<NumericVector, MatrixT, IntegerVector>(al, Q, P, xi, qv, di);
  auto model = CF1<NumericVector, GPH<NumericVector, MatrixT, IntegerVector>>(alpha, rate, gph);
  
  auto tdat = as<NumericVector>(data["time"]);
  auto gdat = as<IntegerVector>(data["counts"]);
  auto idat = as<IntegerVector>(data["indicators"]);
  double maxtime = as<double>(data["maxtime"]);
  int glast = as<int>(data["last"]);
  auto m = tdat.length();
  auto dat = PHGroupSample<NumericVector,IntegerVector,IntegerVector>(tdat, gdat, idat, maxtime, glast);
  
  auto eres = GPHEres<std::vector<double>, MatrixT>(
    std::vector<double>(n),
    std::vector<double>(n),
    std::vector<double>(n),
    H);
  auto work = GPHWorkSpace1<std::vector<double>>(m);
  
  auto opts = EMOptions();
  opts.maxiter = maxiter;
  opts.atol = atol;
  opts.rtol = rtol;
  opts.steps = steps;
  opts.verbose = verbose;
  opts.ufactor = ufactor;
  opts.poisson_eps = eps;
  
  emfit(model, dat, opts, eres, work);
  
  return List::create(
    Named("alpha") = alpha,
    Named("rate") = rate,
    Named("iter") = opts.iter,
    Named("aerror") = opts.aerror,
    Named("rerror") = opts.rerror,
    Named("llf") = opts.llf,
    Named("convergence") = opts.status == Convergence);
}

// [[Rcpp::export]]
List emfit_cf1_group_poi(double omega,
                         NumericVector alpha,
                         NumericVector rate,
                         List data,
                         List options,
                         S4 Q0,
                         S4 P0,
                         S4 H0) {
  using MatrixT = S4matrix<CSCMatrixT>;
  auto Q = MatrixT(Q0);
  auto P = MatrixT(P0);
  auto H = MatrixT(H0);
  
  auto maxiter = as<int>(options["maxiter"]);
  auto atol = as<double>(options["abstol"]);
  auto rtol = as<double>(options["reltol"]);
  auto verbose = as<bool>(options["em.verbose"]);
  auto steps = as<int>(options["steps"]);
  auto ufactor = as<double>(options["uniform.factor"]);
  auto eps = as<double>(options["poisson.eps"]);
  
  int n = alpha.length();
  IntegerVector di(n);
  diag(Q, di);
  copy(Q, P);
  double qv = unif(P, di, ufactor);
  NumericVector al(n);
  copy(alpha, al);
  NumericVector xi(n);
  xi[n-1] = rate[n-1];
  auto gph = GPH<NumericVector, MatrixT, IntegerVector>(al, Q, P, xi, qv, di);
  auto gphpoi = GPHPoi<GPH<NumericVector, MatrixT, IntegerVector>>(gph, omega);
  auto model = CF1<NumericVector, GPHPoi<GPH<NumericVector, MatrixT, IntegerVector>>>(alpha, rate, gphpoi);
  
  auto tdat = as<NumericVector>(data["time"]);
  auto gdat = as<IntegerVector>(data["counts"]);
  auto idat = as<IntegerVector>(data["indicators"]);
  double maxtime = as<double>(data["maxtime"]);
  int glast = as<int>(data["last"]);
  auto m = tdat.length();
  auto dat = PHGroupSample<NumericVector,IntegerVector,IntegerVector>(tdat, gdat, idat, maxtime, glast);
  
  auto eres = GPHEres<std::vector<double>, MatrixT>(
    std::vector<double>(n),
    std::vector<double>(n),
    std::vector<double>(n),
    H);
  auto work = GPHWorkSpace1<std::vector<double>>(m);
  
  auto opts = EMOptions();
  opts.maxiter = maxiter;
  opts.atol = atol;
  opts.rtol = rtol;
  opts.steps = steps;
  opts.verbose = verbose;
  opts.ufactor = ufactor;
  opts.poisson_eps = eps;
  
  emfit(model, dat, opts, eres, work);
  
  return List::create(
    Named("omega") = model.gph.omega,
    Named("alpha") = alpha,
    Named("rate") = rate,
    Named("iter") = opts.iter,
    Named("aerror") = opts.aerror,
    Named("rerror") = opts.rerror,
    Named("llf") = opts.llf,
    Named("convergence") = opts.status == Convergence);
}

/*** R
alpha <- c(0.2, 0.6, 0.2)
Q <- rbind(
  c(-2.0, 2.0, 0.0),
  c(0.0, -5.0, 5.0),
  c(0.0, 0.0, -8.0))
xi <- c(0.0, 0.0, 8.0)
rate <- c(2.0, 5.0, 8.0)
dat <- list(time=c(1,2,1,3,4), weights=c(1,1,1,1,1), maxtime=4)
options <- list(maxiter=10, abstol=1.0e-3, reltol=1.0e-6,
                steps=1, em.verbose=TRUE, uniform.factor=1.01,
                poisson.eps=1.0e-8)
matclass <- "dgCMatrix"
result <- emfit_cf1_wtime(alpha, rate, dat, options, as(Q, matclass), as(Q, matclass), as(Q, matclass))
print(result)

alpha <- c(0.2, 0.6, 0.2)
Q <- rbind(
  c(-2.0, 2.0, 0.0),
  c(0.0, -5.0, 5.0),
  c(0.0, 0.0, -8.0))
xi <- c(0.0, 0.0, 8.0)
rate <- c(2.0, 5.0, 8.0)
dat <- list(time=c(1,2,1,3,4), counts=c(1,3,-1,2,4), indicators=c(0,0,0,1,0), last=10, maxtime=4)
options <- list(maxiter=10, abstol=1.0e-3, reltol=1.0e-6,
                steps=1, em.verbose=TRUE, uniform.factor=1.01,
                poisson.eps=1.0e-8)
matclass <- "dgCMatrix"
result <- emfit_cf1_group(alpha, rate, dat, options, as(Q, matclass), as(Q, matclass), as(Q, matclass))
print(result)

alpha <- c(0.2, 0.6, 0.2)
Q <- rbind(
  c(-2.0, 2.0, 0.0),
  c(0.0, -5.0, 5.0),
  c(0.0, 0.0, -8.0))
xi <- c(0.0, 0.0, 8.0)
rate <- c(2.0, 5.0, 8.0)
omega <- 10
dat <- list(time=c(1,2,1,3,4), counts=c(1,3,-1,2,4), indicators=c(0,0,0,1,0), last=10, maxtime=4)
options <- list(maxiter=10, abstol=1.0e-3, reltol=1.0e-6,
                steps=1, em.verbose=TRUE, uniform.factor=1.01,
                poisson.eps=1.0e-8)
matclass <- "dgCMatrix"
result <- emfit_cf1_group_poi(omega, alpha, rate, dat, options, as(Q, matclass), as(Q, matclass), as(Q, matclass))
print(result)
*/

