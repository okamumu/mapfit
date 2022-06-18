#include <Rcpp.h>
using namespace Rcpp;

#include "map_gen.h"
#include "emfit.h"

// [[Rcpp::export]]
List emfit_mapgen_group(
    NumericVector alpha,
    NumericVector xi,
    NumericMatrix D0,
    NumericMatrix D1,
    List data,
    List options,
    NumericMatrix P0,
    NumericMatrix P1,
    NumericMatrix H0,
    NumericMatrix H1,
    NumericMatrix en0,
    NumericMatrix en1) {
  using MatrixT = NumericMatrix;
  // auto Q = MatrixT(Q0);
  // auto P = MatrixT(P0);
  // auto H = MatrixT(H0);
  
  auto maxiter = as<int>(options["maxiter"]);
  auto atol = as<double>(options["abstol"]);
  auto rtol = as<double>(options["reltol"]);
  auto verbose = as<bool>(options["em.verbose"]);
  auto steps = as<int>(options["steps"]);
  auto ufactor = as<double>(options["uniform.factor"]);
  auto eps = as<double>(options["poisson.eps"]);

  int n = alpha.length();
  IntegerVector di(n);
  diag(D0, di);
  copy(D0, P0);
  copy(D1, P1);
  double qv = unif(P0, di, ufactor);
  scal(1.0/qv, P1);
  auto model = MAP<NumericVector,MatrixT,IntegerVector>(alpha, xi, D0, D1, P0, P1, di, qv);
  
  auto tdat = as<NumericVector>(data["time"]);
  auto gdat = as<IntegerVector>(data["counts"]);
  auto idat = as<IntegerVector>(data["indicators"]);
  double maxtime = as<double>(data["maxtime"]);
  int maxcount = as<int>(data["maxcount"]);
  auto m = tdat.length();
  auto dat = MAPGroupSample<NumericVector,IntegerVector,IntegerVector>(tdat, gdat, idat, maxtime, maxcount);
  
  auto eres = MAPEres<NumericVector,MatrixT>(
    NumericVector(n),
    NumericVector(n),
    en0, en1);
  auto work = MAPWorkSpace<MatrixT>(H0, H1);
  
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
    Named("xi") = xi,
    Named("D0") = D0,
    Named("D1") = D1,
    Named("iter") = opts.iter,
    Named("aerror") = opts.aerror,
    Named("rerror") = opts.rerror,
    Named("llf") = opts.llf,
    Named("convergence") = opts.status == Convergence);
}

/*** R
alpha <- c(0.4, 0.6)
xi <- c(1,1)
D0 <- rbind(
  c(-4, 2),
  c(1, -5)
)
D1 <- rbind(
  c(2, 0),
  c(1, 3)
)
P0 <- matrix(0, 2, 2)
P1 <- matrix(0, 2, 2)
H0 <- matrix(0, 2, 2)
H1 <- matrix(0, 2, 2)
en0 <- matrix(0, 2, 2)
en1 <- matrix(0, 2, 2)
dat <- list(time=c(1,2,1,3,4), counts=c(1,3,0,2,4), indicators=c(0,0,0,0,0), maxtime=4, maxcount=4)
options <- list(maxiter=10, abstol=1.0e-3, reltol=1.0e-6,
                steps=1, em.verbose=TRUE, uniform.factor=1.01,
                poisson.eps=1.0e-8)
result <- emfit_mapgen_group(alpha, xi, D0, D1, dat, options,
                             P0, P1, H0, H1, en0, en1)
print(result)
*/

