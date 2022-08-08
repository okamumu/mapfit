#include <Rcpp.h>
using namespace Rcpp;

#include "phase_gen_estep_leftright.h"
#include "phase_gen_mstep.h"
#include "emfit.h"

 // [[Rcpp::export]]
 List emfit_gph_leftright(NumericVector alpha,
                          S4 Q0,
                          NumericVector xi,
                          List data,
                          List options,
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
   auto model = GPH<NumericVector, MatrixT, IntegerVector>(alpha, Q, P, xi, qv, di);
   
   auto tdat = as<NumericVector>(data["intervals"]);
   auto nu = as<IntegerVector>(data["nu"]);
   double maxtime = as<double>(data["maxinterval"]);
   auto m = tdat.length();
   auto dat = PHLeftRightSample<NumericVector,IntegerVector>(tdat, nu, maxtime);
   
   auto eres = GPHEres<std::vector<double>, MatrixT>(
     std::vector<double>(n),
     std::vector<double>(n),
     std::vector<double>(n),
     H);
   auto work = GPHWorkSpaceLeftRight(m, n);
   
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
     Named("Q") = Q0,
     Named("xi") = xi,
     Named("iter") = opts.iter,
     Named("aerror") = opts.aerror,
     Named("rerror") = opts.rerror,
     Named("llf") = opts.llf,
     Named("convergence") = opts.status == Convergence);
 }
 
/*** R
alpha <- c(0.2, 0.6, 0.2)
Q <- rbind(
  c(-2.0, 1.0, 0.0),
  c(2.0, -5.0, 1.0),
  c(3.0, 2.0, -8.0)) / 10.0
xi <- c(1.0, 2.0, 3.0) / 10.0
dat <- list(intervals=c(0.0, 0.0, 1.0, 0.2, 0.8, 2.5), nu=c(3, 3, 0, 1, 3, 1),
            maxinterval=2.5)
options <- list(maxiter=10, abstol=1.0e-3, reltol=1.0e-6,
                steps=1, em.verbose=TRUE, uniform.factor=1.01,
                poisson.eps=1.0e-8)
result <- try(emfit_gph_leftright(alpha, as(Q, "dgCMatrix"), xi, dat, options, as(Q, "dgCMatrix"), as(Q, "dgCMatrix")))
print(result)
*/

