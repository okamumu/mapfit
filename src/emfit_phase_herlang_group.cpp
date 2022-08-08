#include <Rcpp.h>
using namespace Rcpp;

#include "phase_herlang_estep_group.h"
#include "phase_herlang_mstep.h"
#include "emfit.h"

// [[Rcpp::export]]
List emfit_herlang_group(NumericVector alpha,
                         IntegerVector shape,
                         NumericVector rate,
                         List data,
                         List options) {
  auto model = HErlang<NumericVector, IntegerVector>(alpha, shape, rate);
  auto n = model.size();
  auto tdat = as<NumericVector>(data["intervals"]);
  auto gdat = as<IntegerVector>(data["counts"]);
  auto idat = as<IntegerVector>(data["instants"]);
  double maxtime = as<double>(data["maxinterval"]);
  int glast = as<int>(data["lastcount"]);
  auto m = tdat.length();
  auto dat = PHGroupSample<NumericVector,IntegerVector,IntegerVector>(tdat, gdat, idat, maxtime, glast);
  auto eres = HErlangEres<std::vector<double>>(std::vector<double>(n), std::vector<double>(n));
  auto work = HErlangWorkSpaceGroup(m, n);
  
  auto maxiter = as<int>(options["maxiter"]);
  auto atol = as<double>(options["abstol"]);
  auto rtol = as<double>(options["reltol"]);
  auto verbose = as<bool>(options["em.verbose"]);
  auto steps = as<int>(options["steps"]);
  
  auto opts = EMOptions();
  opts.maxiter = maxiter;
  opts.atol = atol;
  opts.rtol = rtol;
  opts.steps = steps;
  opts.verbose = verbose;
  
  emfit(model, dat, opts, eres, work);
  
  return List::create(
    Named("alpha") = model.alpha,
    Named("rate") = model.rate,
    Named("shape") = model.shape,
    Named("iter") = opts.iter,
    Named("aerror") = opts.aerror,
    Named("rerror") = opts.rerror,
    Named("llf") = opts.llf,
    Named("convergence") = opts.status == Convergence);
}

/*** R
alpha <- c(0.4, 0.6)
rate <- c(1.0, 2.0)
shape <- c(1, 2)
dat <- list(intervals=c(1,2,1,3,4), counts=c(1,3,-1,2,4), instants=c(0,0,0,1,0), lastcount=10, maxinterval=4)
options <- list(maxiter=10, abstol=1.0e-3, reltol=1.0e-6, steps=5, em.verbose=TRUE)
result <- emfit_herlang_group(alpha, shape, rate, dat, options)
print(result)
*/
