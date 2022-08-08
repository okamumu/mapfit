#include <Rcpp.h>
using namespace Rcpp;

#include "../phase_herlang_estep_wtime.h"
#include "../phase_herlang_mstep.h"

// [[Rcpp::export]]
void test_estep_wtime(NumericVector alpha,
                      IntegerVector shape,
                      NumericVector rate,
                      List data) {
  int n = 2;
  int m = 5;
  HErlangWorkSpaceWTime work(m, n);
  auto tdat = as<NumericVector>(data["time"]);
  auto wdat = as<NumericVector>(data["weights"]);
  double maxtime = 10.0;
  auto eres = HErlangEres<std::vector<double>>(std::vector<double>(n), std::vector<double>(n));
  double llf;
  auto dat = PHWeightSample<NumericVector,NumericVector>(tdat, wdat, maxtime);
  auto model = HErlang<NumericVector, IntegerVector>(alpha, shape, rate);
  llf = estep(model, dat, eres, nullptr, work);
  Rcout << llf << std::endl;
  mstep(eres, model, nullptr);
  llf = estep(model, dat, eres, nullptr, work);
  Rcout << llf << std::endl;
}

/*** R
alpha <- c(0.4, 0.6)
rate <- c(1.0, 2.0)
shape <- c(1, 2)
dat <- list(size=5, time=c(1,2,1,3,4), weights=c(1,3,4,2,4))
test_estep_wtime(alpha, shape, rate, dat)

*/
