#ifndef MAPFIT_PHASE_ERLANG_ESTEP_WTIME_H
#define MAPFIT_PHASE_ERLANG_ESTEP_WTIME_H

#include <Rcpp.h>
#include "traits.h"
#include "gamma.h"
#include "blas.h"
#include "phase_data.h"
#include "phase_models.h"

/**
 Description: estep for Erlang-PH with weighted time and group/truncated data
 
 alpha    (in): initial vector
 shape    (in): shape parameter vector
 rate     (in): rate parameter vector
 tdat     (in): interarrival time
 wdat     (in): weights for interarrivals
 gdat     (in): # of arrivals (-1 means NA)
 gdatlast (in): # of arrivals in [lasttime, infinity] (-1 means NA)
 idat     (in): indicator whether an arrival occurs at the last instant
 etotal  (out): expected # of arrivals
 eb      (out): expected # of starts
 ew      (out): expected sojourn time?
 
 return value -> llf (log-likelihood)
 */

#define TDAT(k) (tdat[(k)-1])
#define WDAT(k) (wdat[(k)-1])

struct HErlangWorkSpaceWTime {
  std::vector<std::vector<double>> perl0;
  std::vector<std::vector<double>> perl1;
  
  HErlangWorkSpaceWTime(int m, int n) :
    perl0(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    perl1(std::vector<std::vector<double>>(m+1, std::vector<double>(n)))
  {}
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename OptionT, typename WorkSpace>
double estep(
    const HErlang<T1,T2>& model,
    const PHWeightSample<T3,T4>& data,
    HErlangEres<T5>& eres,
    OptionT,
    WorkSpace& work) {

  const int n = model.size();
  const int m = data.size();
  const double* alpha = stride_vector_traits<T1,double>::value(model.alpha);
  const double* rate = stride_vector_traits<T1,double>::value(model.rate);
  const int* shape = stride_vector_traits<T2,int>::value(model.shape);
  const double* tdat = stride_vector_traits<T3,double>::value(data.time);
  const double* wdat = stride_vector_traits<T4,double>::value(data.weights);
  
  // work
  // std::vector<std::vector<double>> perl0(m+1, std::vector<double>(n));
  // std::vector<std::vector<double>> perl1(m+1, std::vector<double>(n));
  std::vector<std::vector<double>>& perl0(work.perl0);
  std::vector<std::vector<double>>& perl1(work.perl1);
  
  double scale, tmp;
  double llf = 0.0;
  
  // set Erlang pdf
  tmp = 0.0;
  for (int k=1; k<=m; k++) {
    tmp += TDAT(k);
    for (int i=0; i<n; i++) {
      perl0[k][i] = gam::erlang_pdf(shape[i], rate[i], tmp);
      perl1[k][i] = tmp * perl0[k][i];
    }
  }
  
  eres.etotal = 0.0;
  fill(eres.eb, 0.0);
  fill(eres.ew, 0.0);
  for (int k=1; k<=m; k++) {
    scale = dot(model.alpha, perl0[k]);
    axpy(WDAT(k)/scale, perl0[k], eres.eb);
    axpy(WDAT(k)/scale, perl1[k], eres.ew);
    llf += WDAT(k) * log(scale);
    eres.etotal += WDAT(k);
  }

  double* eb = stride_vector_traits<T5>::value(eres.eb);
  double* ew = stride_vector_traits<T5>::value(eres.ew);
  for (int i=0; i<n; i++) {
    eb[i] *= alpha[i];
    ew[i] *= alpha[i];
  }
  return llf;
}

#endif
