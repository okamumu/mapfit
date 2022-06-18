#ifndef MAPFIT_PHASE_ERLANG_H
#define MAPFIT_PHASE_ERLANG_H

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

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename OptionT, typename WorkSpace>
double estep(
    const HErlang<T1,T2>& model,
    const PHWeightSample<T3,T4>& data,
    HErlangEres<T5>& eres,
    OptionT,
    WorkSpace& work) {

  using trait1 = stride_vector_traits<typename WorkSpace::Type>;
  using trait2 = stride_vector_traits<T1>;
  using trait3 = stride_vector_traits<T2,int>;
  using trait4 = stride_vector_traits<T5>;
  const double* alpha = trait2::value(model.alpha);
  const double* rate = trait2::value(model.rate);
  const int* shape = trait3::value(model.shape);
  double* eb = trait4::value(eres.eb);
  double* ew = trait4::value(eres.ew);
  
  const int n = model.size();
  const int m = data.size();
  const double* tdat = stride_vector_traits<T3,double>::value(data.time);
  const double* wdat = stride_vector_traits<T4,double>::value(data.weights);
  
  double scale, tmp;
  double llf = 0.0;
  
  // set Erlang pdf
  tmp = 0.0;
  for (int k=0; k<m; k++) {
    tmp += tdat[k];
    for (int i=0; i<n; i++) {
      double* p0 = trait1::value(work.perl0[k]);
      double* p1 = trait1::value(work.perl1[k]);
      p0[i] = gam::erlang_pdf(shape[i], rate[i], tmp);
      p1[i] = tmp * p0[i];
    }
  }
  
  eres.etotal = 0.0;
  fill(eres.eb, 0.0);
  fill(eres.ew, 0.0);
  for (int k=0; k<m; k++) {
    scale = dot(model.alpha, work.perl0[k]);
    axpy(wdat[k]/scale, work.perl0[k], eres.eb);
    axpy(wdat[k]/scale, work.perl1[k], eres.ew);
    llf += wdat[k] * log(scale);
    eres.etotal += wdat[k];
  }
  for (int i=0; i<n; i++) {
    eb[i] *= alpha[i];
    ew[i] *= alpha[i];
  }
  return llf;
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename OptionT, typename WorkSpace>
double estep(
    const HErlang<T1,T2>& model,
    const PHGroupSample<T3,T4,T5>& data,
    HErlangEres<T6>& eres,
    OptionT,
    WorkSpace& work) {
  
  using trait1 = stride_vector_traits<typename WorkSpace::Type>;
  using trait2 = stride_vector_traits<T1>;
  using trait3 = stride_vector_traits<T2,int>;
  using trait4 = stride_vector_traits<T6>;
  const double* alpha = trait2::value(model.alpha);
  const double* rate = trait2::value(model.rate);
  const int* shape = trait3::value(model.shape);
  double* eb = trait4::value(eres.eb);
  double* ew = trait4::value(eres.ew);
  
  int n = model.size();
  int m = data.size();
  const double* tdat = stride_vector_traits<T3,double>::value(data.time);
  const int* gdat = stride_vector_traits<T4,int>::value(data.counts);
  const int* idat = stride_vector_traits<T5,int>::value(data.indicators);
  int gdatlast = data.last;
  
  
  // set Erlang pdf and cdf
  double tmp = 0.0;
  fill(work.cerl0[0], 0.0);
  fill(work.cerl1[0], 0.0);
  for (int k=0; k<m; k++) {
    tmp += tdat[k];
    for (int i=0; i<n; i++) {
      double* p0 = trait1::value(work.perl0[k]);
      double* p1 = trait1::value(work.perl1[k]);
      double* c0 = trait1::value(work.cerl0[k+1]);
      double* c1 = trait1::value(work.cerl1[k+1]);
      p0[i] = gam::erlang_pdf(shape[i], rate[i], tmp);
      p1[i] = tmp * p0[i];
      c0[i] = gam::erlang_cdf(shape[i], rate[i], tmp);
      c1[i] = (shape[i] / rate[i]) * gam::erlang_cdf(shape[i]+1, rate[i], tmp);
    }
  }
  fill(work.cerl0[m+1], 1.0);
  double* c1 = trait1::value(work.cerl1[m+1]);
  for (int i=0; i<n; i++) {
    c1[i] = shape[i] / rate[i];
  }
  
  // estep
  std::vector<double> tmpv0(n);
  std::vector<double> tmpv1(n);
  double scale;
  double nn = 0.0;
  double uu = 0.0;
  double llf = 0.0;
  eres.etotal = 0.0;
  fill(eres.eb, 0.0);
  fill(eres.ew, 0.0);
  for (int k=0; k<m; k++) {
    if (gdat[k] >= 0 && tdat[k] != 0.0) {
      // tmpv0 = cerl0[k+1] - cerl0[k];
      // tmpv1 = cerl1[k+1] - cerl1[k];
      copy(work.cerl0[k+1], tmpv0);
      copy(work.cerl1[k+1], tmpv1);
      axpy(-1.0, work.cerl0[k], tmpv0);
      axpy(-1.0, work.cerl1[k], tmpv1);
      scale = dot(model.alpha, tmpv0);
      nn += gdat[k];
      uu += scale;
      axpy(gdat[k]/scale, tmpv0, eres.eb);
      axpy(gdat[k]/scale, tmpv1, eres.ew);
      llf += gdat[k] * log(scale) - gam::lfact(gdat[k]);
    }
    if (idat[k] == 1) {
      scale = dot(model.alpha, work.perl0[k]);
      nn += 1.0;
      axpy(1.0/scale, work.perl0[k], eres.eb);
      axpy(1.0/scale, work.perl1[k], eres.ew);
      llf += log(scale);
    }
  }
  if (gdatlast >= 0) {
    // tmpv0 = cerl0[m+1] - cerl0[m];
    // tmpv1 = cerl1[m+1] - cerl1[m];
    copy(work.cerl0[m+1], tmpv0);
    copy(work.cerl1[m+1], tmpv1);
    axpy(-1.0, work.cerl0[m], tmpv0);
    axpy(-1.0, work.cerl1[m], tmpv1);
    scale = dot(model.alpha, tmpv0);
    nn += gdatlast;
    uu += scale;
    axpy(gdatlast/scale, tmpv0, eres.eb);
    axpy(gdatlast/scale, tmpv1, eres.ew);
    llf += gdatlast * log(scale) - gam::lfact(gdatlast);
  }
  
  for (int k=0; k<m; k++) {
    if (gdat[k] == -1) {
      // tmpv0 = cerl0[k] - cerl0[k-1];
      // tmpv1 = cerl1[k] - cerl1[k-1];
      copy(work.cerl0[k+1], tmpv0);
      copy(work.cerl1[k+1], tmpv1);
      axpy(-1.0, work.cerl0[k], tmpv0);
      axpy(-1.0, work.cerl1[k], tmpv1);
      axpy(nn/uu, tmpv0, eres.eb);
      axpy(nn/uu, tmpv1, eres.ew);
    }
  }
  if (gdatlast == -1) {
    // tmpv0 = cerl0[m+1] - cerl0[m];
    // tmpv1 = cerl1[m+1] - cerl1[m];
    copy(work.cerl0[m+1], tmpv0);
    copy(work.cerl1[m+1], tmpv1);
    axpy(-1.0, work.cerl0[m], tmpv0);
    axpy(-1.0, work.cerl1[m], tmpv1);
    axpy(nn/uu, tmpv0, eres.eb);
    axpy(nn/uu, tmpv1, eres.ew);
  }
  llf += gam::lgamma(nn + 1.0) - nn * log(uu);
  
  eres.etotal = nn / uu;
  for (int i=0; i<n; i++) {
    eb[i] *= alpha[i];
    ew[i] *= alpha[i];
  }
  return llf;
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename OptionT, typename WorkSpace>
double estep(
    const HErlangPoi<T1,T2>& model,
    const PHGroupSample<T3,T4,T5>& data,
    HErlangEres<T6>& eres,
    OptionT,
    WorkSpace& work) {
  
  using trait1 = stride_vector_traits<typename WorkSpace::Type>;
  using trait2 = stride_vector_traits<T1>;
  using trait3 = stride_vector_traits<T2,int>;
  using trait4 = stride_vector_traits<T6>;
  const double* alpha = trait2::value(model.alpha);
  const double* rate = trait2::value(model.rate);
  const int* shape = trait3::value(model.shape);
  double* eb = trait4::value(eres.eb);
  double* ew = trait4::value(eres.ew);
  
  int n = model.size();
  int m = data.size();
  const double* tdat = stride_vector_traits<T3,double>::value(data.time);
  const int* gdat = stride_vector_traits<T4,int>::value(data.counts);
  const int* idat = stride_vector_traits<T5,int>::value(data.indicators);
  int gdatlast = data.last;
  double omega = model.omega;


  // set Erlang pdf and cdf
  double tmp = 0.0;
  fill(work.cerl0[0], 0.0);
  fill(work.cerl1[0], 0.0);
  for (int k=0; k<m; k++) {
    tmp += tdat[k];
    for (int i=0; i<n; i++) {
      double* p0 = trait1::value(work.perl0[k]);
      double* p1 = trait1::value(work.perl1[k]);
      double* c0 = trait1::value(work.cerl0[k+1]);
      double* c1 = trait1::value(work.cerl1[k+1]);
      p0[i] = gam::erlang_pdf(shape[i], rate[i], tmp);
      p1[i] = tmp * p0[i];
      c0[i] = gam::erlang_cdf(shape[i], rate[i], tmp);
      c1[i] = (shape[i] / rate[i]) * gam::erlang_cdf(shape[i]+1, rate[i], tmp);
    }
  }
  fill(work.cerl0[m+1], 1.0);
  double* c1 = trait1::value(work.cerl1[m+1]);
  for (int i=0; i<n; i++) {
    c1[i] = shape[i] / rate[i];
  }

  // estep
  std::vector<double> tmpv0(n);
  std::vector<double> tmpv1(n);
  double scale;
  double nn = 0.0;
  double uu = 0.0;
  double llf = 0.0;
  eres.etotal = 0.0;
  fill(eres.eb, 0.0);
  fill(eres.ew, 0.0);
  for (int k=0; k<m; k++) {
    if (gdat[k] >= 0 && tdat[k] != 0.0) {
      // tmpv0 = cerl0[k+1] - cerl0[k];
      // tmpv1 = cerl1[k+1] - cerl1[k];
      copy(work.cerl0[k+1], tmpv0);
      copy(work.cerl1[k+1], tmpv1);
      axpy(-1.0, work.cerl0[k], tmpv0);
      axpy(-1.0, work.cerl1[k], tmpv1);
      scale = dot(model.alpha, tmpv0);
      nn += gdat[k];
      uu += scale;
      axpy(gdat[k]/scale, tmpv0, eres.eb);
      axpy(gdat[k]/scale, tmpv1, eres.ew);
      llf += gdat[k] * log(scale) - gam::lfact(gdat[k]);
    }
    if (idat[k] == 1) {
      scale = dot(model.alpha, work.perl0[k]);
      nn += 1.0;
      axpy(1.0/scale, work.perl0[k], eres.eb);
      axpy(1.0/scale, work.perl1[k], eres.ew);
      llf += log(scale);
    }
  }
  if (gdatlast >= 0) {
    // tmpv0 = cerl0[m+1] - cerl0[m];
    // tmpv1 = cerl1[m+1] - cerl1[m];
    copy(work.cerl0[m+1], tmpv0);
    copy(work.cerl1[m+1], tmpv1);
    axpy(-1.0, work.cerl0[m], tmpv0);
    axpy(-1.0, work.cerl1[m], tmpv1);
    scale = dot(model.alpha, tmpv0);
    nn += gdatlast;
    uu += scale;
    axpy(gdatlast/scale, tmpv0, eres.eb);
    axpy(gdatlast/scale, tmpv1, eres.ew);
    llf += gdatlast * log(scale) - gam::lfact(gdatlast);
  }

  for (int k=0; k<m; k++) {
    if (gdat[k] == -1) {
      // tmpv0 = cerl0[k] - cerl0[k-1];
      // tmpv1 = cerl1[k] - cerl1[k-1];
      copy(work.cerl0[k+1], tmpv0);
      copy(work.cerl1[k+1], tmpv1);
      axpy(-1.0, work.cerl0[k], tmpv0);
      axpy(-1.0, work.cerl1[k], tmpv1);
      axpy(omega, tmpv0, eres.eb);
      axpy(omega, tmpv1, eres.ew);
    }
  }
  if (gdatlast == -1) {
    // tmpv0 = cerl0[m+1] - cerl0[m];
    // tmpv1 = cerl1[m+1] - cerl1[m];
    copy(work.cerl0[m+1], tmpv0);
    copy(work.cerl1[m+1], tmpv1);
    axpy(-1.0, work.cerl0[m], tmpv0);
    axpy(-1.0, work.cerl1[m], tmpv1);
    axpy(omega, tmpv0, eres.eb);
    axpy(omega, tmpv1, eres.ew);
  }
  llf += nn * log(omega) - omega * uu;

  eres.etotal = nn + omega * (1.0 - uu);
  for (int i=0; i<n; i++) {
    eb[i] *= alpha[i];
    ew[i] *= alpha[i];
  }
  return llf;
}

template <typename V, typename I, typename V2, typename OptionT>
void mstep(const HErlangEres<V2>& eres,
           HErlang<V,I>& model,
           OptionT) {

  using trait2 = stride_vector_traits<V>;
  using trait3 = stride_vector_traits<I,int>;
  using trait4 = stride_vector_traits<V2>;
  double* alpha = trait2::value(model.alpha);
  double* rate = trait2::value(model.rate);
  int* shape = trait3::value(model.shape);
  const double* eb = trait4::value(eres.eb);
  const double* ew = trait4::value(eres.ew);

  int n = model.size();
  copy(eres.eb, model.alpha);
  scal(1.0/eres.etotal, model.alpha);
  for (int i=0; i<n; i++) {
    rate[i] = shape[i] * eb[i] / ew[i];
  }
}

template <typename V, typename I, typename V2, typename OptionT>
void mstep(const HErlangEres<V2>& eres,
           HErlangPoi<V,I>& model,
           OptionT) {

  using trait2 = stride_vector_traits<V>;
  using trait3 = stride_vector_traits<I,int>;
  using trait4 = stride_vector_traits<V2>;
  double* alpha = trait2::value(model.alpha);
  double* rate = trait2::value(model.rate);
  int* shape = trait3::value(model.shape);
  const double* eb = trait4::value(eres.eb);
  const double* ew = trait4::value(eres.ew);

  int n = model.size();
  copy(eres.eb, model.alpha);
  scal(1.0/eres.etotal, model.alpha);
  model.omega = eres.etotal;
  for (int i=0; i<n; i++) {
    rate[i] = shape[i] * eb[i] / ew[i];
  }
}

#endif
