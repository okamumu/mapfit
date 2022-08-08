#ifndef MAPFIT_PHASE_ERLANG_GROUP_H
#define MAPFIT_PHASE_ERLANG_GROUP_H

#include <Rcpp.h>
#include "traits.h"
#include "gamma.h"
#include "blas.h"
#include "phase_data.h"
#include "phase_models.h"

#define TDAT(k) (tdat[(k)-1])
#define GDAT(k) (gdat[(k)-1])
#define IDAT(k) (idat[(k)-1])

struct HErlangWorkSpaceGroup {
  std::vector<std::vector<double>> perl0;
  std::vector<std::vector<double>> perl1;
  std::vector<std::vector<double>> cerl0;
  std::vector<std::vector<double>> cerl1;
  
  HErlangWorkSpaceGroup(int m, int n) :
    perl0(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    perl1(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    cerl0(std::vector<std::vector<double>>(m+2, std::vector<double>(n))),
    cerl1(std::vector<std::vector<double>>(m+2, std::vector<double>(n)))
  {}
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename OptionT, typename WorkSpace>
double estep(
    const HErlang<T1,T2>& model,
    const PHGroupSample<T3,T4,T5>& data,
    HErlangEres<T6>& eres,
    OptionT,
    WorkSpace& work) {
  
  int n = model.size();
  int m = data.size();
  const double* alpha = stride_vector_traits<T1,double>::value(model.alpha);
  const double* rate = stride_vector_traits<T1,double>::value(model.rate);
  const int* shape = stride_vector_traits<T2,int>::value(model.shape);
  const double* tdat = stride_vector_traits<T3,double>::value(data.time);
  const int* gdat = stride_vector_traits<T4,int>::value(data.counts);
  const int* idat = stride_vector_traits<T5,int>::value(data.indicators);
  int gdatlast = data.last;

  // workspace
  std::vector<double> tmpv0(n);
  std::vector<double> tmpv1(n);
  // std::vector<std::vector<double>> perl0(m+1, std::vector<double>(n));
  // std::vector<std::vector<double>> perl1(m+1, std::vector<double>(n));
  // std::vector<std::vector<double>> cerl0(m+2, std::vector<double>(n));
  // std::vector<std::vector<double>> cerl1(m+2, std::vector<double>(n));
  std::vector<std::vector<double>>& perl0(work.perl0);
  std::vector<std::vector<double>>& perl1(work.perl1);
  std::vector<std::vector<double>>& cerl0(work.cerl0);
  std::vector<std::vector<double>>& cerl1(work.cerl1);
  
  // set Erlang pdf and cdf
  double tmp = 0.0;
  fill(cerl0[0], 0.0);
  fill(cerl1[0], 0.0);
  for (int k=1; k<=m; k++) {
    tmp += TDAT(k);
    for (int i=0; i<n; i++) {
      perl0[k][i] = gam::erlang_pdf(shape[i], rate[i], tmp);
      perl1[k][i] = tmp * perl0[k][i];
      cerl0[k][i] = gam::erlang_cdf(shape[i], rate[i], tmp);
      cerl1[k][i] = (shape[i] / rate[i]) * gam::erlang_cdf(shape[i]+1, rate[i], tmp);
    }
  }
  fill(cerl0[m+1], 1.0);
  for (int i=0; i<n; i++) {
    cerl1[m+1][i] = shape[i] / rate[i];
  }
  
  // estep
  double scale;
  double nn = 0.0;
  double uu = 0.0;
  double llf = 0.0;
  eres.etotal = 0.0;
  fill(eres.eb, 0.0);
  fill(eres.ew, 0.0);
  for (int k=1; k<=m; k++) {
    if (GDAT(k) >= 0 && TDAT(k) != 0.0) {
      // tmpv0 = cerl0[k] - cerl0[k-1];
      // tmpv1 = cerl1[k] - cerl1[k-1];
      copy(cerl0[k], tmpv0);
      copy(cerl1[k], tmpv1);
      axpy(-1.0, cerl0[k-1], tmpv0);
      axpy(-1.0, cerl1[k-1], tmpv1);
      scale = dot(model.alpha, tmpv0);
      nn += GDAT(k);
      uu += scale;
      axpy(GDAT(k)/scale, tmpv0, eres.eb);
      axpy(GDAT(k)/scale, tmpv1, eres.ew);
      llf += GDAT(k) * log(scale) - gam::lfact(GDAT(k));
    }
    if (IDAT(k) == 1) {
      scale = dot(model.alpha, perl0[k]);
      nn += 1.0;
      axpy(1.0/scale, perl0[k], eres.eb);
      axpy(1.0/scale, perl1[k], eres.ew);
      llf += log(scale);
    }
  }
  if (gdatlast >= 0) {
    // tmpv0 = cerl0[m+1] - cerl0[m];
    // tmpv1 = cerl1[m+1] - cerl1[m];
    copy(cerl0[m+1], tmpv0);
    copy(cerl1[m+1], tmpv1);
    axpy(-1.0, cerl0[m], tmpv0);
    axpy(-1.0, cerl1[m], tmpv1);
    scale = dot(model.alpha, tmpv0);
    nn += gdatlast;
    uu += scale;
    axpy(gdatlast/scale, tmpv0, eres.eb);
    axpy(gdatlast/scale, tmpv1, eres.ew);
    llf += gdatlast * log(scale) - gam::lfact(gdatlast);
  }
  
  for (int k=1; k<=m; k++) {
    if (GDAT(k) == -1) {
      // tmpv0 = cerl0[k] - cerl0[k-1];
      // tmpv1 = cerl1[k] - cerl1[k-1];
      copy(cerl0[k], tmpv0);
      copy(cerl1[k], tmpv1);
      axpy(-1.0, cerl0[k-1], tmpv0);
      axpy(-1.0, cerl1[k-1], tmpv1);
      axpy(nn/uu, tmpv0, eres.eb);
      axpy(nn/uu, tmpv1, eres.ew);
    }
  }
  if (gdatlast == -1) {
    // tmpv0 = cerl0[m+1] - cerl0[m];
    // tmpv1 = cerl1[m+1] - cerl1[m];
    copy(cerl0[m+1], tmpv0);
    copy(cerl1[m+1], tmpv1);
    axpy(-1.0, cerl0[m], tmpv0);
    axpy(-1.0, cerl1[m], tmpv1);
    axpy(nn/uu, tmpv0, eres.eb);
    axpy(nn/uu, tmpv1, eres.ew);
  }
  llf += gam::lgamma(nn + 1.0) - nn * log(uu);
  
  double* eb = stride_vector_traits<T6>::value(eres.eb);
  double* ew = stride_vector_traits<T6>::value(eres.ew);
  eres.etotal = nn / uu;
  for (int i=0; i<n; i++) {
    eb[i] *= alpha[i];
    ew[i] *= alpha[i];
  }
  return llf;
}

#endif
