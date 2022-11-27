#ifndef MAPFIT_ESTEP_PHASE_GROUP_GPH
#define MAPFIT_ESTEP_PHASE_GROUP_GPH

#include <Rcpp.h>
#include "traits.h"
#include "poisson.h"
#include "gamma.h"
#include "blas.h"
#include "unif.h"
#include "phase_data.h"
#include "phase_models.h"

#define TDAT(k) (tdat[(k)-1])
#define GDAT(k) (gdat[(k)-1])
#define IDAT(k) (idat[(k)-1])

struct GPHWorkSpaceGroup {
  std::vector<std::vector<double>> barvf;
  std::vector<std::vector<double>> barvb;
  std::vector<std::vector<double>> vb;
  std::vector<std::vector<double>> vc;
  
  GPHWorkSpaceGroup(int m, int n) :
    barvf(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    barvb(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    vb(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    vc(std::vector<std::vector<double>>(m+1, std::vector<double>(n)))
  {}
};

template <typename T0, typename T1, typename T2, typename T4,
          typename T5, typename T6, typename OptionT>
double llf(
    const GPH<T1,T2,T0>& model,
    const PHGroupSample<T4,T5,T6>& data,
    OptionT& poisson_eps) noexcept {
  
  const int m = data.size();
  const double* tdat = stride_vector_traits<T4,double>::value(data.time);
  const int* gdat = stride_vector_traits<T5,int>::value(data.counts);
  const int* idat = stride_vector_traits<T6,int>::value(data.indicators);
  const int gdatlast = data.last;
  const double tmax = data.maxtime;
  
  int n = model.size();
  double qv = model.qv;
  std::vector<double> vone(n, 1.0);

  // work
  int right = poi::rightbound(qv*tmax, poisson_eps) + 1;
  std::vector<double> prob(right+1);
  std::vector<double> barvb0(n);
  std::vector<double> barvb1(n);
  std::vector<double> tmpvb(n);
  std::vector<double> tmpv(n);

  double scale;
  double llf = 0.0;
  int nn = 0;
  double uu = 0.0;
  
  // forward & backward
  copy(vone, barvb0);
  for (int k=1; k<=m; k++) {
    // barvf[k] = barvf[k-1] * exp(T * tdat[k])
    // barvb[k] = exp(T * tdat[k]) * barvb[k-1]
    
    int right = poi::rightbound(qv*TDAT(k), poisson_eps) + 1;
    double weight = poi::pmf(qv*TDAT(k), 0, right, prob);
    
    fill(barvb1, 0.0);
    copy(barvb0, tmpvb);
    axpy(prob[0], tmpvb, barvb1);
    for (int u=1; u<=right; u++) {
      gemv(NOTRANS{}, 1.0, model.P, tmpvb, 0.0, tmpv);
      copy(tmpv, tmpvb);
      axpy(prob[u], tmpvb, barvb1);
    }
    scal(1.0/weight, barvb1);
    
    
    // tmpvb = barvb[k-1] - barvb[k]
    copy(barvb0, tmpvb);
    axpy(-1.0, barvb1, tmpvb);
    
    if (GDAT(k) >= 0 && TDAT(k) > 0.0) {
      double tmp = dot(model.alpha, tmpvb);
      llf += GDAT(k) * log(tmp) - gam::lfact(GDAT(k));
      nn += GDAT(k);
      uu += tmp;
    }
    if (IDAT(k) == 1) {
      // vb[k] = (-ph.T) * barvb[k]
      gemv(NOTRANS{}, -1.0, model.Q, barvb1, 0.0, tmpvb);
      double tmp = dot(model.alpha, tmpvb);
      llf += log(tmp);
      nn += 1;
    }
    
    copy(barvb1, barvb0);
  }
  // for the interval [t_m, infinity)
  if (gdatlast >= 0) {
    double tmp = dot(model.alpha, barvb1);
    llf += gdatlast * log(tmp) - gam::lfact(gdatlast);
    nn += gdatlast;
    uu += tmp;
  }
  llf += gam::lfact(nn) - nn * log(uu);
  
  return llf;
}

template <typename T0, typename T1, typename T2, typename T4,
          typename T5, typename T6, typename T7, typename T8,
          typename OptionT, typename WorkSpace>
double estep(
    const GPH<T1,T2,T0>& model,
    const PHGroupSample<T4,T5,T6>& data,
    GPHEres<T7,T8>& eres,
    OptionT& options,
    WorkSpace& work) noexcept {
  int n = model.size();
  std::vector<double> baralpha(n);
  gesv(TRANS{}, -1.0, model.Q, model.alpha, baralpha);
  return estep_group(model, baralpha, data, eres, options, work);
}

template <typename T0, typename T1, typename T2, typename T4,
          typename T5, typename T6, typename T7, typename T8,
          typename T9,
          typename OptionT, typename WorkSpace>
double estep_group(
    const GPH<T1,T2,T0>& model,
    const T9& baralpha,
    const PHGroupSample<T4,T5,T6>& data,
    GPHEres<T7,T8>& eres,
    OptionT& options,
    WorkSpace& work) noexcept {
  
  const int m = data.size();
  const double* tdat = stride_vector_traits<T4,double>::value(data.time);
  const int* gdat = stride_vector_traits<T5,int>::value(data.counts);
  const int* idat = stride_vector_traits<T6,int>::value(data.indicators);
  const int gdatlast = data.last;
  const double tmax = data.maxtime;
  
  int n = model.size();
  double qv = model.qv;
  std::vector<double> vone(n, 1.0);

  // work
  int right = poi::rightbound(qv*tmax, options.poisson_eps) + 1;
  std::vector<double> prob(right+1);
  std::vector<std::vector<double>> vx(right+1, std::vector<double>(n, 0));
  std::vector<double> tmpvf(n);
  std::vector<double> tmpvb(n);
  std::vector<double> tmpv(n);
  // std::vector<std::vector<double>> barvf(m+1, std::vector<double>(n));
  // std::vector<std::vector<double>> barvb(m+1, std::vector<double>(n));
  // std::vector<std::vector<double>> vb(m+1, std::vector<double>(n));
  // std::vector<std::vector<double>> vc(m+1, std::vector<double>(n));
  std::vector<double> wg(m+2);
  std::vector<double> wp(m+2);

  std::vector<std::vector<double>>& barvf(work.barvf);
  std::vector<std::vector<double>>& barvb(work.barvb);
  std::vector<std::vector<double>>& vb(work.vb);
  std::vector<std::vector<double>>& vc(work.vc);

  double scale;
  double llf = 0.0;
  int nn = 0;
  double uu = 0.0;

  eres.etotal = 0.0;
  fill(eres.eb, 0.0);
  fill(eres.ey, 0.0);
  fill(eres.ez, 0.0);
  fill(eres.en, 0.0);
  
  // forward & backward
  copy(baralpha, barvf[0]);
  copy(vone, barvb[0]);
  copy(model.xi, vb[0]);
  for (int k=1; k<=m; k++) {
    // barvf[k] = barvf[k-1] * exp(T * tdat[k])
    // barvb[k] = exp(T * tdat[k]) * barvb[k-1]

    int right = poi::rightbound(qv*TDAT(k), options.poisson_eps) + 1;
    double weight = poi::pmf(qv*TDAT(k), 0, right, prob);

    fill(barvf[k], 0.0);
    fill(barvb[k], 0.0);
    copy(barvf[k-1], tmpvf);
    copy(barvb[k-1], tmpvb);
    axpy(prob[0], tmpvf, barvf[k]);
    axpy(prob[0], tmpvb, barvb[k]);
    for (int u=1; u<=right; u++) {
      gemv(TRANS{}, 1.0, model.P, tmpvf, 0.0, tmpv);
      copy(tmpv, tmpvf);
      gemv(NOTRANS{}, 1.0, model.P, tmpvb, 0.0, tmpv);
      copy(tmpv, tmpvb);
      axpy(prob[u], tmpvf, barvf[k]);
      axpy(prob[u], tmpvb, barvb[k]);
    }
    scal(1.0/weight, barvf[k]);
    scal(1.0/weight, barvb[k]);

    // vb[k] = (-ph.T) * barvb[k]
    gemv(NOTRANS{}, -1.0, model.Q, barvb[k], 0.0, vb[k]);

    // tmpvf = barvf[k-1] - barvf[k]
    // tmpvb = barvb[k-1] - barvb[k]
    copy(barvf[k-1], tmpvf);
    copy(barvb[k-1], tmpvb);
    axpy(-1.0, barvf[k], tmpvf);
    axpy(-1.0, barvb[k], tmpvb);

    if (GDAT(k) >= 0 && TDAT(k) > 0.0) {
      double tmp = dot(model.alpha, tmpvb);
      llf += GDAT(k) * log(tmp) - gam::lfact(GDAT(k));
      nn += GDAT(k);
      uu += tmp;
      wg[k] = GDAT(k) / tmp;
      axpy(wg[k], tmpvb, eres.eb);
      axpy(wg[k], tmpvf, eres.ey);
    }
    if (IDAT(k) == 1) {
      double tmp = dot(model.alpha, vb[k]);
      llf += log(tmp);
      nn += 1;
      wp[k] = 1 / tmp;
      axpy(wp[k], vb[k], eres.eb);
      gemv(TRANS{}, -wp[k], model.Q, barvf[k], 1.0, eres.ey);
    } else {
      wp[k] = 0.0;
    }
    
  }
  // for the interval [t_m, infinity)
  if (gdatlast >= 0) {
    double tmp = dot(model.alpha, barvb[m]);
    llf += gdatlast * log(tmp) - gam::lfact(gdatlast);
    nn += gdatlast;
    uu += tmp;
    wg[m+1] = gdatlast / tmp;
    axpy(wg[m+1], barvb[m], eres.eb);
    axpy(wg[m+1], barvf[m], eres.ey);
  }
  for (int k=1; k<=m; k++) {
    if (GDAT(k) == -1) {
      wg[k] = nn / uu;
      copy(barvf[k-1], tmpvf);
      copy(barvb[k-1], tmpvb);
      axpy(-1.0, barvf[k], tmpvf);
      axpy(-1.0, barvb[k], tmpvb);
      axpy(wg[k], tmpvb, eres.eb);
      axpy(wg[k], tmpvf, eres.ey);
    }
  }
  if (gdatlast == -1) {
    wg[m+1] = nn / uu;
    axpy(wg[m+1], barvb[m], eres.eb);
    axpy(wg[m+1], barvf[m], eres.ey);
  }
  llf += gam::lfact(nn) - nn * log(uu);

  // compute vectors and convolution

  fill(vc[m], 0.0);
  axpy(wg[m+1] - wg[m], baralpha, vc[m]);
  if (IDAT(m) == 1) {
    axpy(wp[m], model.alpha, vc[m]);
  }
  for (int k=m-1; k>=1; k--) {
    // vc[k] = vc[k+1] * exp(T * tdat[k+1]) + (wg[k+1] - wg[k]) * baralpha + I(idat[k]==1) (wp[k] * alpha)
    int right = poi::rightbound(qv*TDAT(k+1), options.poisson_eps) + 1;
    double weight = poi::pmf(qv*TDAT(k+1), 0, right, prob);

    fill(vc[k], 0.0);
    copy(vc[k+1], tmpvf);
    axpy(prob[0], tmpvf, vc[k]);
    for (int u=1; u<=right; u++) {
      gemv(TRANS{}, 1.0, model.P, tmpvf, 0.0, tmpv);
      copy(tmpv, tmpvf);
      axpy(prob[u], tmpvf, vc[k]);
    }
    scal(1.0/weight, vc[k]);
    axpy(wg[k+1]-wg[k], baralpha, vc[k]);
    if (IDAT(k) == 1) {
      axpy(wp[k], model.alpha, vc[k]);
    }
  }

  for (int k=1; k<=m; k++) {
    // compute convolution integral
    // int_0^tdat[k] exp(T* s) * vb[k-1] * vc[k] * exp(T(tdat[k]-s)) ds
    int right = poi::rightbound(qv*TDAT(k), options.poisson_eps) + 1;
    double weight = poi::pmf(qv*TDAT(k), 0, right, prob);

    fill(vx[right], 0.0);
    axpy(prob[right], vb[k-1], vx[right]);
    for (int l=right-1; l>=1; l--) {
      gemv(NOTRANS{}, 1.0, model.P, vx[l+1], 0.0, vx[l]);
      axpy(prob[l], vb[k-1], vx[l]);
    }

    ger(NOTRANS{}, 1.0/(qv*weight), vc[k], vx[1], eres.en);
    for (int l=1; l<=right-1; l++) {
      gemv(TRANS{}, 1.0, model.P, vc[k], 0.0, tmpv);
      copy(tmpv, vc[k]);
      ger(NOTRANS{}, 1.0/(qv*weight), vc[k], vx[l+1], eres.en);
    }
    ger(NOTRANS{}, wg[k+1]-wg[k], baralpha, barvb[k], eres.en);
  }
  ger(NOTRANS{}, wg[1], baralpha, barvb[0], eres.en);

  const double* alpha = vector_traits<T1>::value(model.alpha);
  const double* Q = vector_traits<T2>::value(model.Q);
  const double* xi = vector_traits<T1>::value(model.xi);
  const int* diag = vector_traits<T0,int>::value(model.diag);
  double* eb = vector_traits<T7>::value(eres.eb);
  double* ey = vector_traits<T7>::value(eres.ey);
  double* ez = vector_traits<T7>::value(eres.ez);
  double* en = vector_traits<T8>::value(eres.en);

  eres.etotal = nn / uu;
  for (int i=0; i<n; i++) {
    eb[i] *= alpha[i];
    ey[i] *= xi[i];
    ez[i] = en[diag[i]];
  }
  for (int i=0; i<vector_traits<T8>::size(eres.en); i++) {
    en[i] *= Q[i];
  }
  return llf;
}

#endif

