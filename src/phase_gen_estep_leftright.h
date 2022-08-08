#ifndef MAPFIT_PHASE_ESTEP_LEFTRIGHT_GPH
#define MAPFIT_PHASE_ESTEP_LEFTRIGHT_GPH

#include <Rcpp.h>
#include "traits.h"
#include "poisson.h"
#include "gamma.h"
#include "blas.h"
#include "unif.h"
#include "phase_data.h"
#include "phase_models.h"

#define TDAT(k) (tdat[(k)-1])
#define NU(k)   (nu[(k)-1])

struct GPHWorkSpaceLeftRight {
  std::vector<std::vector<double>> barvf;
  std::vector<std::vector<double>> barvb;
  std::vector<std::vector<double>> vb;
  std::vector<std::vector<double>> vc;
  
  GPHWorkSpaceLeftRight(int m, int n) :
    barvf(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    barvb(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    vb(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    vc(std::vector<std::vector<double>>(m+1, std::vector<double>(n)))
  {}
};

template <typename T0, typename T1, typename T2, typename T4,
          typename T5, typename OptionT>
double llf(
    const GPH<T1,T2,T0>& model,
    const PHLeftRightSample<T4,T5>& data,
    OptionT& poisson_eps) noexcept {
  
  const int m = data.size();
  const double* tdat = stride_vector_traits<T4,double>::value(data.time);
  const int* nu = stride_vector_traits<T5,int>::value(data.nu);
  const double tmax = data.maxtime;
  
  int n = model.size();
  double qv = model.qv;
  std::vector<double> vone(n, 1.0);

  // work
  int right = poi::rightbound(qv*tmax, poisson_eps) + 1;
  std::vector<double> prob(right+1);
  std::vector<std::vector<double>> vx(right+1, std::vector<double>(n, 0));
  std::vector<double> tmpvb(n);
  std::vector<double> tmpv(n);
  std::vector<double> barvb0(n);
  std::vector<double> barvb1(n);
  
  double llf = 0.0;

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
    
    // vb[k] = (-ph.T) * barvb[k]
    gemv(NOTRANS{}, -1.0, model.Q, barvb1, 0.0, tmpvb);
    
    if (NU(k) == 3) { // left truncation time
      double tmp = dot(model.alpha, barvb1);
      llf -= log(tmp);
    } else if (NU(k) == 1) { // right censoring time
      double tmp = dot(model.alpha, barvb1);
      llf += log(tmp);
    } else if (NU(k) == 0) { // observed failure time
      double tmp = dot(model.alpha, tmpvb);
      llf += log(tmp);
    }
    copy(barvb1, barvb0);
  }
  
  return llf;
}

template <typename T0, typename T1, typename T2, typename T4,
          typename T5, typename T7, typename T8,
          typename OptionT, typename WorkSpace>
double estep(
    const GPH<T1,T2,T0>& model,
    const PHLeftRightSample<T4,T5>& data,
    GPHEres<T7,T8>& eres,
    OptionT& options,
    WorkSpace& work) noexcept {
  
  const int m = data.size();
  const double* tdat = stride_vector_traits<T4,double>::value(data.time);
  const int* nu = stride_vector_traits<T5,int>::value(data.nu);
  const double tmax = data.maxtime;

  int n = model.size();
  double qv = model.qv;
  std::vector<double> baralpha(n);
  std::vector<double> vone(n, 1.0);
  gesv(TRANS{}, -1.0, model.Q, model.alpha, baralpha);

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
  std::vector<double> wb(m+2);

  std::vector<std::vector<double>>& barvf(work.barvf);
  std::vector<std::vector<double>>& barvb(work.barvb);
  std::vector<std::vector<double>>& vb(work.vb);
  std::vector<std::vector<double>>& vc(work.vc);

  double llf = 0.0;
  double nn = 0;

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

    if (NU(k) == 3) { // left truncation time
      double tmp = dot(model.alpha, barvb[k]);
      llf -= log(tmp);
      wb[k] = 1.0 / tmp;
      nn += wb[k];
      axpy(wb[k], vone, eres.eb);
      axpy(-wb[k], barvb[k], eres.eb);
      axpy(wb[k], baralpha, eres.ey);
      axpy(-wb[k], barvf[k], eres.ey);
    } else if (NU(k) == 1) { // right censoring time
      double tmp = dot(model.alpha, barvb[k]);
      llf += log(tmp);
      wb[k] = 1.0 / tmp;
      axpy(wb[k], barvb[k], eres.eb);
      axpy(wb[k], barvf[k], eres.ey);
    } else if (NU(k) == 0) { // observed failure time
      double tmp = dot(model.alpha, vb[k]);
      llf += log(tmp);
      wb[k] = 1 / tmp;
      axpy(wb[k], vb[k], eres.eb);
      gemv(TRANS{}, -wb[k], model.Q, barvf[k], 1.0, eres.ey);
    }
  }

  // compute vectors and convolution

  fill(vc[m], 0.0);
  if (NU(m) == 3) {
    axpy(-wb[m], baralpha, vc[m]);
  } else if (NU(m) == 1) {
    axpy(wb[m], baralpha, vc[m]);
  } else if (NU(m) == 0) {
    axpy(wb[m], model.alpha, vc[m]);
  }
  for (int k=m-1; k>=1; k--) {
    // vc[k] = vc[k+1] * exp(T * tdat[k+1]) + ...
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
    if (NU(k) == 3) {
      axpy(-wb[k], baralpha, vc[k]);
    } else if (NU(k) == 1) {
      axpy(wb[k], baralpha, vc[k]);
    } else if (NU(k) == 0) {
      axpy(wb[k], model.alpha, vc[k]);
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
    if (NU(k) == 3) {
      ger(NOTRANS{}, wb[k], baralpha, vone, eres.en);
      ger(NOTRANS{}, -wb[k], baralpha, barvb[k], eres.en);
    } else if (NU(k) == 1) {
      ger(NOTRANS{}, wb[k], baralpha, barvb[k], eres.en);
    }
  }

  const double* alpha = vector_traits<T1>::value(model.alpha);
  const double* Q = vector_traits<T2>::value(model.Q);
  const double* xi = vector_traits<T1>::value(model.xi);
  const int* diag = vector_traits<T0,int>::value(model.diag);
  double* eb = vector_traits<T7>::value(eres.eb);
  double* ey = vector_traits<T7>::value(eres.ey);
  double* ez = vector_traits<T7>::value(eres.ez);
  double* en = vector_traits<T8>::value(eres.en);

  eres.etotal = nn;
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

