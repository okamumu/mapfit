#ifndef MAPFIT_ESTEP_PHASE_WTIME_GPH
#define MAPFIT_ESTEP_PHASE_WTIME_GPH

#include <Rcpp.h>
#include "traits.h"
#include "poisson.h"
#include "gamma.h"
#include "blas.h"
#include "unif.h"
#include "phase_data.h"
#include "phase_models.h"

#define TDAT(k) (tdat[(k)-1])
#define WDAT(k) (wdat[(k)-1])

struct GPHWorkSpaceWTime {
  std::vector<std::vector<double>> vf;
  std::vector<std::vector<double>> vb;
  std::vector<std::vector<double>> vc;
  
  GPHWorkSpaceWTime(int m, int n) :
    vf(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    vb(std::vector<std::vector<double>>(m+1, std::vector<double>(n))),
    vc(std::vector<std::vector<double>>(m+1, std::vector<double>(n)))
  {}
};

template <typename T0, typename T1, typename T2, typename T4, typename T5,
          typename OptionT>
double llf(
    const GPH<T1,T2,T0>& model,
    const PHWeightSample<T4,T5>& data,
    OptionT& poisson_eps) noexcept {
  
  const int m = data.size();
  const double* tdat = vector_traits<T4,double>::value(data.time);
  const double* wdat = vector_traits<T5,double>::value(data.weights);
  double tmax = data.maxtime;
  
  int n = model.size();
  double qv = model.qv;
  
  // alloc
  int right = poi::rightbound(qv*tmax, poisson_eps) + 1;
  std::vector<double> prob(right+1);
  std::vector<double> tmpv(n);
  std::vector<double> vf0(n);
  std::vector<double> vf1(n);

  double llf = 0.0;

  // forward & backward
  copy(model.alpha, vf0);
  for (int k=1; k<=m; k++) {
    int right = poi::rightbound(qv*TDAT(k), poisson_eps) + 1;
    double weight = poi::pmf(qv*TDAT(k), 0, right, prob);
    
    // mexp::mexp_unifvec(sci::mat::T, P, qv, r, poi, weight, vf0, vf1, 0.0);
    fill(vf1, 0.0);
    axpy(prob[0], vf0, vf1);
    for (int u=1; u<=right; u++) {
      gemv(TRANS{}, 1.0, model.P, vf0, 0.0, tmpv);
      copy(tmpv, vf0);
      axpy(prob[u], vf0, vf1);
    }
    scal(1.0/weight, vf1);
    llf += WDAT(k) * log(dot(vf1, model.xi));
    copy(vf1, vf0);

  }
  
  return llf;
}

template <typename T0, typename T1, typename T2, typename T4, typename T5, typename T6, typename T7,
          typename OptionT, typename WorkSpace>
double estep(
    const GPH<T1,T2,T0>& model,
    const PHWeightSample<T4,T5>& data,
    GPHEres<T6,T7>& eres,
    OptionT& options,
    WorkSpace& work) noexcept {
  
  const int m = data.size();
  const double* tdat = vector_traits<T4,double>::value(data.time);
  const double* wdat = vector_traits<T5,double>::value(data.weights);
  double tmax = data.maxtime;
  
  int n = model.size();
  double qv = model.qv;
  
  // alloc
  int right = poi::rightbound(qv*tmax, options.poisson_eps) + 1;
  std::vector<double> prob(right+1);
  std::vector<std::vector<double>> vx(right+1, std::vector<double>(n));
  // std::vector<std::vector<double>> vf(m+1, std::vector<double>(n));
  // std::vector<std::vector<double>> vb(m+1, std::vector<double>(n));
  // std::vector<std::vector<double>> vc(m+1, std::vector<double>(n));
  std::vector<double> blf(m+1);
  std::vector<double> tmpv(n);
  std::vector<double> xtmp(n);

  std::vector<std::vector<double>>& vf(work.vf);
  std::vector<std::vector<double>>& vb(work.vb);
  std::vector<std::vector<double>>& vc(work.vc);
  
  double scale;
  double llf = 0.0;
  double tllf = 0.0;
  eres.etotal = 0.0;
  fill(eres.eb, 0.0);
  fill(eres.ey, 0.0);
  fill(eres.ez, 0.0);
  fill(eres.en, 0.0);
  
  // forward & backward
  copy(model.alpha, vf[0]);
  copy(model.xi, vb[0]);
  for (int k=1; k<=m; k++) {
    int right = poi::rightbound(qv*TDAT(k), options.poisson_eps) + 1;
    double weight = poi::pmf(qv*TDAT(k), 0, right, prob);
    
    // mexp::mexp_unifvec(sci::mat::T, P, qv, r, poi, weight, vf[k-1], vf[k], 0.0);
    fill(vf[k], 0.0);
    copy(vf[k-1], xtmp);
    axpy(prob[0], xtmp, vf[k]);
    for (int u=1; u<=right; u++) {
      gemv(TRANS{}, 1.0, model.P, xtmp, 0.0, tmpv);
      copy(tmpv, xtmp);
      axpy(prob[u], xtmp, vf[k]);
    }
    scal(1.0/weight, vf[k]);
    scale = dot(vf[k], model.xi);
    scal(1.0/scale, vf[k]);
    axpy(WDAT(k), vf[k], eres.ey);

    blf[k] = scale;
    // mexp::mexp_unifvec(sci::mat::N, P, qv, r, poi, weight, vb[k-1], vb[k], 0.0);
    fill(vb[k], 0.0);
    copy(vb[k-1], xtmp);
    axpy(prob[0], xtmp, vb[k]);
    for (int u=1; u<=right; u++) {
      gemv(NOTRANS{}, 1.0, model.P, xtmp, 0.0, tmpv);
      copy(tmpv, xtmp);
      axpy(prob[u], xtmp, vb[k]);
    }
    scal(1.0/weight, vb[k]);
    scale = dot(model.alpha, vb[k]);
    scal(1.0/scale, vb[k]);
    axpy(WDAT(k), vb[k], eres.eb);
    
    eres.etotal += WDAT(k);
    tllf += log(blf[k]);
    llf += WDAT(k) * tllf;
  }
  
  fill(vc[m], 0.0);
  axpy(WDAT(m)/blf[m], model.alpha, vc[m]);
  for (int k=m-1; k>=1; k--) {
    int right = poi::rightbound(qv*TDAT(k+1), options.poisson_eps) + 1;
    double weight = poi::pmf(qv*TDAT(k+1), 0, right, prob);
    
    // mexp::mexp_unifvec(sci::mat::T, P, qv, r, poi, weight, vc[k+1], vc[k], 0.0);
    fill(vc[k], 0.0);
    copy(vc[k+1], xtmp);
    axpy(prob[0], xtmp, vc[k]);
    for (int u=1; u<=right; u++) {
      gemv(TRANS{}, 1.0, model.P, xtmp, 0.0, tmpv);
      copy(tmpv, xtmp);
      axpy(prob[u], xtmp, vc[k]);
    }
    scal(1.0/(weight*blf[k]), vc[k]);
    axpy(WDAT(k)/blf[k], model.alpha, vc[k]);
  }
  
  for (int k=1; k<=m; k++) {
    int right = poi::rightbound(qv*TDAT(k), options.poisson_eps) + 1;
    double weight = poi::pmf(qv*TDAT(k), 0, right, prob);
    
    // mexp::mexpc_unif(sci::mat::T, sci::mat::N, P, qv, r, poi, weight,
    //                  vc[k], vb[k-1], vb[k-1], en);
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
  }
  
  const double* alpha = vector_traits<T1>::value(model.alpha);
  const double* Q = vector_traits<T2>::value(model.Q);
  const double* xi = vector_traits<T1>::value(model.xi);
  const int* diag = vector_traits<T0,int>::value(model.diag);
  double* eb = vector_traits<T6>::value(eres.eb);
  double* ey = vector_traits<T6>::value(eres.ey);
  double* ez = vector_traits<T6>::value(eres.ez);
  double* en = vector_traits<T7>::value(eres.en);
  
  for (int i=0; i<n; i++) {
    eb[i] *= alpha[i];
    ey[i] *= xi[i];
    ez[i] = en[diag[i]];
  }
  for (int i=0; i<vector_traits<T7>::size(eres.en); i++) {
    en[i] *= Q[i];
  }

  return llf;
}

#endif

