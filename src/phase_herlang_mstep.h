#ifndef MAPFIT_PHASE_ERLANG_MSTEP_H
#define MAPFIT_PHASE_ERLANG_MSTEP_H

#include <Rcpp.h>
#include "traits.h"
#include "blas.h"
#include "phase_data.h"
#include "phase_models.h"

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
