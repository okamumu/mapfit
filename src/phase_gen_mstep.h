#ifndef MAPFIT_MSTEP_PHASE_GPH
#define MAPFIT_MSTEP_PHASE_GPH

#include <Rcpp.h>
#include <vector>

#include "traits.h"
#include "unif.h"
#include "phase_data.h"
#include "phase_models.h"

namespace _mstep_ {

template <typename IVecT, typename EresVecT, typename EresMatT, typename VecT, typename MatT,
          typename OptionT>
void mstep(const GPHEres<EresVecT,EresMatT>& eres,
           GPH<VecT,MatT,IVecT>& model,
           OptionT& options, DenseMatrixT) noexcept{
  const int n = model.size();
  const double* eb = vector_traits<EresVecT>::value(eres.eb);
  const double* ey = vector_traits<EresVecT>::value(eres.ey);
  const double* ez = vector_traits<EresVecT>::value(eres.ez);
  const double* en = dense_matrix_traits<EresMatT>::value(eres.en);
  const int ld = dense_matrix_traits<EresMatT>::ld(eres.en);
  double* alpha = vector_traits<VecT>::value(model.alpha);
  double* xi = vector_traits<VecT>::value(model.xi);
  double* Q = dense_matrix_traits<MatT>::value(model.Q);
  const int ldq = dense_matrix_traits<MatT>::ld(model.Q);
  const int* d = vector_traits<IVecT,int>::value(model.diag);
  std::vector<double> tmpv(n, 0.0);
  
  for (int j=0; j<n; j++) {
    for (int i=0; i<n; i++) {
      if (i != j) {
        Q[j*ldq+i] = en[j*ld+i] / ez[i];
        tmpv[i] += Q[j*ldq+i];
      }
    }
  }
  for (int i=0; i<n; i++) {
    alpha[i] = eb[i] / eres.etotal;
    xi[i] = ey[i] / ez[i];
    tmpv[i] += xi[i];
    Q[d[i]] = -tmpv[i];
  }
}

template <typename IVecT, typename EresVecT, typename EresMatT, typename VecT, typename MatT,
          typename OptionT>
void mstep(const GPHEres<EresVecT,EresMatT>& eres,
           GPH<VecT,MatT,IVecT>& model,
           OptionT& options, CSRMatrixT) {
  const int n = model.size();
  const double* eb = vector_traits<EresVecT>::value(eres.eb);
  const double* ey = vector_traits<EresVecT>::value(eres.ey);
  const double* ez = vector_traits<EresVecT>::value(eres.ez);
  const double* en = csr_matrix_traits<EresMatT>::value(eres.en);
  
  double* alpha = vector_traits<VecT>::value(model.alpha);
  double* xi = vector_traits<VecT>::value(model.xi);
  double* Q = csr_matrix_traits<MatT>::value(model.Q);
  const int* rowptr = csr_matrix_traits<MatT>::rowptr(model.Q);
  const int* colind = csr_matrix_traits<MatT>::colind(model.Q);
  const int base = csr_matrix_traits<MatT>::base(model.Q);
  const int* d = vector_traits<IVecT,int>::value(model.diag);
  std::vector<double> tmpv(n, 0.0);
  
  for (int i=0; i<n; i++) {
    for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
      int j = colind[z] - base;
      if (i != j) {
        Q[z] = en[z] / ez[i];
        tmpv[i] += Q[z];
      }
    }
  }
  for (int i=0; i<n; i++) {
    alpha[i] = eb[i] / eres.etotal;
    xi[i] = ey[i] / ez[i];
    tmpv[i] += xi[i];
    Q[d[i]] = -tmpv[i];
  }
}

template <typename IVecT, typename EresVecT, typename EresMatT, typename VecT, typename MatT,
          typename OptionT>
void mstep(const GPHEres<EresVecT,EresMatT>& eres,
           GPH<VecT,MatT,IVecT>& model,
           OptionT& options, CSCMatrixT) {
  const int n = model.size();
  const double* eb = vector_traits<EresVecT>::value(eres.eb);
  const double* ey = vector_traits<EresVecT>::value(eres.ey);
  const double* ez = vector_traits<EresVecT>::value(eres.ez);
  const double* en = csc_matrix_traits<EresMatT>::value(eres.en);

  double* alpha = vector_traits<VecT>::value(model.alpha);
  double* xi = vector_traits<VecT>::value(model.xi);
  double* Q = csc_matrix_traits<MatT>::value(model.Q);
  const int* colptr = csc_matrix_traits<MatT>::colptr(model.Q);
  const int* rowind = csc_matrix_traits<MatT>::rowind(model.Q);
  const int base = csc_matrix_traits<MatT>::base(model.Q);
  const int* d = vector_traits<IVecT,int>::value(model.diag);
  std::vector<double> tmpv(n, 0.0);
  
  for (int j=0; j<n; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      if (i != j) {
        Q[z] = en[z] / ez[i];
        tmpv[i] += Q[z];
      }
    }
  }
  for (int i=0; i<n; i++) {
    alpha[i] = eb[i] / eres.etotal;
    xi[i] = ey[i] / ez[i];
    tmpv[i] += xi[i];
    Q[d[i]] = -tmpv[i];
  }
}

template <typename IVecT, typename EresVecT, typename EresMatT, typename VecT, typename MatT,
          typename OptionT>
void mstep(const GPHEres<EresVecT,EresMatT>& eres,
           GPH<VecT,MatT,IVecT>& model,
           OptionT& options, COOMatrixT) {
  const int n = model.size();
  const double* eb = vector_traits<EresVecT>::value(eres.eb);
  const double* ey = vector_traits<EresVecT>::value(eres.ey);
  const double* ez = vector_traits<EresVecT>::value(eres.ez);
  const double* en = coo_matrix_traits<EresMatT>::value(eres.en);
  
  double* alpha = vector_traits<VecT>::value(model.alpha);
  double* xi = vector_traits<VecT>::value(model.xi);
  double* Q = coo_matrix_traits<MatT>::value(model.Q);
  const int* rowind = coo_matrix_traits<MatT>::rowind(model.Q);
  const int* colind = coo_matrix_traits<MatT>::colind(model.Q);
  const int base = coo_matrix_traits<MatT>::base(model.Q);
  const int nnz = coo_matrix_traits<MatT>::nnz(model.Q);
  const int* d = vector_traits<IVecT,int>::value(model.diag);
  std::vector<double> tmpv(n, 0.0);
  
  for (int z=0; z<nnz; z++) {
    int i = rowind[z] - base;
    int j = colind[z] - base;
    if (i != j) {
      Q[z] = en[z] / ez[i];
      tmpv[i] += Q[z];
    }
  }
  for (int i=0; i<n; i++) {
    alpha[i] = eb[i] / eres.etotal;
    xi[i] = ey[i] / ez[i];
    tmpv[i] += xi[i];
    Q[d[i]] = -tmpv[i];
  }
}

}

template <typename IVecT, typename EresVecT, typename EresMatT, typename VecT, typename MatT,
          typename OptionT>
void mstep(const GPHEres<EresVecT,EresMatT>& eres,
           GPH<VecT,MatT,IVecT>& model,
           OptionT& options) {
  _mstep_::mstep(eres, model, options, typename matrix_category<MatT>::type{});
  // unif
  copy(model.Q, model.P);
  double qv = unif(model.P, model.diag, options.ufactor);
  model.qv = qv;
}

template <typename EresVecT, typename EresMatT, typename GPHT, typename OptionT>
void mstep(const GPHEres<EresVecT,EresMatT>& eres,
           GPHPoi<GPHT>& model,
           OptionT& options) noexcept {
  mstep(eres, model.gph, options);
  model.omega = eres.etotal;
}

#endif

