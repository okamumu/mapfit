#ifndef MAPFIT_MODEL_TRAITS_H
#define MAPFIT_MODEL_TRAITS_H

#include "traits.h"

template <typename VectorT, typename IVectorT>
struct HErlang {
  VectorT alpha;
  IVectorT shape;
  VectorT rate;
  
  HErlang(const VectorT& _alpha,
          const IVectorT& _shape,
          const VectorT& _rate)
    : alpha(_alpha), shape(_shape), rate(_rate) {}

  inline int size() const {
    return vector_traits<VectorT>::size(alpha);
  }
};

template <typename VectorT, typename IVectorT>
struct HErlangPoi {
  VectorT alpha;
  IVectorT shape;
  VectorT rate;
  double omega;
  
  HErlangPoi(const VectorT& _alpha,
          const IVectorT& _shape,
          const VectorT& _rate,
          double _omega)
    : alpha(_alpha), shape(_shape), rate(_rate), omega(_omega) {}
  
  inline int size() const {
    return vector_traits<VectorT>::size(alpha);
  }
};

template <typename VectorT>
struct HErlangEres {
  double etotal;
  VectorT eb;
  VectorT ew;
  
  inline
    HErlangEres(const VectorT& _eb, const VectorT& _ew) : etotal(0), eb(_eb), ew(_ew) {}
};

template <typename VectorT, typename MatrixT, typename IVectorT>
struct GPH {
  VectorT alpha;
  MatrixT Q;
  MatrixT P;
  VectorT xi;
  double qv;
  IVectorT diag;

  GPH(const VectorT& _alpha,
      const MatrixT& _Q,
      const MatrixT& _P,
      const VectorT& _xi,
      double _qv,
      const IVectorT& _diag)
    : alpha(_alpha), Q(_Q), P(_P), xi(_xi), qv(_qv), diag(_diag) {}
  
  inline int size() const {
    return vector_traits<VectorT>::size(alpha);
  }
};

template <typename GPHT>
struct GPHPoi {
  GPHT gph;
  double omega;
  
  GPHPoi(const GPHT& _gph, double _omega)
    : gph(_gph), omega(_omega) {}
};

template <typename VectorT, typename GPHT>
struct CF1 {
  VectorT alpha;
  VectorT rate;
  GPHT gph;

  CF1(const VectorT& _alpha, const VectorT& _rate, const GPHT& _gph)
    : alpha(_alpha), rate(_rate), gph(_gph) {}
};

template <typename VectorT, typename MatrixT>
struct GPHEres {
  double etotal;
  VectorT eb;
  VectorT ey;
  VectorT ez;
  MatrixT en;
  
  inline
    GPHEres(const VectorT& _eb, const VectorT& _ey, const VectorT& _ez, const MatrixT& _en)
    : etotal(0),
      eb(_eb),
      ey(_ey),
      ez(_ez),
      en(_en) {}
};

#endif