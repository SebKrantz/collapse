// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>
inline bool isnaNUM(typename Rcpp::traits::storage_type<RTYPE>::type x) {
  return x != x;
}

template <int RTYPE>
inline bool isnaOTH(typename Rcpp::traits::storage_type<RTYPE>::type x) {
  return x == Vector<RTYPE>::get_na();
}

template <int RTYPE>
Vector<RTYPE> na_rmImpl(const Vector<RTYPE>& x) { // better than logical vector !!!
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;
  int l = x.size(), count = 0;
  for(int i = 0; i != l; ++i) if(!isnanT(x[i])) ++count;
  Vector<RTYPE> y = no_init_vector(count);
  DUPLICATE_ATTRIB(y, x); // doesn't matter where it is, here or at the end !!
  int j = 0;
  SEXP nam = Rf_getAttrib(x, R_NamesSymbol);
  if(Rf_isNull(nam)) {
    for(int i = 0; i != count; ++i) if(!isnanT(x[i])) y[j++] = x[i];
  } else {
    CharacterVector name = nam;
    CharacterVector names = no_init_vector(count);
    for(int i = 0; i != count; ++i) {
      if(isnanT(x[i])) continue;
      y[j] = x[i];
      names[j++] = name[i];
    }
    Rf_setAttrib(y, R_NamesSymbol, names);
    // y.attr("names") = names;
  }
  // auto it = std::remove_copy_if(x.begin(), x.end(), y.begin(), isnanT);
  // y.erase(it, y.end());
  // SEXP nam = Rf_getAttrib(x, R_NamesSymbol);
  // if(nam != R_NilValue) {
  //   CharacterVector names = nam;
  //   out.attr("names") = CharacterVector::create(names[seq(1,y.size())]);
  // }
  return y;
}

template <>
Vector<CPLXSXP> na_rmImpl(const Vector<CPLXSXP>& x) {
  stop("Not supported SEXP type!");
}

template <>
Vector<VECSXP> na_rmImpl(const Vector<VECSXP>& x) {
  stop("Not supported SEXP type!");
}

template <>
Vector<RAWSXP> na_rmImpl(const Vector<RAWSXP>& x) {
  stop("Not supported SEXP type!");
}

template <>
Vector<EXPRSXP> na_rmImpl(const Vector<EXPRSXP>& x) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP na_rm(SEXP x) {
  RCPP_RETURN_VECTOR(na_rmImpl, x);
}

