// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// const Vector & implementation ??
template <int RTYPE>
LogicalVector fvaryingCppImpl(Vector<RTYPE> x, int ng, IntegerVector g, bool any_group) { // int na_option,

  int l = x.size();
  typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
  auto isnanT = (RTYPE == REALSXP) ? [](typename Rcpp::traits::storage_type<RTYPE>::type x) { return x != x; } :
    [](typename Rcpp::traits::storage_type<RTYPE>::type x) { return x == Vector<RTYPE>::get_na(); };

    if(ng == 0) { // Note: Does not return NA if all NA... can be checked with fNobs ...
      int j = l-1;
      storage_t vi = x[j];
      while(isnanT(vi) && j!=0) vi = x[--j]; // think better about na option !! // if(na_option == 1)
      if(j != 0) for(int i = j; i--; ) if(!isnanT(x[i]) && x[i] != vi) return LogicalVector(1, true);
      return LogicalVector(1);
    } else { // with groups
      Vector<RTYPE> valg(ng, Vector<RTYPE>::get_na());
      if(any_group) {
        for(int i = 0; i != l; ++i) {
          if(isnanT(x[i])) continue; // Fastest ?? // think about na option !!
          if(isnanT(valg[g[i]-1])) {
            valg[g[i]-1] = x[i];
          } else { // slightly better than else if without brackets !!
            if(x[i] != valg[g[i]-1]) return LogicalVector(1, true);
          }
        }
        return LogicalVector(1);
      } else {
        LogicalVector varyg(ng, NA_LOGICAL);  // make optional .... varyg(ng);
        // int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(isnanT(x[i])) continue; // Fastest ?? // think about na option !!
          int gi = g[i]-1; // slightly faster !!
          if(isnanT(valg[gi])) {
            valg[gi] = x[i];
            varyg[gi] = false;
          } else {
            if(!varyg[gi] && x[i] != valg[gi]) { // !varyg[gi] && // best ??
              varyg[gi] = true; // fastest ?? !varyg[g[i]-1] && x[i] != valg[g[i]-1]
              // ++ngs; // Omitting this is faster for most datasets -> most are ordered !! (i.e. PRIO Grid 1.27 vs. 1.14 seconds)
              // if(ngs == ng) break;
            }
          }
        }
        Rf_setAttrib(varyg, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
        return varyg;
      }
    }
}

template <>
LogicalVector fvaryingCppImpl(Vector<CPLXSXP> x, int ng, IntegerVector g, bool any_group) { // , int na_option
  stop("Not supported SEXP type!");
}

template <>
LogicalVector fvaryingCppImpl(Vector<VECSXP> x, int ng, IntegerVector g, bool any_group) { // , int na_option
  stop("Not supported SEXP type!");
}

template <>
LogicalVector fvaryingCppImpl(Vector<RAWSXP> x, int ng, IntegerVector g, bool any_group) { // , int na_option
  stop("Not supported SEXP type!");
}

template <>
LogicalVector fvaryingCppImpl(Vector<EXPRSXP> x, int ng, IntegerVector g, bool any_group) { // , int na_option
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
LogicalVector fvaryingCpp(const SEXP& x, int ng = 0, const IntegerVector& g = 0, bool any_group = true){ // , int na_option = 1
  RCPP_RETURN_VECTOR(fvaryingCppImpl, x, ng, g, any_group);
}



template <int RTYPE>
SEXP fvaryingmCppImpl(Matrix<RTYPE> x, int ng, IntegerVector g, bool any_group, bool drop) { // , int na_option = 1
  int col = x.ncol();
  LogicalMatrix out = (ng == 0 || any_group) ? no_init_matrix(1, col) : no_init_matrix(ng, col);
  for(int j = col; j--; ) out(_, j) = fvaryingCppImpl<RTYPE>(x(_, j), ng, g, any_group);
  if(drop) {
    Rf_setAttrib(out, R_DimSymbol, R_NilValue);
    Rf_setAttrib(out, R_NamesSymbol, colnames(x));
  } else {
    colnames(out) = colnames(x);
  }
  return out;
}

template <>
SEXP fvaryingmCppImpl(Matrix<CPLXSXP> x, int ng, IntegerVector g, bool any_group, bool drop) { // , int na_option
  stop("Not supported SEXP type!");
}

template <>
SEXP fvaryingmCppImpl(Matrix<VECSXP> x, int ng, IntegerVector g, bool any_group, bool drop) { // , int na_option
  stop("Not supported SEXP type!");
}

template <>
SEXP fvaryingmCppImpl(Matrix<RAWSXP> x, int ng, IntegerVector g, bool any_group, bool drop) { // , int na_option
  stop("Not supported SEXP type!");
}

template <>
SEXP fvaryingmCppImpl(Matrix<EXPRSXP> x, int ng, IntegerVector g, bool any_group, bool drop) { // , int na_option
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP fvaryingmCpp(const SEXP& x, int ng = 0, const IntegerVector& g = 0, bool any_group = true, bool drop = true){ // , int na_option = 1
  RCPP_RETURN_MATRIX(fvaryingmCppImpl, x, ng, g, any_group, drop);
}


// [[Rcpp::export]]
SEXP fvaryinglCpp(const List& x, int ng = 0, const IntegerVector& g = 0, bool any_group = true, bool drop = true) { // , int na_option = 1
  int l = x.size();
  List out(l);
  for(int j = l; j--; ) {
    switch(TYPEOF(x[j])) {
    case REALSXP:
      out[j] = fvaryingCppImpl<REALSXP>(x[j], ng, g, any_group);
      break;
    case INTSXP:
      out[j] = fvaryingCppImpl<INTSXP>(x[j], ng, g, any_group);
      break;
    case STRSXP:
      out[j] = fvaryingCppImpl<STRSXP>(x[j], ng, g, any_group);
      break;
    case LGLSXP:
      out[j] = fvaryingCppImpl<LGLSXP>(x[j], ng, g, any_group);
      break;
    default: stop("Not supported SEXP type !");
    }
  }
  if(drop) {
    LogicalVector outl = no_init_vector(l);
    for(int i = l; i--; ) outl[i] = out[i];
    Rf_setAttrib(outl, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
    return outl;
  } else {
    DUPLICATE_ATTRIB(out, x);
    if(ng == 0 || any_group) out.attr("row.names") = 1;
    else out.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng);
    return out;
  }
}
