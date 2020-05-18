// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// const Vector & implementation ? na_option ?
template <int RTYPE>
LogicalVector varyingCppImpl(Vector<RTYPE> x, int ng, IntegerVector g, bool any_group) {

  int l = x.size();
  typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
  auto isnanT = (RTYPE == REALSXP) ? [](storage_t x) { return x != x; } :
    [](storage_t x) { return x == Vector<RTYPE>::get_na(); };

    if(ng == 0) { // Note: Does not return NA if all NA... can be checked with fNobs ...
      int j = l-1;
      storage_t vi = x[j];
      while(isnanT(vi) && j!=0) vi = x[--j];
      if(j != 0) for(int i = j; i--; ) if(!isnanT(x[i]) && x[i] != vi) return LogicalVector(1, true);
      return LogicalVector(1);
    } else { // with groups
      Vector<RTYPE> valg(ng, Vector<RTYPE>::get_na());
      if(any_group) {
        for(int i = 0; i != l; ++i) {
          if(isnanT(x[i])) continue; // Fastest ?
          if(isnanT(valg[g[i]-1])) {
            valg[g[i]-1] = x[i];
          } else { // slightly better than else if without brackets
            if(x[i] != valg[g[i]-1]) return LogicalVector(1, true);
          }
        }
        return LogicalVector(1);
      } else {
        LogicalVector varyg(ng, NA_LOGICAL);  // make optional .... varyg(ng);
        // int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(isnanT(x[i])) continue; // Fastest ?
          int gi = g[i]-1; // slightly faster
          if(isnanT(valg[gi])) {
            valg[gi] = x[i];
            varyg[gi] = false;
          } else {
            if(!varyg[gi] && x[i] != valg[gi]) { // best ?
              varyg[gi] = true; // fastest ?
              // ++ngs; // Omitting this is faster for most datasets -> most are ordered ! (i.e. PRIO Grid 1.27 vs. 1.14 seconds)
              // if(ngs == ng) break;
            }
          }
        }
        // Rf_setAttrib(varyg, R_NamesSymbol, R_NilValue);
        return varyg;
      }
    }
}

template <>
LogicalVector varyingCppImpl(Vector<CPLXSXP> x, int ng, IntegerVector g, bool any_group) {
  stop("Not supported SEXP type!");
}

template <>
LogicalVector varyingCppImpl(Vector<VECSXP> x, int ng, IntegerVector g, bool any_group) {
  stop("Not supported SEXP type!");
}

template <>
LogicalVector varyingCppImpl(Vector<RAWSXP> x, int ng, IntegerVector g, bool any_group) {
  stop("Not supported SEXP type!");
}

template <>
LogicalVector varyingCppImpl(Vector<EXPRSXP> x, int ng, IntegerVector g, bool any_group) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
LogicalVector varyingCpp(const SEXP& x, int ng = 0, const IntegerVector& g = 0, bool any_group = true){
  RCPP_RETURN_VECTOR(varyingCppImpl, x, ng, g, any_group);
}



template <int RTYPE>
SEXP varyingmCppImpl(Matrix<RTYPE> x, int ng, IntegerVector g, bool any_group, bool drop) {
  int col = x.ncol();
  LogicalMatrix out = (ng == 0 || any_group) ? no_init_matrix(1, col) : no_init_matrix(ng, col);
  for(int j = col; j--; ) out(_, j) = varyingCppImpl<RTYPE>(x(_, j), ng, g, any_group);
  if(drop && any_group) {
    Rf_setAttrib(out, R_DimSymbol, R_NilValue);
    Rf_setAttrib(out, R_NamesSymbol, colnames(x));
  } else {
    colnames(out) = colnames(x);
  }
  return out;
}

template <>
SEXP varyingmCppImpl(Matrix<CPLXSXP> x, int ng, IntegerVector g, bool any_group, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP varyingmCppImpl(Matrix<VECSXP> x, int ng, IntegerVector g, bool any_group, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP varyingmCppImpl(Matrix<RAWSXP> x, int ng, IntegerVector g, bool any_group, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP varyingmCppImpl(Matrix<EXPRSXP> x, int ng, IntegerVector g, bool any_group, bool drop) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP varyingmCpp(const SEXP& x, int ng = 0, const IntegerVector& g = 0, bool any_group = true, bool drop = true){
  RCPP_RETURN_MATRIX(varyingmCppImpl, x, ng, g, any_group, drop);
}


// [[Rcpp::export]]
SEXP varyinglCpp(const List& x, int ng = 0, const IntegerVector& g = 0, bool any_group = true, bool drop = true) {
  int l = x.size();
  List out(l);
  for(int j = l; j--; ) {
    switch(TYPEOF(x[j])) {
    case REALSXP:
      out[j] = varyingCppImpl<REALSXP>(x[j], ng, g, any_group);
      break;
    case INTSXP:
      out[j] = varyingCppImpl<INTSXP>(x[j], ng, g, any_group);
      break;
    case STRSXP:
      out[j] = varyingCppImpl<STRSXP>(x[j], ng, g, any_group);
      break;
    case LGLSXP:
      out[j] = varyingCppImpl<LGLSXP>(x[j], ng, g, any_group);
      break;
    default: stop("Not supported SEXP type !");
    }
  }
  if(drop && any_group) {
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
