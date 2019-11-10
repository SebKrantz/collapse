// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Rather call it fN  or fnobs ??
template <int RTYPE>
IntegerVector fnobsCppImpl(Vector<RTYPE> x, int ng, IntegerVector g) {
  
  int l = x.size();
  
  if (ng == 0) {
      int n = 0;
      if(TYPEOF(x) == REALSXP) {
        for(int i = 0; i != l; ++i) if(x[i] == x[i]) ++n; // This loop is faster !!
      } else {
        for(int i = 0; i != l; ++i) if(x[i] != Vector<RTYPE>::get_na()) ++n;
      }
     return IntegerVector::create(n);
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    IntegerVector n(ng); 
      if(TYPEOF(x) == REALSXP) { 
        for(int i = 0; i != l; ++i) if(x[i] == x[i]) ++n[g[i]-1];
      } else {
        for(int i = 0; i != l; ++i) if(x[i] != Vector<RTYPE>::get_na()) ++n[g[i]-1];
      }
    if(Rf_getAttrib(x, R_ClassSymbol) == R_NilValue) {
      SHALLOW_DUPLICATE_ATTRIB(n, x); 
    } else {
      n.attr("label") = x.attr("label");
    }
    return n;
  }
}

template <>
IntegerVector fnobsCppImpl(Vector<CPLXSXP> x, int ng, IntegerVector) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector fnobsCppImpl(Vector<VECSXP> x, int ng, IntegerVector) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector fnobsCppImpl(Vector<RAWSXP> x, int ng, IntegerVector) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector fnobsCppImpl(Vector<EXPRSXP> x, int ng, IntegerVector) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
IntegerVector fnobsCpp(SEXP x, int ng = 0, IntegerVector g = 0){
  RCPP_RETURN_VECTOR(fnobsCppImpl, x, ng, g);
}
