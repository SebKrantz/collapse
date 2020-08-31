// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>
LogicalVector ganyNACppImpl(Vector<RTYPE> x, int ng, IntegerVector g) {
  int l = x.size();
  if(l != g.size()) stop("length(x) must match length(g)");
  LogicalVector out(ng);

  if(RTYPE == REALSXP) { // numeric vector: all logical operations on NA/NaN evaluate to false, except != which is true.
    for(int i = 0; i < l; ++i) if(x[i] != x[i]) out[g[i]-1] = true;
  } else { // other vectors
    for(int i = 0; i < l; ++i) if(x[i] == Vector<RTYPE>::get_na()) out[g[i]-1] = true;
  }

  return out;
}

// disabling other types
template <>
LogicalVector ganyNACppImpl(Vector<CPLXSXP> x, int ng, IntegerVector) {
  stop("Not supported SEXP type!");
}

template <>
LogicalVector ganyNACppImpl(Vector<VECSXP> x, int ng, IntegerVector) {
  stop("Not supported SEXP type!");
}

template <>
LogicalVector ganyNACppImpl(Vector<RAWSXP> x, int ng, IntegerVector) {
  stop("Not supported SEXP type!");
}

template <>
LogicalVector ganyNACppImpl(Vector<EXPRSXP> x, int ng, IntegerVector) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
LogicalVector ganyNACpp(const SEXP& x, int ng = 0, const IntegerVector& g = 0){
  RCPP_RETURN_VECTOR(ganyNACppImpl, x, ng, g);
}
