// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>
SEXP fnobsmCppImpl(Matrix<RTYPE> x, int ng, IntegerVector g, bool drop) { 
  int l = x.nrow(), col = x.ncol(); 
  
  if(ng == 0) { // Fastest loops ?? Is perhaps an iterator metter suited?? or a STD algorithm to count ??
    IntegerVector nobs = no_init_vector(col);
      if(TYPEOF(x) == REALSXP) {
        for(int j = col; j--; ) { // fastest loop???????????????
          MatrixColumn<RTYPE> column = x( _ , j); 
          int ni = 0; // fastest way??????????????? 
          for(int i = 0; i != l; ++i) if(column[i] == column[i]) ++ni;
          nobs[j] = ni;
        }
      } else {
        for(int j = col; j--; ) {
          MatrixColumn<RTYPE> column = x( _ , j); 
          int ni = 0;
          for(int i = 0; i != l; ++i) if(column[i] != Vector<RTYPE>::get_na()) ++ni;
          nobs[j] = ni;
        }
      }
    if(drop) nobs.attr("names") = colnames(x); 
    else {
      nobs.attr("dim") = Dimension(1, col);
      colnames(nobs) = colnames(x); 
    }
    return nobs;
  } else { // with groups 
    if(g.size() != l) stop("length(g) must match nrow(X)");
    IntegerMatrix nobs(ng, col); // init best solution ??
      if(TYPEOF(x) == REALSXP) {
        for(int j = col; j--; ) { 
          MatrixColumn<RTYPE> column = x( _ , j); 
          IntegerMatrix::Column nj = nobs( _ , j); 
          for(int i = 0; i != l; ++i) if(column[i] == column[i]) ++nj[g[i]-1];
        }
      } else {
        for(int j = col; j--; ) { 
          MatrixColumn<RTYPE> column = x( _ , j); 
          IntegerMatrix::Column nj = nobs( _ , j); 
          for(int i = 0; i != l; ++i) if(column[i] != Vector<RTYPE>::get_na()) ++nj[g[i]-1];
        }
      }
    colnames(nobs) = colnames(x);
    return nobs;
  }
}

template <>
SEXP fnobsmCppImpl(Matrix<CPLXSXP> x, int ng, IntegerVector g, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fnobsmCppImpl(Matrix<VECSXP> x, int ng, IntegerVector g, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fnobsmCppImpl(Matrix<RAWSXP> x, int ng, IntegerVector g, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fnobsmCppImpl(Matrix<EXPRSXP> x, int ng, IntegerVector g, bool drop) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP fnobsmCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool drop = true){
  RCPP_RETURN_MATRIX(fnobsmCppImpl, x, ng, g, drop);
}

