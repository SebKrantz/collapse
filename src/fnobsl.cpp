#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fnobslCpp(List x, int ng = 0, IntegerVector g = 0, bool drop = true) {
  int l = x.size();
  
  if (ng == 0) {
    NumericVector nobs = no_init_vector(l); 
      for(int j = l; j--; ) { // fastest loop ???
      // for(int j = 0; j != l; ++j) { Not sure, could be faster !!!
        int ni = 0;
        switch(TYPEOF(x[j])) { // Faster than using iterator ?? // https://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
        case REALSXP: {
          NumericVector column = x[j];
          int k = column.size();
          for(int i = 0; i != k; ++i) if(!std::isnan(column[i])) ++ni;
          // for(int i = 0; i != column.size(); ++i) if(!std::isnan(column[i])) ++ni; // Note: Column size function called repeatedly is very slow!!!
          break;    
        }
        case INTSXP: {
          IntegerVector column = x[j];
          int k = column.size();
          for(int i = 0; i != k; ++i) if(column[i] != NA_INTEGER) ++ni; 
          break;    
        }
        case STRSXP: {
          CharacterVector column = x[j];
          int k = column.size();
          for(int i = 0; i != k; ++i) if(column[i] != NA_STRING) ++ni;
          break;    
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          int k = column.size();
          for(int i = 0; i != k; ++i) if(column[i] != NA_LOGICAL) ++ni;
          break;    
        }
        default: {
          stop("incompatible SEXP encountered;");
          break;
        }
        }
        nobs[j] = ni;
      }
    if(drop) { 
      nobs.attr("names") = x.attr("names");
      return nobs;
    } else {
      List out(l);
      for(int j = l; j--; ) {
        out[j] = nobs[j];
        if(Rf_getAttrib(x[j], R_ClassSymbol) == R_NilValue) {
          SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
        } else {
          Rf_setAttrib(out[j], wrap("label"), Rf_getAttrib(x[j], wrap("label")));
        }
      }
      DUPLICATE_ATTRIB(out, x);
      out.attr("row.names") = 1;
      return out;
    }
  } else { // With groups !!
    List nobs(l); 
    for(int j = l; j--; ) { // fastest loop ???
      IntegerVector ni(ng);
      switch(TYPEOF(x[j])) { 
      case REALSXP: {
        NumericVector column = x[j];
        int k = column.size();
        for(int i = 0; i != k; ++i) if(!std::isnan(column[i])) ++ni[g[i]-1];
        break;    
      }
      case INTSXP: {
        IntegerVector column = x[j];
        int k = column.size();
        for(int i = 0; i != k; ++i) if(column[i] != NA_INTEGER) ++ni[g[i]-1]; 
        break;    
      }
      case STRSXP: {
        CharacterVector column = x[j];
        int k = column.size();
        for(int i = 0; i != k; ++i) if(column[i] != NA_STRING) ++ni[g[i]-1];
        break;    
      }
      case LGLSXP: {
        LogicalVector column = x[j];
        int k = column.size();
        for(int i = 0; i != k; ++i) if(column[i] != NA_LOGICAL) ++ni[g[i]-1];
        break;    
      }
      default: {
        stop("incompatible SEXP encountered;");
        break;
      }
      }
      if(Rf_getAttrib(x[j], R_ClassSymbol) == R_NilValue) {
        SHALLOW_DUPLICATE_ATTRIB(ni, x[j]);
      } else {
        Rf_setAttrib(ni, wrap("label"), Rf_getAttrib(x[j], wrap("label")));
      }
      nobs[j] = ni;
    }
    DUPLICATE_ATTRIB(nobs, x);
    nobs.attr("row.names") = NumericVector::create(NA_REAL, -ng);
    return nobs;
  }
}
