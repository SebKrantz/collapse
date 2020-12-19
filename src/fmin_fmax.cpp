#include <Rcpp.h>
using namespace Rcpp;

template <typename F> // Faster non constant references here ! (for some weird reason)
NumericVector fminmaxCppImpl(NumericVector x, int ng = 0, IntegerVector g = 0,
                             bool narm = true, F FUN = [](double a, double b) { return a > b; }, double init = R_PosInf) {
  int l = x.size();
  if(l < 2) return x; // Prevents seqfault for numeric(0) #101

  if(ng == 0) {
    if(narm) {
      int j = l-1;
      double min = x[j];
      while(std::isnan(min) && j!=0) min = x[--j];
      if(j != 0) for(int i = j; i--; ) {
        if(FUN(min, x[i])) min = x[i];
      }
      return Rf_ScalarReal(min);
    } else {
      double min = x[0];
      for(int i = 0; i != l; ++i) {
        if(std::isnan(x[i])) {
          min = x[i];
          break;
        } else {
          if(FUN(min, x[i])) min = x[i];
        }
      }
      return Rf_ScalarReal(min);
    }
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    if(narm) {
      NumericVector min(ng, NA_REAL); // Other way ?
      for(int i = l; i--; ) { // adding if isnan(x[i]) before is not faster
        if(FUN(min[g[i]-1], x[i]) || std::isnan(min[g[i]-1])) min[g[i]-1] = x[i];  // fastest
      }
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, min);
      return min;
    } else {
      NumericVector min(ng, init); // INFINITY // DBL_MAX // good? -> yes, same value but better output !
      int ngs = 0;
      for(int i = 0; i != l; ++i) {
        if(std::isnan(x[i])) {
          if(!std::isnan(min[g[i]-1])) {
            min[g[i]-1] = x[i];
            ++ngs;
            if(ngs == ng) break;
          }
        } else {
          if(FUN(min[g[i]-1], x[i])) min[g[i]-1] = x[i];
        }
      }
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, min);
      return min;
    }
  }
}


// [[Rcpp::export]]
NumericVector fminmaxCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0,
                         bool narm = true, int ret = 1) {
  if(ret == 1) return fminmaxCppImpl(x, ng, g, narm, [](double a, double b) { return a > b; }, R_PosInf);
  return fminmaxCppImpl(x, ng, g, narm, [](double a, double b) { return a < b; }, R_NegInf);
}


template <typename F>
SEXP fminmaxmCppImpl(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0,
                     bool narm = true, bool drop = true, F FUN = [](double a, double b) { return a > b; }, double init = R_PosInf) {
  int l = x.nrow(), col = x.ncol();

  if(ng == 0) {
    NumericVector min = no_init_vector(col); // Initialize faster -> Nope
    if(narm) {
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        int k = l-1;
        double minj = column[k];
        while(std::isnan(minj) && k!=0) minj = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(FUN(minj, column[i])) minj = column[i]; // continue here
        }
        min[j] = minj;
      }
    } else {
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        double minj = column[0];
        for(int i = 0; i != l; ++i) {
          if(std::isnan(column[i])) {
            minj = column[i];
            break;
          } else {
            if(FUN(minj, column[i])) minj = column[i];
          }
        }
        min[j] = minj;
      }
    }
    if(drop) Rf_setAttrib(min, R_NamesSymbol, colnames(x));
    else {
      Rf_dimgets(min, Dimension(1, col));
      colnames(min) = colnames(x);
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, min);
    }
    return min;
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    NumericMatrix min = no_init_matrix(ng, col);
    if(narm) {
      std::fill(min.begin(), min.end(), NA_REAL);
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column minj = min( _ , j);
        for(int i = l; i--; ) {
          if(!std::isnan(column[i])) { // Keeping this is faster !
            if(FUN(minj[g[i]-1], column[i]) || std::isnan(minj[g[i]-1])) minj[g[i]-1] = column[i];
          }
        }
      }
    } else {
      std::fill(min.begin(), min.end(), init);
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column minj = min( _ , j);
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(column[i])) {
            if(!std::isnan(minj[g[i]-1])) {
              minj[g[i]-1] = column[i];
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            if(FUN(minj[g[i]-1], column[i])) minj[g[i]-1] = column[i];
          }
        }
      }
    }
    colnames(min) = colnames(x);
    if(!Rf_isObject(x)) Rf_copyMostAttrib(x, min);
    return min;
  }
}

// [[Rcpp::export]]
SEXP fminmaxmCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0,
                 bool narm = true, bool drop = true, int ret = 1) {
  if(ret == 1) return fminmaxmCppImpl(x, ng, g, narm, drop, [](double a, double b) { return a > b; }, R_PosInf);
  return fminmaxmCppImpl(x, ng, g, narm, drop, [](double a, double b) { return a < b; }, R_NegInf);
}


template <typename F>
SEXP fminmaxlCppImpl(const List& x, int ng = 0, const IntegerVector& g = 0,
                     bool narm = true, bool drop = true, F FUN = [](double a, double b) { return a > b; }, double init = R_PosInf) {
  int l = x.size();

  if (ng == 0) {
    NumericVector min = no_init_vector(l); // good and fast here !
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        int k = column.size()-1;
        double mini = column[k];
        while(std::isnan(mini) && k!=0) mini = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(FUN(mini, column[i])) mini = column[i];
        }
        min[j] = mini;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        double mini = column[0];
        int row = column.size();
        for(int i = 0; i != row; ++i) {
          if(std::isnan(column[i])) {
            mini = column[i];
            break;
          } else {
            if(FUN(mini, column[i])) mini = column[i];
          }
        }
        min[j] = mini;
      }
    }
    if(drop) {
      Rf_setAttrib(min, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
      return min;
    } else {
      List out(l);
      for(int j = l; j--; ) {
        out[j] = min[j];
        SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
      }
      DUPLICATE_ATTRIB(out, x);
      Rf_setAttrib(out, R_RowNamesSymbol, Rf_ScalarInteger(1));
      return out;
    }
  } else { // With groups
    List min(l);
    int gss = g.size();
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector minj(ng, NA_REAL);
        for(int i = gss; i--; ) {
          if(!std::isnan(column[i])) { // Keeping this is faster !
            if(FUN(minj[g[i]-1], column[i]) || std::isnan(minj[g[i]-1])) minj[g[i]-1] = column[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(minj, column);
        min[j] = minj;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector minj(ng, init);
        int ngs = 0;
        for(int i = 0; i != gss; ++i) {
          if(std::isnan(column[i])) {
            if(!std::isnan(minj[g[i]-1])) {
              minj[g[i]-1] = column[i];
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            if(FUN(minj[g[i]-1], column[i])) minj[g[i]-1] = column[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(minj, column);
        min[j] = minj;
      }
    }
    DUPLICATE_ATTRIB(min, x);
    Rf_setAttrib(min, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -ng));
    return min;
  }
}

// [[Rcpp::export]]
SEXP fminmaxlCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
                 bool narm = true, bool drop = true, int ret = 1) {
  if(ret == 1) return fminmaxlCppImpl(x, ng, g, narm, drop, [](double a, double b) { return a > b; }, R_PosInf);
  return fminmaxlCppImpl(x, ng, g, narm, drop, [](double a, double b) { return a < b; }, R_NegInf);
}
