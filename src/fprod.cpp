#include <Rcpp.h>
using namespace Rcpp;

// Note: For weighted computations the model code is fmean.cpp.

// [[Rcpp::export]]
NumericVector fprodCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0,
                       const SEXP& w = R_NilValue, bool narm = true) {
  int l = x.size();
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

 if(Rf_isNull(w)) { // No weights
  if(ng == 0) {
    if(narm) {
      int j = l-1;
      long double prod = x[j];
      while(std::isnan(prod) && j!=0) prod = x[--j];
      if(j != 0) for(int i = j; i--; ) {
        if(!std::isnan(x[i])) prod *= x[i]; // Fastest ?
      }
      return Rf_ScalarReal((double)prod);
    } else {
      long double prod = 1;
      for(int i = 0; i != l; ++i) {
        if(std::isnan(x[i])) {
          prod = x[i];
          break;
        } else {
          prod *= x[i];
        }
      }
      return Rf_ScalarReal((double)prod);
    }
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    if(narm) {
      NumericVector prod(ng, NA_REAL); // Other way ?
      for(int i = l; i--; ) {
        if(!std::isnan(x[i])) { // faster way to code this ? -> Not Bad at all
          if(std::isnan(prod[g[i]-1])) prod[g[i]-1] = x[i];
          else prod[g[i]-1] *= x[i];
        }
      }
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, prod);
      return prod;
    } else {
      NumericVector prod(ng, 1.0); // good? -> yes
      int ngs = 0;
      for(int i = 0; i != l; ++i) {
        if(std::isnan(x[i])) {
          if(!std::isnan(prod[g[i]-1])) {
            prod[g[i]-1] = x[i];
            ++ngs;
            if(ngs == ng) break;
          }
        } else {
          prod[g[i]-1] *= x[i];
        }
      }
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, prod);
      return prod;
    }
  }
 } else { // With weights
   NumericVector wg = w;
   if(l != wg.size()) stop("length(w) must match length(x)");
   if(ng == 0) {
     if(narm) {
       int j = l-1;
       while((std::isnan(x[j]) || std::isnan(wg[j])) && j!=0) --j;
       long double prod = x[j]*wg[j];
       if(j != 0) for(int i = j; i--; ) {
         if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
         prod *= x[i]*wg[i];
       }
       return Rf_ScalarReal((double)prod);
     } else {
       long double prod = 1;
       for(int i = 0; i != l; ++i) {
         if(std::isnan(x[i]) || std::isnan(wg[i])) {
           prod = x[i]+wg[i];
           break;
         } else {
           prod *= x[i]*wg[i];
         }
       }
       return Rf_ScalarReal((double)prod);
     }
   } else { // with groups
     if(g.size() != l) stop("length(g) must match nrow(X)");
     if(narm) {
       NumericVector prod(ng, NA_REAL);
       for(int i = l; i--; ) {
         if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
         if(std::isnan(prod[g[i]-1])) prod[g[i]-1] = x[i]*wg[i];
         else prod[g[i]-1] *= x[i]*wg[i];
       }
       if(!Rf_isObject(x)) Rf_copyMostAttrib(x, prod);
       return prod;
     } else {
       NumericVector prod(ng, 1.0);
       int ngs = 0;
       for(int i = 0; i != l; ++i) {
         if(std::isnan(x[i]) || std::isnan(wg[i])) {
           if(!std::isnan(prod[g[i]-1])) {
             prod[g[i]-1] = x[i]+wg[i];
             ++ngs;
             if(ngs == ng) break;
           }
         } else {
           prod[g[i]-1] *= x[i]*wg[i];
         }
       }
       if(!Rf_isObject(x)) Rf_copyMostAttrib(x, prod);
       return prod;
     }
   }
 }
}




// [[Rcpp::export]]
SEXP fprodmCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0,
               const SEXP& w = R_NilValue, bool narm = true, bool drop = true) {
  int l = x.nrow(), col = x.ncol();

 if(Rf_isNull(w)) { // No weights
  if(ng == 0) {
    NumericVector prod = no_init_vector(col); // Initialize faster -> Nope
    if(narm) {
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        int k = l-1;
        long double prodj = column[k];
        while(std::isnan(prodj) && k!=0) prodj = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(!std::isnan(column[i])) prodj *= column[i];
        }
        prod[j] = (double)prodj;
      }
    } else {
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        long double prodj = 1;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(column[i])) {
            prodj = column[i];
            break;
          } else {
            prodj *= column[i];
          }
        }
        prod[j] = (double)prodj;
      }
    }
    if(drop) Rf_setAttrib(prod, R_NamesSymbol, colnames(x));
    else {
      Rf_dimgets(prod, Dimension(1, col));
      colnames(prod) = colnames(x);
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, prod);
    }
    return prod;
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    if(narm) {
      NumericMatrix prod = no_init_matrix(ng, col);
      std::fill(prod.begin(), prod.end(), NA_REAL);
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column prodj = prod( _ , j);
        for(int i = l; i--; ) {
          if(!std::isnan(column[i])) {
            if(std::isnan(prodj[g[i]-1])) prodj[g[i]-1] = column[i];
            else prodj[g[i]-1] *= column[i];
          }
        }
      }
      colnames(prod) = colnames(x);
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, prod);
      return prod;
    } else {
      NumericMatrix prod = no_init_matrix(ng, col); // no init numerically unstable
      std::fill(prod.begin(), prod.end(), 1.0);
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column prodj = prod( _ , j);
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(column[i])) {
            if(!std::isnan(prodj[g[i]-1])) {
              prodj[g[i]-1] = column[i];
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            prodj[g[i]-1] *= column[i];
          }
        }
      }
      colnames(prod) = colnames(x);
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, prod);
      return prod;
    }
  }
 } else { // With weights
   NumericVector wg = w;
   if(l != wg.size()) stop("length(w) must match nrow(X)");
   if(ng == 0) {
     NumericVector prod = no_init_vector(col);
     if(narm) {
       for(int j = col; j--; ) {
         NumericMatrix::ConstColumn column = x( _ , j);
         int k = l-1;
         while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
         long double prodj = column[k]*wg[k];
         if(k != 0) for(int i = k; i--; ) {
           if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
           prodj *= column[i]*wg[i];
         }
         prod[j] = (double)prodj;
       }
     } else {
       for(int j = col; j--; ) {
         NumericMatrix::ConstColumn column = x( _ , j);
         long double prodj = 1;
         for(int i = 0; i != l; ++i) {
           if(std::isnan(column[i]) || std::isnan(wg[i])) {
             prodj = column[i]+wg[i];
             break;
           } else {
             prodj *= column[i]*wg[i];
           }
         }
         prod[j] = (double)prodj;
       }
     }
     if(drop) Rf_setAttrib(prod, R_NamesSymbol, colnames(x));
     else {
       Rf_dimgets(prod, Dimension(1, col));
       colnames(prod) = colnames(x);
       if(!Rf_isObject(x)) Rf_copyMostAttrib(x, prod);
     }
     return prod;
   } else { // with groups
     if(g.size() != l) stop("length(g) must match nrow(X)");
     if(narm) {
       NumericMatrix prod = no_init_matrix(ng, col);
       std::fill(prod.begin(), prod.end(), NA_REAL);
       for(int j = col; j--; ) {
         NumericMatrix::ConstColumn column = x( _ , j);
         NumericMatrix::Column prodj = prod( _ , j);
         for(int i = l; i--; ) {
           if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
           if(std::isnan(prodj[g[i]-1])) {
             prodj[g[i]-1] = column[i]*wg[i];
           } else {
             prodj[g[i]-1] *= column[i]*wg[i];
           }
         }
       }
       colnames(prod) = colnames(x);
       if(!Rf_isObject(x)) Rf_copyMostAttrib(x, prod);
       return prod;
     } else {
       NumericMatrix prod = no_init_matrix(ng, col);
       std::fill(prod.begin(), prod.end(), 1.0);
       for(int j = col; j--; ) {
         NumericMatrix::ConstColumn column = x( _ , j);
         NumericMatrix::Column prodj = prod( _ , j);
         int ngs = 0;
         for(int i = 0; i != l; ++i) {
           if(std::isnan(column[i]) || std::isnan(wg[i])) {
             if(!std::isnan(prodj[g[i]-1])) {
               prodj[g[i]-1] = column[i]+wg[i];
               ++ngs;
               if(ngs == ng) break;
             }
           } else {
             prodj[g[i]-1] *= column[i]*wg[i];
           }
         }
       }
       colnames(prod) = colnames(x);
       if(!Rf_isObject(x)) Rf_copyMostAttrib(x, prod);
       return prod;
     }
   }
 }
}




// [[Rcpp::export]]
SEXP fprodlCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
               const SEXP& w = R_NilValue, bool narm = true, bool drop = true) {
  int l = x.size();

 if(Rf_isNull(w)) { // No weights
  if(ng == 0) {
    NumericVector prod(l);
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        int k = column.size()-1;
        long double prodi = column[k];
        while(std::isnan(prodi) && k!=0) prodi = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(!std::isnan(column[i])) prodi *= column[i];
        }
        prod[j] = (double)prodi;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        long double prodi = 1;
        int row = column.size();
        for(int i = 0; i != row; ++i) {
          if(std::isnan(column[i])) {
            prodi = column[i];
            break;
          } else {
            prodi *= column[i];
          }
        }
        prod[j] = (double)prodi;
      }
    }
    if(drop) {
      Rf_setAttrib(prod, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
      return prod;
    } else {
      List out(l);
      for(int j = l; j--; ) {
        out[j] = prod[j];
        SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      Rf_setAttrib(out, R_RowNamesSymbol, Rf_ScalarInteger(1));
      return out;
    }
  } else { // With groups
    List prod(l);
    int gss = g.size();
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector prodj(ng, NA_REAL);
        for(int i = gss; i--; ) {
          if(!std::isnan(column[i])) {
            if(std::isnan(prodj[g[i]-1])) prodj[g[i]-1] = column[i];
            else prodj[g[i]-1] *= column[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(prodj, column);
        prod[j] = prodj;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector prodj(ng, 1.0);
        int ngs = 0;
        for(int i = 0; i != gss; ++i) {
          if(std::isnan(column[i])) {
            if(!std::isnan(prodj[g[i]-1])) {
              prodj[g[i]-1] = column[i];
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            prodj[g[i]-1] *= column[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(prodj, column);
        prod[j] = prodj;
      }
    }
    SHALLOW_DUPLICATE_ATTRIB(prod, x);
    Rf_setAttrib(prod, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -ng));
    return prod;
  }
 } else { // With weights
   NumericVector wg = w;
   int wgs = wg.size();
   if (ng == 0) {
     NumericVector prod(l);
     if(narm) {
       for(int j = l; j--; ) {
         NumericVector column = x[j];
         if(column.size() != wgs) stop("length(w) must match nrow(X)");
         int k = wgs-1;
         while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
         long double prodi = column[k]*wg[k];
         if(k != 0) for(int i = k; i--; ) {
           if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
           prodi *= column[i]*wg[i];
         }
         prod[j] = (double)prodi;
       }
     } else {
       for(int j = l; j--; ) {
         NumericVector column = x[j];
         if(column.size() != wgs) stop("length(w) must match nrow(X)");
         long double prodi = 1;
         for(int i = 0; i != wgs; ++i) {
           if(std::isnan(column[i]) || std::isnan(wg[i])) {
             prodi = column[i]+wg[i];
             break;
           } else {
             prodi *= column[i]*wg[i];
           }
         }
         prod[j] = (double)prodi;
       }
     }
     if(drop) {
       Rf_setAttrib(prod, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
       return prod;
     } else {
       List out(l);
       for(int j = l; j--; ) {
         out[j] = prod[j];
         SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
       }
       SHALLOW_DUPLICATE_ATTRIB(out, x);
       Rf_setAttrib(out, R_RowNamesSymbol, Rf_ScalarInteger(1));
       return out;
     }
   } else { // With groups
     List prod(l);
     int gss = g.size();
     if(wgs != gss) stop("length(w) must match length(g)");
     if(narm) {
       for(int j = l; j--; ) {
         NumericVector column = x[j];
         if(gss != column.size()) stop("length(g) must match nrow(X)");
         NumericVector prodj(ng, NA_REAL);
         for(int i = gss; i--; ) {
           if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
           if(std::isnan(prodj[g[i]-1])) {
             prodj[g[i]-1] = column[i]*wg[i];
           } else {
             prodj[g[i]-1] *= column[i]*wg[i];
           }
         }
         SHALLOW_DUPLICATE_ATTRIB(prodj, column);
         prod[j] = prodj;
       }
     } else {
       for(int j = l; j--; ) {
         NumericVector column = x[j];
         if(gss != column.size()) stop("length(g) must match nrow(X)");
         NumericVector prodj(ng, 1.0);
         int ngs = 0;
         for(int i = 0; i != gss; ++i) {
           if(std::isnan(column[i]) || std::isnan(wg[i])) {
             if(!std::isnan(prodj[g[i]-1])) {
               prodj[g[i]-1] = column[i]+wg[i];
               ++ngs;
               if(ngs == ng) break;
             }
           } else {
             prodj[g[i]-1] *= column[i]*wg[i];
           }
         }
         SHALLOW_DUPLICATE_ATTRIB(prodj, column);
         prod[j] = prodj;
       }
     }
     SHALLOW_DUPLICATE_ATTRIB(prod, x);
     Rf_setAttrib(prod, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -ng));
     return prod;
   }
 }
}
