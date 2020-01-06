#include <Rcpp.h>
using namespace Rcpp;
// // [[Rcpp::interfaces(cpp)]] // not useful!!

// [[Rcpp::export]] // // [[Rcpp::export(.fmaxCpp)]]
NumericVector fmaxCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0,
                      bool narm = true) {
  int l = x.size();

  if(ng == 0) {
    if(narm) {
      int j = l-1;
      double max = x[j];
      while(std::isnan(max) && j!=0) max = x[--j];
      if(j != 0) for(int i = j; i--; ) {
        if(max < x[i]) max = x[i];
      }
      return NumericVector::create(max);
    } else {
      double max = x[0];
      for(int i = 0; i != l; ++i) {
        if(std::isnan(x[i])) {
          max = x[i];
          break;
        } else {
          if(max < x[i]) max = x[i];
        }
      }
      return NumericVector::create(max);
    }
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    if(narm) {
      NumericVector max(ng, NA_REAL); // Other way ??
      for(int i = l; i--; ) { // adding if isnan(x[i]) before is not faster !!!
        if(max[g[i]-1] < x[i] || std::isnan(max[g[i]-1])) max[g[i]-1] = x[i];  // fastest !!
      }
      DUPLICATE_ATTRIB(max, x);
      return max;
    } else {
      NumericVector max(ng, R_NegInf); // -INFINITY // good?? -> yes
      int ngs = 0;
      for(int i = 0; i != l; ++i) {
        if(std::isnan(x[i])) {
          if(!std::isnan(max[g[i]-1])) {
            max[g[i]-1] = x[i];
            ++ngs;
            if(ngs == ng) break;
          }
        } else {
          if(max[g[i]-1] < x[i]) max[g[i]-1] = x[i];
        }
      }
      DUPLICATE_ATTRIB(max, x);
      return max;
    }
  }
}




// [[Rcpp::export]]
SEXP fmaxmCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0,
              bool narm = true, bool drop = true) {
  int l = x.nrow(), col = x.ncol();

  if(ng == 0) {
    NumericVector max = no_init_vector(col); // Initialize faster -> Nope !!!
    if(narm) {
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        int k = l-1;
        double maxj = column[k];
        while(std::isnan(maxj) && k!=0) maxj = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(maxj < column[i]) maxj = column[i];
        }
        max[j] = maxj;
      }
    } else {
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        double maxj = column[0];
        for(int i = 0; i != l; ++i) {
          if(std::isnan(column[i])) {
            maxj = column[i];
            break;
          } else {
            if(maxj < column[i]) maxj = column[i];
          }
        }
        max[j] = maxj;
      }
    }
    if(drop) max.attr("names") = colnames(x);
    else {
      max.attr("dim") = Dimension(1, col);
      colnames(max) = colnames(x);
    }
    return max;
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    NumericMatrix max = no_init_matrix(ng, col);
    if(narm) {
      std::fill(max.begin(), max.end(), NA_REAL);
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column maxj = max( _ , j);
        for(int i = l; i--; ) {
          if(!std::isnan(column[i])) { // Keeping this is faster !!!!
            if(maxj[g[i]-1] < column[i] || std::isnan(maxj[g[i]-1])) maxj[g[i]-1] = column[i];
          }
        }
      }
    } else {
      std::fill(max.begin(), max.end(), R_NegInf); // -INFINITY
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column maxj = max( _ , j);
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(column[i])) {
            if(!std::isnan(maxj[g[i]-1])) {
              maxj[g[i]-1] = column[i];
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            if(maxj[g[i]-1] < column[i]) maxj[g[i]-1] = column[i];
          }
        }
      }
    }
    colnames(max) = colnames(x);
    return max;
  }
}




// [[Rcpp::export]]
SEXP fmaxlCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
              bool narm = true, bool drop = true) {
  int l = x.size();

  if (ng == 0) {
    NumericVector max = no_init_vector(l); // Good and faster !!
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        int k = column.size()-1;
        double maxi = column[k];
        while(std::isnan(maxi) && k!=0) maxi = column[--k];
        if(k != 0) for(int i = k; i--; ) { // This is fastest, checking isnan(column[i]) does not give extra speed!!
          if(maxi < column[i]) maxi = column[i];
        }
        max[j] = maxi;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        double maxi = column[0];
        int row = column.size();
        for(int i = 0; i != row; ++i) {
          if(std::isnan(column[i])) {
            maxi = column[i];
            break;
          } else {
            if(maxi < column[i]) maxi = column[i];
          }
        }
        max[j] = maxi;
      }
    }
    if(drop) {
      max.attr("names") = x.attr("names");
      return max;
    } else {
      List out(l);
      for(int j = l; j--; ) {
        out[j] = max[j];
        SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
      }
      DUPLICATE_ATTRIB(out, x);
      out.attr("row.names") = 1;
      return out;
    }
  } else { // With groups !!
    List max(l);
    int gss = g.size();
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector maxj(ng, NA_REAL);
        for(int i = gss; i--; ) {
          if(!std::isnan(column[i])) { // Keeping this is faster !!!!
            if(maxj[g[i]-1] < column[i] || std::isnan(maxj[g[i]-1])) maxj[g[i]-1] = column[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(maxj, column);
        max[j] = maxj;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector maxj(ng, R_NegInf); // -INFINITY
        int ngs = 0;
        for(int i = 0; i != gss; ++i) {
          if(std::isnan(column[i])) {
            if(!std::isnan(maxj[g[i]-1])) {
              maxj[g[i]-1] = column[i];
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            if(maxj[g[i]-1] < column[i]) maxj[g[i]-1] = column[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(maxj, column);
        max[j] = maxj;
      }
    }
    DUPLICATE_ATTRIB(max, x);
    max.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng); // NumericVector::create(NA_REAL, -ng);
    return max;
  }
}
