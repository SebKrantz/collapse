#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix BWmCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0, 
                     const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, 
                     bool narm = true, bool option = false, bool B = false) { 
  int l = x.nrow(), col = x.ncol(); 
  NumericMatrix out = no_init_matrix(l, col);
  
  if (Rf_isNull(w)) { // No weights !!
    if(ng == 0) { 
      if(!B && option) stop("For this return option a grouping vector needs to be supplied");
      if(narm) { 
        for(int j = col; j--; ) { // Instead Am(j,_) you can use Am.row(j).
          NumericMatrix::ConstColumn column = x( _ , j); 
          NumericMatrix::Column outj = out( _ , j); 
          int k = l-1, nj = 1;
          long double sumj = column[k]; 
          while(std::isnan(sumj) && k!=0) sumj = column[--k];
          if(k != 0) for(int i = k; i--; ) {
            if(std::isnan(column[i])) continue;
            sumj += column[i]; 
            ++nj;
          }
          sumj = sumj/nj;
          if(B) {
            if(option) std::fill(outj.begin(), outj.end(), (double)sumj); 
            else {
              for(int i = 0; i != l; ++i) {
                if(std::isnan(column[i])) outj[i] = column[i];
                else outj[i] = sumj; 
              }
            }
          } else {
            outj = column - sumj; 
          }
        }
      } else {
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j);
          NumericMatrix::Column outj = out( _ , j); 
          long double sumj = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(column[i])) {
              sumj = column[i]; 
              break;
            } else { 
              sumj += column[i];
            }
          }
          sumj = sumj/l;
          if(B) {
            std::fill(outj.begin(), outj.end(), (double)sumj); 
          } else {
            outj = column - sumj; 
          }
        }
      }
    } else { // with groups 
      if(g.size() != l) stop("length(g) must match nrow(X)");
      if(narm) {
        for(int j = col; j--; ) { 
          NumericMatrix::ConstColumn column = x( _ , j); 
          NumericMatrix::Column outj = out( _ , j);
          std::vector<double> sumj(ng, NA_REAL); // faster than NumericVector ?? 
          int nj[ng]; // use vector also ?? 
          for(int i = l; i--; ) { 
            if(!std::isnan(column[i])) { 
              if(std::isnan(sumj[g[i]-1])) {
                sumj[g[i]-1] = column[i];
                nj[g[i]-1] = 1;
              } else {
                sumj[g[i]-1] += column[i];
                ++nj[g[i]-1];
              }
            }
          }
          if(B) {
            for(int i = ng; i--; ) sumj[i] /= nj[i];
            if(option) {
              for(int i = 0; i != l; ++i) outj[i] = sumj[g[i]-1];
            } else {
              for(int i = 0; i != l; ++i) {
                if(std::isnan(column[i])) outj[i] = column[i];
                else outj[i] = sumj[g[i]-1]; 
              }
            }
          } else {
            if(!option) {
              for(int i = ng; i--; ) sumj[i] /= nj[i];
              for(int i = 0; i != l; ++i) outj[i] = column[i] - sumj[g[i]-1]; 
            } else {
              int on = 0;
              double osum = 0;
              for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
                if(std::isnan(sumj[i])) continue; // solves the issue !!
                osum += sumj[i];
                on += nj[i];
                sumj[i] /= nj[i]; 
              }
              osum = osum/on;
              for(int i = 0; i != l; ++i) outj[i] = column[i] - sumj[g[i]-1] + osum;
            }
          }
        }
      } else {
        if(Rf_isNull(gs)) {
          int gsv[ng], memsize = sizeof(int)*ng;
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j); 
            NumericMatrix::Column outj = out( _ , j); 
            std::vector<double> sumj(ng); // better than array or NumericVector ??
            memset(gsv, 0, memsize); 
            int ngs = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(column[i])) { 
                if(!std::isnan(sumj[g[i]-1])) {
                  sumj[g[i]-1] = column[i]; 
                  ++ngs;
                  if(ngs == ng) break;
                }
              } else { 
                sumj[g[i]-1] += column[i];
                ++gsv[g[i]-1];
              }
            }
            if(B) {
              for(int i = ng; i--; ) sumj[i] /= gsv[i];
              for(int i = 0; i != l; ++i) outj[i] = sumj[g[i]-1];
            } else {
              if(!option) {
                for(int i = ng; i--; ) sumj[i] /= gsv[i];
                for(int i = 0; i != l; ++i) outj[i] = column[i] - sumj[g[i]-1]; 
              } else {
                int on = 0;
                double osum = 0;
                for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
                  if(std::isnan(sumj[i])) continue; // solves the issue !!
                  osum += sumj[i];
                  on += gsv[i];
                  sumj[i] /= gsv[i]; 
                }
                osum = osum/on;
                for(int i = 0; i != l; ++i) outj[i] = column[i] - sumj[g[i]-1] + osum;
              }
            }
          }      
        } else {
          IntegerVector gsv = gs;
          if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups"); 
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j); 
            NumericMatrix::Column outj = out( _ , j);
            std::vector<double> sumj(ng);
            int ngs = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(column[i])) { 
                if(!std::isnan(sumj[g[i]-1])) {
                  sumj[g[i]-1] = column[i]; 
                  ++ngs;
                  if(ngs == ng) break;
                }
              } else { 
                sumj[g[i]-1] += column[i];
              }
            }
            if(B) {
              for(int i = ng; i--; ) sumj[i] /= gsv[i];
              for(int i = 0; i != l; ++i) outj[i] = sumj[g[i]-1];
            } else {
              if(!option) {
                for(int i = ng; i--; ) sumj[i] /= gsv[i];
                for(int i = 0; i != l; ++i) outj[i] = column[i] - sumj[g[i]-1]; 
              } else {
                int on = 0;
                double osum = 0;
                for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
                  if(std::isnan(sumj[i])) continue; // solves the issue !!
                  osum += sumj[i];
                  on += gsv[i];
                  sumj[i] /= gsv[i]; 
                }
                osum = osum/on;
                for(int i = 0; i != l; ++i) outj[i] = column[i] - sumj[g[i]-1] + osum;
              }
            }
          }
        }
      }
    }
  } else { // With weights
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match nrow(X)");
    if(ng == 0) { 
      if(!B && option) stop("For this return option a grouping vector needs to be supplied");
      if(narm) { 
        for(int j = col; j--; ) { // Instead Am(j,_) you can use Am.row(j).
          NumericMatrix::ConstColumn column = x( _ , j);
          NumericMatrix::Column outj = out( _ , j);
          int k = l-1;
          while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k; 
          long double sumj = column[k]*wg[k], sumwj = wg[k];
          if(k != 0) for(int i = k; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
            sumj += column[i]*wg[i];
            sumwj += wg[i];
          }
          sumj = sumj/sumwj;
          if(B) {
            if(option) std::fill(outj.begin(), outj.end(), (double)sumj); 
            else {
              for(int i = 0; i != l; ++i) {
                if(std::isnan(column[i])) outj[i] = column[i];
                else outj[i] = sumj; 
              }
            }
          } else {
            outj = column - sumj; 
          }
        }
      } else {
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j);
          NumericMatrix::Column outj = out( _ , j);
          long double sumj = 0, sumwj = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) {
              sumj = column[i]+wg[i]; 
              break;
            } else { 
              sumj += column[i]*wg[i];
              sumwj += wg[i];
            }
          }
          sumj = sumj/sumwj;
          if(B) {
            std::fill(outj.begin(), outj.end(), (double)sumj); 
          } else {
            outj = column - sumj; 
          }
        }
      }
    } else { // with groups 
      if(g.size() != l) stop("length(g) must match nrow(X)");
      if(narm) {
        for(int j = col; j--; ) { 
          NumericMatrix::ConstColumn column = x( _ , j); 
          NumericMatrix::Column outj = out( _ , j);
          std::vector<double> sumj(ng, NA_REAL); // best ?? 
          double sumwj[ng]; // also make std::vector ?? 
          for(int i = l; i--; ) { 
            if(std::isnan(column[i]) || std::isnan(wg[i])) continue; 
            if(std::isnan(sumj[g[i]-1])) {
              sumj[g[i]-1] = column[i]*wg[i];
              sumwj[g[i]-1] = wg[i];
            } else {
              sumj[g[i]-1] += column[i]*wg[i]; 
              sumwj[g[i]-1] += wg[i];
            }
          }
          if(B) {
            for(int i = ng; i--; ) sumj[i] /= sumwj[i];
            if(option) {
              for(int i = 0; i != l; ++i) outj[i] = sumj[g[i]-1];
            } else {
              for(int i = 0; i != l; ++i) {
                if(std::isnan(column[i])) outj[i] = column[i];
                else outj[i] = sumj[g[i]-1]; 
              }
            }
          } else {
            if(!option) {
              for(int i = ng; i--; ) sumj[i] /= sumwj[i];
              for(int i = 0; i != l; ++i) outj[i] = column[i] - sumj[g[i]-1]; 
            } else {
              double osum = 0, osumw = 0;
              for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
                if(std::isnan(sumj[i])) continue; // solves the issue !!
                osum += sumj[i];
                osumw += sumwj[i];
                sumj[i] /= sumwj[i]; 
              }
              osum = osum/osumw;
              for(int i = 0; i != l; ++i) outj[i] = column[i] - sumj[g[i]-1] + osum;
            }
          }
        }
      } else {
        double sumj[ng], sumwj[ng]; 
        int memsize = sizeof(double)*ng;
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j); 
          NumericMatrix::Column outj = out( _ , j);
          memset(sumj, 0, memsize);
          memset(sumwj, 0, memsize);
          int ngs = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) { 
              if(!std::isnan(sumj[g[i]-1])) {
                sumj[g[i]-1] = sumwj[g[i]-1] = column[i]+wg[i]; // or NA_REAL ?? -> Nope, good !!
                ++ngs;
                if(ngs == ng) break;
              }
            } else {
              sumj[g[i]-1] += column[i]*wg[i];
              sumwj[g[i]-1] += wg[i];
            }
          }
          if(B) {
            for(int i = ng; i--; ) sumj[i] /= sumwj[i];
            for(int i = 0; i != l; ++i) outj[i] = sumj[g[i]-1];
          } else {
            if(!option) {
              for(int i = ng; i--; ) sumj[i] /= sumwj[i];
              for(int i = 0; i != l; ++i) outj[i] = column[i] - sumj[g[i]-1]; 
            } else {
              double osum = 0, osumw = 0;
              for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
                if(std::isnan(sumj[i])) continue; // solves the issue !!
                osum += sumj[i];
                osumw += sumwj[i];
                sumj[i] /= sumwj[i]; 
              }
              osum = osum/osumw;
              for(int i = 0; i != l; ++i) outj[i] = column[i] - sumj[g[i]-1] + osum;
            }
          }
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  return out;
}







