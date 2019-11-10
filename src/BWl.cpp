#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List BWlCpp(const List& x, int ng = 0, const IntegerVector& g = 0, 
            const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, 
            bool narm = true, bool option = false, bool B = false) {
  
  int l = x.size();
  List out(l);
  
  if (Rf_isNull(w)) { // No weights !!
    if (ng == 0) {
      if(!B && option) stop("For this return option a grouping vector needs to be supplied");
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          int row = column.size();
          int k = row-1, ni = 1;
          long double sumi = column[k];
          while(std::isnan(sumi) && k!=0) sumi = column[--k];
          if(k != 0) for(int i = k; i--; ) {
            if(std::isnan(column[i])) continue;
            sumi += column[i];
            ++ni;
          }
          sumi = sumi/ni;
          if(B) {
            if(option) out[j] = rep((double)sumi, row); // good ??
            else {
              NumericVector outj = no_init_vector(row);
              for(int i = 0; i != row; ++i) {
                if(std::isnan(column[i])) outj[i] = column[i];
                else outj[i] = sumi; 
              }
              out[j] = outj;
            }
          } else {
            out[j] = column - sumi; 
          }
          SHALLOW_DUPLICATE_ATTRIB(out[j], column); // good ??
        }
      } else {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          long double sumi = 0;
          int row = column.size();
          for(int i = 0; i != row; ++i) {
            if(std::isnan(column[i])) {
              sumi = column[i];
              break;
            } else {
              sumi += column[i];
            }
          }
          sumi = sumi/row;
          if(B) {
            out[j] = rep((double)sumi, row); 
          } else {
            out[j] = column - sumi; 
          }
          SHALLOW_DUPLICATE_ATTRIB(out[j], column);
        }
      }
    } else { // With groups !!
      int gss = g.size();
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          int row = column.size();
          if(gss != row) stop("length(g) must match nrow(X)");
          std::vector<double> sumj(ng, NA_REAL);
          std::vector<int>  nj(ng, 1); 
          for(int i = row; i--; ) {
            if(!std::isnan(column[i])) { 
              if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
              else { 
                sumj[g[i]-1] += column[i];
                ++nj[g[i]-1];
              }
            }
          }
          NumericVector outj = no_init_vector(row);
          if(B) {
            for(int i = ng; i--; ) sumj[i] /= nj[i];
            if(option) {
              for(int i = 0; i != row; ++i) outj[i] = sumj[g[i]-1];
            } else {
              for(int i = 0; i != row; ++i) {
                if(std::isnan(column[i])) outj[i] = column[i];
                else outj[i] = sumj[g[i]-1]; 
              }
            }
          } else {
            if(!option) {
              for(int i = ng; i--; ) sumj[i] /= nj[i];
              for(int i = 0; i != row; ++i) outj[i] = column[i] - sumj[g[i]-1]; 
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
              for(int i = 0; i != row; ++i) outj[i] = column[i] - sumj[g[i]-1] + osum;
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(outj, column);
          out[j] = outj;
        }
      } else {
        if(Rf_isNull(gs)) {
          int gsv[ng], memsize = sizeof(int)*ng;
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            int row = column.size();
            if(gss != row) stop("length(g) must match nrow(X)");
            std::vector<double> sumj(ng); 
            memset(gsv, 0, memsize);
            int ngs = 0;
            for(int i = 0; i != row; ++i) {
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
            NumericVector outj = no_init_vector(row);
            if(B) {
              for(int i = ng; i--; ) sumj[i] /= gsv[i];
              for(int i = 0; i != row; ++i) outj[i] = sumj[g[i]-1];
            } else {
              if(!option) {
                for(int i = ng; i--; ) sumj[i] /= gsv[i];
                for(int i = 0; i != row; ++i) outj[i] = column[i] - sumj[g[i]-1]; 
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
                for(int i = 0; i != row; ++i) outj[i] = column[i] - sumj[g[i]-1] + osum;
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        } else {
          IntegerVector gsv = gs;
          if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups"); 
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            int row = column.size();
            if(gss != row) stop("length(g) must match nrow(X)");
            NumericVector sumj(ng); //  = no_init_vector //  Not initializing seems to be numerically unstable !!!!
            int ngs = 0;
            for(int i = 0; i != row; ++i) {
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
            NumericVector outj = no_init_vector(row);
            if(B) {
              for(int i = ng; i--; ) sumj[i] /= gsv[i];
              for(int i = 0; i != row; ++i) outj[i] = sumj[g[i]-1];
            } else {
              if(!option) {
                for(int i = ng; i--; ) sumj[i] /= gsv[i];
                for(int i = 0; i != row; ++i) outj[i] = column[i] - sumj[g[i]-1]; 
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
                for(int i = 0; i != row; ++i) outj[i] = column[i] - sumj[g[i]-1] + osum;
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        }
      }
    }
  } else { // With weights 
    NumericVector wg = w; 
    int wgs = wg.size();
    if (ng == 0) {
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          int row = column.size();
          if(row != wgs) stop("length(w) must match nrow(X)"); 
          int k = row-1;
          while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k; 
          long double sumi = column[k]*wg[k], sumwi = wg[k];
          if(k != 0) for(int i = k; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
            sumi += column[i]*wg[i];
            sumwi += wg[i];
          }
          sumi = sumi/sumwi;
          if(B) {
            if(option) out[j] = rep((double)sumi, row); 
            else {
              NumericVector outj = no_init_vector(row);
              for(int i = 0; i != row; ++i) {
                if(std::isnan(column[i])) outj[i] = column[i];
                else outj[i] = sumi; 
              }
              out[j] = outj;
            }
          } else {
            out[j] = column - sumi; 
          }
          SHALLOW_DUPLICATE_ATTRIB(out[j], column); // good like this ?? 
        }
      } else {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          int row = column.size();
          if(row != wgs) stop("length(w) must match nrow(X)"); 
          long double sumi = 0, sumwi = 0;
          for(int i = 0; i != row; ++i) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) {
              sumi = column[i]+wg[i];
              break;
            } else {
              sumi += column[i]*wg[i];
              sumwi += wg[i];
            }
          }
          sumi = sumi/sumwi;
          if(B) {
            out[j] = rep((double)sumi, row); 
          } else {
            out[j] = column - sumi; 
          }
          SHALLOW_DUPLICATE_ATTRIB(out[j], column);
        }
      }
    } else { // With groups !!
      int gss = g.size();
      if(wgs != gss) stop("length(w) must match length(g)");
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          int row = column.size();
          if(gss != row) stop("length(g) must match nrow(X)");
          std::vector<double> sumj(ng, NA_REAL);
          double sumwj[ng];
          for(int i = row; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) continue; 
            if(std::isnan(sumj[g[i]-1])) {
              sumj[g[i]-1] = column[i]*wg[i];
              sumwj[g[i]-1] = wg[i];
            } else {
              sumj[g[i]-1] += column[i]*wg[i]; 
              sumwj[g[i]-1] += wg[i];
            }
          }
          NumericVector outj = no_init_vector(row);
          if(B) {
            for(int i = ng; i--; ) sumj[i] /= sumwj[i];
            if(option) {
              for(int i = 0; i != row; ++i) outj[i] = sumj[g[i]-1];
            } else {
              for(int i = 0; i != row; ++i) {
                if(std::isnan(column[i])) outj[i] = column[i];
                else outj[i] = sumj[g[i]-1]; 
              }
            }
          } else {
            if(!option) {
              for(int i = ng; i--; ) sumj[i] /= sumwj[i];
              for(int i = 0; i != row; ++i) outj[i] = column[i] - sumj[g[i]-1]; 
            } else {
              double osum = 0, osumw = 0;
              for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
                if(std::isnan(sumj[i])) continue; // solves the issue !!
                osum += sumj[i];
                osumw += sumwj[i];
                sumj[i] /= sumwj[i]; 
              }
              osum = osum/osumw;
              for(int i = 0; i != row; ++i) outj[i] = column[i] - sumj[g[i]-1] + osum;
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(outj, column);
          out[j] = outj;
        }
      } else {
        double sumj[ng], sumwj[ng];
        int memsize = sizeof(double)*ng;
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          int row = column.size();
          if(gss != row) stop("length(g) must match nrow(X)");
          memset(sumj, 0, memsize);
          memset(sumwj, 0, memsize);
          int ngs = 0;
          for(int i = 0; i != row; ++i) {
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
          NumericVector outj = no_init_vector(row);
          if(B) {
            for(int i = ng; i--; ) sumj[i] /= sumwj[i];
            for(int i = 0; i != row; ++i) outj[i] = sumj[g[i]-1];
          } else {
            if(!option) {
              for(int i = ng; i--; ) sumj[i] /= sumwj[i];
              for(int i = 0; i != row; ++i) outj[i] = column[i] - sumj[g[i]-1]; 
            } else {
              double osum = 0, osumw = 0;
              for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
                if(std::isnan(sumj[i])) continue; // solves the issue !!
                osum += sumj[i];
                osumw += sumwj[i];
                sumj[i] /= sumwj[i]; 
              }
              osum = osum/osumw;
              for(int i = 0; i != row; ++i) outj[i] = column[i] - sumj[g[i]-1] + osum;
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(outj, column);
          out[j] = outj;
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  return out;
}
