#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fmeanlCpp(const List& x, int ng = 0, const IntegerVector& g = 0, const SEXP& gs = R_NilValue,
               const SEXP& w = R_NilValue, bool narm = true, bool drop = true) {
  int l = x.size();
  
  if(Rf_isNull(w)) { // No weights !!
    if(ng == 0) {
      NumericVector sum(l); // not initializing not faster WIth NWDI (35 instead of 32 milliseconds)
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          int k = column.size()-1, ni = 1;
          long double sumi = column[k]; // long double gives 45 instead of 35 milliseconds !!!
          while(std::isnan(sumi) && k!=0) sumi = column[--k];
          if(k != 0) for(int i = k; i--; ) {
            if(std::isnan(column[i])) continue;
            sumi += column[i];
            ++ni;
          }
          sumi = sumi/ni;
          sum[j] = (double)sumi; 
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
          sum[j] = (double)sumi;
        }
      }
      if(drop) {
        sum.attr("names") = x.attr("names");
        return sum;
      } else {
        List out(l);
        for(int j = l; j--; ) {
          out[j] = sum[j];
          SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
        }
        DUPLICATE_ATTRIB(out, x);
        out.attr("row.names") = 1;
        return out;
      }
    } else { // With groups !!
      List sum(l);
      int gss = g.size();
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          NumericVector sumj(ng, NA_REAL);
          // IntegerVector nj(ng, 1); // faster than no_init_vector ?? -> good, cannot divide by interger 0!!, also else numerically unstable and no speed loss !!
          std::vector<int>  nj(ng, 1); // better memory allocation, and nearly same speed as integer array -> doesn't work because sets all byte to 1 -> https://stackoverflow.com/questions/14761015/memset-an-array-to-1
          for(int i = gss; i--; ) {
            if(!std::isnan(column[i])) { // faster way to code this ??? -> Not Bad at all, 54.. millisec for WDIM
              if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
              else { 
                sumj[g[i]-1] += column[i];
                ++nj[g[i]-1];
              }
            }
          }
          for(int i = ng; i--; ) sumj[i] /= nj[i];
          SHALLOW_DUPLICATE_ATTRIB(sumj, column);
          sum[j] = sumj;
        }
      } else {
        if(Rf_isNull(gs)) {
          int gsv[ng], memsize = sizeof(int)*ng;
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector sumj(ng); //  = no_init_vector //  Not initializing seems to be numerically unstable !!!!
            memset(gsv, 0, memsize);
            int ngs = 0;
            for(int i = 0; i != gss; ++i) {
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
            for(int i = ng; i--; ) sumj[i] /= gsv[i];
            SHALLOW_DUPLICATE_ATTRIB(sumj, column);
            sum[j] = sumj;
          }
        } else {
          IntegerVector gsv = gs;
          if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups"); 
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector sumj(ng); //  = no_init_vector //  Not initializing seems to be numerically unstable !!!!
            int ngs = 0;
            for(int i = 0; i != gss; ++i) {
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
            for(int i = ng; i--; ) sumj[i] /= gsv[i];
            SHALLOW_DUPLICATE_ATTRIB(sumj, column);
            sum[j] = sumj;
          }
        }
      }
      DUPLICATE_ATTRIB(sum, x);
      sum.attr("row.names") = NumericVector::create(NA_REAL, -ng);
      return sum;
    }
  } else { // With weights 
    NumericVector wg = w; // wg(w) //  No speed loss ?? -> Yes, and possibly even faster !!
    int wgs = wg.size();
    if (ng == 0) {
      NumericVector sum(l); // not initializing not faster WIth NWDI (35 instead of 32 milliseconds)
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          if(column.size() != wgs) stop("length(w) must match nrow(X)"); // Really necessary ?? 
          int k = wgs-1;
          while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k; 
          long double sumi = column[k]*wg[k], sumwi = wg[k];
          if(k != 0) for(int i = k; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
            sumi += column[i]*wg[i];
            sumwi += wg[i];
          }
          sumi = sumi/sumwi;
          sum[j] = (double)sumi;
        }
      } else {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          if(column.size() != wgs) stop("length(w) must match nrow(X)"); // Really necessary ?? 
          long double sumi = 0, sumwi = 0;
          for(int i = 0; i != wgs; ++i) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) {
              sumi = column[i]+wg[i];
              break;
            } else {
              sumi += column[i]*wg[i];
              sumwi += wg[i];
            }
          }
          sumi = sumi/sumwi;
          sum[j] = (double)sumi;
        }
      }
      if(drop) {
        sum.attr("names") = x.attr("names");
        return sum;
      } else {
        List out(l);
        for(int j = l; j--; ) {
          out[j] = sum[j];
          SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
        }
        DUPLICATE_ATTRIB(out, x);
        out.attr("row.names") = 1;
        return out;
      }
    } else { // With groups !!
      List sum(l);
      int gss = g.size();
      if(wgs != gss) stop("length(w) must match length(g)");
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          NumericVector sumj(ng, NA_REAL);
          // NumericVector sumwj = no_init_vector(ng); // no_init_vector is faster and stable !!! (you only divide by it every round)
          double sumwj[ng];
          for(int i = gss; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) continue; 
            if(std::isnan(sumj[g[i]-1])) {
              sumj[g[i]-1] = column[i]*wg[i];
              sumwj[g[i]-1] = wg[i];
            } else {
              sumj[g[i]-1] += column[i]*wg[i]; 
              sumwj[g[i]-1] += wg[i];
            }
          }
          // sumj = sumj/sumwj;
          for(int i = ng; i--; ) sumj[i] /= sumwj[i];
          SHALLOW_DUPLICATE_ATTRIB(sumj, column);
          sum[j] = sumj;
        }
      } else {
        double sumwj[ng];
        int memsize = sizeof(double)*ng;
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          NumericVector sumj(ng); //  = no_init_vector //  Not initializing seems to be numerically unstable !!!!
          // NumericVector sumwj(ng); // Also here not initializing is numerically unstable 
          memset(sumwj, 0, memsize);
          int ngs = 0;
          for(int i = 0; i != gss; ++i) {
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
          // sumj = sumj/sumwj;
          for(int i = ng; i--; ) sumj[i] /= sumwj[i];
          SHALLOW_DUPLICATE_ATTRIB(sumj, column);
          sum[j] = sumj;
        }
      }
      DUPLICATE_ATTRIB(sum, x);
      sum.attr("row.names") = NumericVector::create(NA_REAL, -ng);
      return sum;
    }
  }
}


// Final Version without weights !!!
// // [[Rcpp::export]]
// SEXP fgmeanlCpp(List x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0,
//                 bool narm = true, bool drop = true) {
//   int l = x.size(), row = 0;
//   
//   if (ng == 0) {
//     NumericVector sum(l); // not initializing not faster WIth NWDI (35 instead of 32 milliseconds)
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         int k = row-1, ni = 1;
//         double sumi = column[k];
//         while(std::isnan(sumi) && k!=0) sumi = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(std::isnan(column[i])) continue;
//             sumi += column[i];
//             ++ni;
//         }
//         sum[j] = sumi/ni;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double sumi = 0;
//         row = column.size();
//         for(int i = 0; i != row; ++i) {
//           if(std::isnan(column[i])) {
//             sumi = column[i];
//             break;
//           } else {
//             sumi += column[i];
//           }
//         }
//         sum[j] = sumi/row;
//       }
//     }
//     if(drop) {
//       sum.attr("names") = x.attr("names");
//       return sum;
//     } else {
//       List out(l);
//       for(int j = l; j--; ) out[j] = sum[j];
//       return out;
//     }
//   } else { // With groups !!
//     List sum(l);
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng, NA_REAL);
//         IntegerVector nj(ng, 1); // faster than no_init_vector ?? -> good, cannot divide by interger 0!!, also else numerically unstable and no speed loss !!
//         for(int i = row; i--; ) {
//           if(!std::isnan(column[i])) { // faster way to code this ??? -> Not Bad at all, 54.. millisec for WDIM
//             if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//             else { 
//               sumj[g[i]-1] += column[i];
//               ++nj[g[i]-1];
//             }
//           }
//         }
//         for(int i = ng; i--; ) sumj[i] /= nj[i];
//         sum[j] = sumj;
//       }
//     } else {
//       if(gs.size() != ng) stop("Vector of group-sizes must match number of groups"); 
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng); //  = no_init_vector //  Not initializing seems to be numerically unstable !!!!
//         int ngs = 0;
//         for(int i = 0; i != row; ++i) {
//           if(std::isnan(column[i])) {
//             if(!std::isnan(sumj[g[i]-1])) {
//               sumj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           } else {
//             sumj[g[i]-1] += column[i];
//           }
//         }
//         for(int i = ng; i--; ) sumj[i] /= gs[i];
//         sum[j] = sumj;
//       }
//     }
//     return sum;
//   }
// }

// Previous Version: With return and before efficient algorithms. 
// // [[Rcpp::export]]
// List fgmeanlCpp(List x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, 
//                 bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.size(); 
//   List sum(l); 
//   
//   if (ng == 0) { // No groups !!
//     if(narm) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double sumi = 0;
//           int ni = 0;
//           int row = column.size();
//           for(int i = row; i--; ) {
//             if(std::isnan(column[i])) continue; // faster way??
//             sumi += column[i];
//             ni++;
//           }
//           sumi /= ni;
//           if(ret == 1) {
//             if(fill) sum[j] = rep(sumi, row);
//             else {
//               NumericVector sgj(row, sumi);
//               for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i]; // 2-times vector subsetting faster??
//               sum[j] = sgj;
//             }
//           } else if (ret == 2) {
//             if(fill) sum[j] = column - sumi; 
//             else {
//               NumericVector sgj(row);
//               for(int i = row; i--; ) {
//                 if(std::isnan(column[i])) sgj[i] = column[i]; // 2-times vector subsetting faster??
//                 else sgj[i] = column[i] - sumi; 
//               }
//               sum[j] = sgj;
//             }
//           } else sum[j] = sumi;
//         }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double sumi = 0;
//         int row = column.size();
//         for(int i = row; i--; ) sumi += column[i];
//         sumi /= row;
//         if(ret == 1) { 
//           sum[j] = rep(sumi, row);
//         } else if (ret == 2) {
//           sum[j] = column - sumi;
//         } else sum[j] = sumi;
//       }
//     }
//   } else { // With groups !!
//     if(narm) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           IntegerVector n(ng);
//           for(int i = row; i--; ) {
//             if(std::isnan(column[i])) continue; // faster way??
//             sumj[g[i]-1] += column[i];  
//             n[g[i]-1]++;
//           }
//           for(int i = ng; i--; ) sumj[i] /= n[i]; 
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             if(fill) for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? 
//             else for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i]; 
//               else sgj[i] = sumj[g[i]-1];
//             }
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row); // This doesn't work vectorized!!
//             if(fill) for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1]; // Fastest way?? 
//             else for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i]; 
//               else sgj[i] = column[i] - sumj[g[i]-1];
//             }
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//     } else { 
//       if(gs.size() == 1) { 
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           IntegerVector n(ng);
//           for(int i = row; i--; ) {
//             sumj[g[i]-1] += column[i]; 
//             n[g[i]-1]++;
//           }
//           for(int i = ng; i--; ) sumj[i] /= n[i];
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? 
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1];
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           for(int i = row; i--; ) sumj[g[i]-1] += column[i];  
//           for(int i = ng; i--; ) sumj[i] /= gs[i];
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? 
//             sum[j] = sgj; 
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1];
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       }
//     }
//   }
//   return sum;
// }

// // // [[Rcpp::plugins(cpp11)]]
// // #include <numeric>
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // Note: it seems that changing the g[i]-1 doesn't make things faster!!!
// // Also: Having gs makes things a lot faster!!
// // Still do boolvector regarding NA's if return == 1 or 2


// Older version: Found out that two times std::isnan is faster than locical vector is fill = false, 
// the above makes the code more parsimonious although slightly slower that this one !!!
// 
// // [[Rcpp::export]]
// List fgmeanlcpp(List x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, bool narm = true, int ret = 0, bool fill = false) { // void
//   int l = x.size();
//   List sum(l);
//   
//   if (ng == 0) { // No groups !!
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double sumi = 0;
//           int ni = 0;
//           int row = column.size();
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             sumi += column[i];
//             ni++;
//           }
//           sumi /= ni;
//           if(ret == 1) { 
//             sum[j] = rep(sumi, row);
//           } else if (ret == 2) {
//             sum[j] = column - sumi; 
//           } else sum[j] = sumi;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double sumi = 0;
//           int ni = 0;
//           int row = column.size();
//           // LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             // isnan[i] = std::isnan(column[i]); // ISNAN
//             // if(isnan[i]) continue; // faster way?? // 50 milliseconds
//             if(std::isnan(column[i])) continue; // 21 milliseconds -> faster than array!!
//             sumi += column[i];
//             ni++;
//           }
//           sumi /= ni;
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row, sumi); // Is this really more efficient that what you had before??
//             // for(int i = row; i--; ) if(isnan[i]) sgj[i] = column[i]; // 125 milliseconds // faster way?? vector subsetting??
//             for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i]; // 89 milliseconds -> faster than array!!
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               //if(isnan[i]) sgj[i] = column[i]; // 122 milliseconds // faster way?? vector subsetting??
//               if(std::isnan(column[i])) sgj[i] = column[i]; // 93 milliseconds -> faster than array!!
//               else sgj[i] = column[i] - sumi; 
//             }
//             sum[j] = sgj;
//           } else sum[j] = sumi;
//         }        
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double sumi = 0;
//         int row = column.size();
//         for(int i = row; i--; ) sumi += column[i];
//         sumi /= row;
//         if(ret == 1) { // possibility of reducing the number of passes??
//           sum[j] = rep(sumi, row);
//         } else if (ret == 2) {
//           sum[j] = column - sumi;
//         } else sum[j] = sumi;
//       }
//     }
//   } else { // With groups !!
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           IntegerVector n(ng);
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//             n[g[i]-1]++;
//           }
//           for(int i = ng; i--; ) sumj[i] /= n[i]; // necessary because n is integer!!
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? better loop through cols??
//             sum[j] = sgj;
//             // sum[j] = sumj[g-1]; // This is slower!!! DOn't use!!
//           } else if (ret == 2) {
//             NumericVector sgj(row); // This doesn't work vectorized!!
//             for(int i = row; i--; )  sgj[i] = column[i] - sumj[g[i]-1];
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           IntegerVector n(ng);
//           // LogicalVector isnan(row); // 117 milliseconds!! -> almost twice as fast as array -> perhaps because of taking one column at a time -> more memory efficient??
//           for(int i = row; i--; ) {
//             // isnan[i] = std::isnan(column[i]);
//             // if(isnan[i]) continue; // faster way??
//             if(std::isnan(column[i])) continue; // 83 milliseconds -> more than 2 times array speed !!
//             sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//             n[g[i]-1]++;
//           }
//           for(int i = ng; i--; ) sumj[i] /= n[i];
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               // if(isnan[i]) sgj[i] = column[i];  // 218 milliseconds -> twice as fast as array !! // faster way?? vector subsetting??
//               if(std::isnan(column[i])) sgj[i] = column[i]; // 181 milliseconds -> faster!!
//               else sgj[i] = sumj[g[i]-1];
//             }
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               // if(isnan[i]) sgj[i] = column[i]; // 220 milliseconds -> twice as fast as array !! // faster way?? vector subsetting??
//               if(std::isnan(column[i])) sgj[i] = column[i]; // 180 milliseconds -> faster !!
//               else sgj[i] = column[i] - sumj[g[i]-1]; 
//             }
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       }
//     } else { // could make this faster by also supplying len!!! ie n = len
//       if(gs.size() == 1) { 
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           IntegerVector n(ng);
//           for(int i = row; i--; ) {
//             sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//             n[g[i]-1]++;
//           }
//           for(int i = ng; i--; ) sumj[i] /= n[i];
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? better loop through cols??
//             sum[j] = sgj;
//             // sum[j] = sumj[g-1]; This is slower!!!
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - sumj[g[i]-1];
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           for(int i = row; i--; ) sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//           for(int i = ng; i--; ) sumj[i] /= gs[i];
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? better loop through cols??
//             sum[j] = sgj;
//             // sum[j] = sumj[g-1]; This is slower!!
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - sumj[g[i]-1];
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       }
//     }
//   }
//   return sum;
// }

// Older version without some efficient vectorizations:
// // [[Rcpp::export]]
// List fgmeanlcppold(List x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, bool narm = true, int ret = 0, bool fill = false) { // void
//   int l = x.size();
//   List sum(l);
//   
//   if (ng == 0) { // No groups
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double sumi = 0;
//           int ni = 0;
//           int row = column.size();
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             sumi += column[i];
//             ni++;
//           }
//           sumi /= ni;
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row); // try std::fill ?? 
//             for(int i = row; i--; )  sgj[i] = sumi; // faster with loop??: for(int i = row; i--; ) sgj[i] = sumj[g[i]-1];
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - sumi;
//             sum[j] = sgj;
//           } else sum[j] = sumi;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double sumi = 0;
//           int ni = 0;
//           int row = column.size();
//           LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             isnan[i] = ISNAN(column[i]);
//             if(isnan[i]) continue; // faster way??
//             sumi += column[i];
//             ni++;
//           }
//           sumi /= ni;
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = sumi; 
//             }
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = column[i] - sumi; 
//             }
//             sum[j] = sgj;
//           } else sum[j] = sumi;
//         }        
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double sumi = 0;
//         int row = column.size();
//         for(int i = row; i--; ) sumi += column[i];
//         sumi /= row;
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = sumi; // faster with loop??: for(int i = row; i--; ) sgj[i] = sumj[g[i]-1];
//           sum[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = column[i] - sumi;
//           sum[j] = sgj;
//         } else sum[j] = sumi;
//       }
//     }
//   } else { // With groups
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           IntegerVector n(ng);
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//             n[g[i]-1]++;
//           }
//           for(int i = ng; i--; ) sumj[i] /= n[i];
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? better loop through cols??
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - sumj[g[i]-1];
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           IntegerVector n(ng);
//           LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             isnan[i] = ISNAN(column[i]);
//             if(isnan[i]) continue; // faster way??
//             sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//             n[g[i]-1]++;
//           }
//           for(int i = ng; i--; ) sumj[i] /= n[i];
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = sumj[g[i]-1]; 
//             }
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = column[i] - sumj[g[i]-1]; 
//             }
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       }
//     } else { // could make this faster by also supplying len!!! ie n = len
//       if(gs.size() == 1) { 
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           IntegerVector n(ng);
//           for(int i = row; i--; ) {
//             sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//             n[g[i]-1]++;
//           }
//           for(int i = ng; i--; ) sumj[i] /= n[i];
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? better loop through cols??
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - sumj[g[i]-1];
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           for(int i = row; i--; ) sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//           for(int i = ng; i--; ) sumj[i] /= gs[i];
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? better loop through cols??
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - sumj[g[i]-1];
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       }
//     }
//   }
//   return sum;
// }


// This version without fill--> same speed as fill = TRUE, fill = FALSE slightly slower. 
// // [[Rcpp::export]]
// List fgmeanlcpp2(List x, IntegerVector g = 0, int ng = 0, bool narm = false, int ret = 0, IntegerVector gs = 0) { // void
//   int l = x.size();
//   List sum(l);
//   
//   if (ng == 0) {
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double sumi = 0;
//         int ni = 0;
//         int row = column.size();
//         for(int i = row; i--; ) {
//           if(ISNAN(column[i])) continue; // faster way??
//           sumi += column[i];
//           ni++;
//         }
//         sum[j] = sumi/ni;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double sumi = 0;
//         int row = column.size();
//         for(int i = row; i--; ) {
//           sumi += column[i];
//         }
//         sum[j] = sumi/row;
//       }
//     }
//   } else {
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         NumericVector sumj(ng);
//         IntegerVector n(ng);
//         for(int i = row; i--; ) {
//           if(ISNAN(column[i])) continue; // faster way??
//           sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//           n[g[i]-1]++;
//         }
//         for(int i = ng; i--; ) { 
//           sumj[i] /= n[i];
//         }
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? better loop through cols??
//           sum[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = column[i] - sumj[g[i]-1];
//           sum[j] = sgj;
//         } else sum[j] = sumj;
//       }
//     } else { // could make this faster by also supplying len!!! ie n = len
//       if(gs.size() == 1) { 
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           IntegerVector n(ng);
//           for(int i = row; i--; ) {
//             sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//             n[g[i]-1]++;
//           }
//           for(int i = ng; i--; ) { 
//             sumj[i] /= n[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? better loop through cols??
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - sumj[g[i]-1];
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector sumj(ng);
//           for(int i = row; i--; ) {
//             sumj[g[i]-1] += column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//           }
//           for(int i = ng; i--; ) { 
//             sumj[i] /= gs[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; // Fastest way?? better loop through cols??
//             sum[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - sumj[g[i]-1];
//             sum[j] = sgj;
//           } else sum[j] = sumj;
//         }
//       }
//     }
//   }
//   return sum;
// }
