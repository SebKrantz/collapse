#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fmeanmCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0, const SEXP& gs = R_NilValue,  
               const SEXP& w = R_NilValue, bool narm = true, bool drop = true) { 
  int l = x.nrow(), col = x.ncol(); 
  
 if(Rf_isNull(w)) { // No weights !! 
  if(ng == 0) { 
    NumericVector sum = no_init_vector(col); //  Initialize faster -> Nope !!!
    if(narm) { 
      for(int j = col; j--; ) { // Instead Am(j,_) you can use Am.row(j).
        NumericMatrix::ConstColumn column = x( _ , j); 
        int k = l-1, nj = 1;
        long double sumj = column[k]; 
        while(std::isnan(sumj) && k!=0) sumj = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(std::isnan(column[i])) continue;
          sumj += column[i]; 
          ++nj;
        }
        sumj = sumj/nj;
        sum[j] = (double)sumj;
      }
    } else {
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
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
        sum[j] = (double)sumj;
      }
    }
    if(drop) sum.attr("names") = colnames(x); // Slight speed loss 31 to 34 milliseconds on WDIM, but doing it in R not faster !!
    else {
      sum.attr("dim") = Dimension(1, col);
      // sum.attr("dimnames") = List::create(R_NilValue,colnames(x)); 
      colnames(sum) = colnames(x); // NEW! faster than R ?? -> yes, good !!!
    }
    return sum;
  } else { // with groups 
    if(g.size() != l) stop("length(g) must match nrow(X)");
    if(narm) {
      NumericMatrix sum = no_init_matrix(ng, col);
      std::fill(sum.begin(), sum.end(), NA_REAL); // fastest ?? or create vector and declare as matrix ??
      // NumericVector sumt(ng*col, NA_REAL); // A tiny speed gain, but not much !! Same memory efficiency !!
      // sumt.attr("dim") = Dimension(ng, col);
      // NumericMatrix sum = as<NumericMatrix>(sumt);
      for(int j = col; j--; ) { 
        NumericMatrix::ConstColumn column = x( _ , j); 
        NumericMatrix::Column sumj = sum( _ , j);
        int nj[ng]; // Numerically stable and faster and more memory efficient than before !!
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
        for(int i = ng; i--; ) sumj[i] /= nj[i];
      }
      colnames(sum) = colnames(x);  // extremely efficient !!
      return sum;
    } else {
      NumericMatrix sum(ng, col); // no init numerically unstable !!!
      if(Rf_isNull(gs)) {
        int gsv[ng], memsize = sizeof(int)*ng;
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j); 
          NumericMatrix::Column sumj = sum( _ , j);
          memset(gsv, 0, memsize); // still a tiny bit faster than std::vector, but both have the same memory efficiency !!
          // std::vector<int> gsv(ng);
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
          for(int i = ng; i--; ) sumj[i] /= gsv[i];
        }      
      } else {
        IntegerVector gsv = gs;
        if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups"); 
          for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j); 
          NumericMatrix::Column sumj = sum( _ , j);
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
          for(int i = ng; i--; ) sumj[i] /= gsv[i];
        }
      }
      colnames(sum) = colnames(x);  // extremely efficient !!
      return sum;
    }
  }
 } else { // With weights
   NumericVector wg = w;
   if(l != wg.size()) stop("length(w) must match nrow(X)");
   if(ng == 0) { 
     NumericVector sum = no_init_vector(col); // Initialize faster -> Nope !!!
     if(narm) { 
       for(int j = col; j--; ) { // Instead Am(j,_) you can use Am.row(j).
         NumericMatrix::ConstColumn column = x( _ , j); 
         int k = l-1;
         while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k; 
         long double sumj = column[k]*wg[k], sumwj = wg[k];
         if(k != 0) for(int i = k; i--; ) {
           if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
           sumj += column[i]*wg[i];
           sumwj += wg[i];
         }
         sumj = sumj/sumwj;
         sum[j] = (double)sumj;
       }
     } else {
       for(int j = col; j--; ) {
         NumericMatrix::ConstColumn column = x( _ , j);
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
         sum[j] = (double)sumj;
       }
     }
     if(drop) sum.attr("names") = colnames(x); // Slight speed loss 31 to 34 milliseconds on WDIM, but doing it in R not faster !!
     else {
       sum.attr("dim") = Dimension(1, col);
       // sum.attr("dimnames") = List::create(R_NilValue,colnames(x)); 
       colnames(sum) = colnames(x); // NEW! faster than R ?? -> yes, good !!!
     }
     return sum;
   } else { // with groups 
     if(g.size() != l) stop("length(g) must match nrow(X)");
     if(narm) {
       NumericMatrix sum = no_init_matrix(ng, col);
       std::fill(sum.begin(), sum.end(), NA_REAL); 
       // NumericMatrix sumw = no_init_matrix(ng, col); // Numerically stable ??????? -> Yes !!
       for(int j = col; j--; ) { 
         NumericMatrix::ConstColumn column = x( _ , j); 
         NumericMatrix::Column sumj = sum( _ , j);
         // NumericMatrix::Column sumwj = sumw( _ , j);
         double sumwj[ng]; // Numerically stable, Slightly faster and a lot more memory efficient!! (but long double is a lot slower)
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
        for(int i = ng; i--; ) sumj[i] /= sumwj[i];
        // sumj = sumj/sumwj;
       }
       colnames(sum) = colnames(x);  // extremely efficient !!
       return sum;
     } else {
       NumericMatrix sum(ng, col); // no init numerically unstable !!!
       // NumericMatrix sumw(ng, col); // also here ?? -> Nope
       double sumwj[ng]; // Also a bit faster and a lot more memory efficient !!
       int memsize = sizeof(double)*ng;
       for(int j = col; j--; ) {
         NumericMatrix::ConstColumn column = x( _ , j); 
         NumericMatrix::Column sumj = sum( _ , j);
         // NumericMatrix::Column sumwj = sumw( _ , j);
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
         for(int i = ng; i--; ) sumj[i] /= sumwj[i];
         // sumj = sumj/sumwj;
       }
       colnames(sum) = colnames(x);  // extremely efficient !!
       return sum;
     }
   }
 }
}


// // Final Version without weights !!!
// // [[Rcpp::export]]
// SEXP fgmeanmCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0,  
//                bool narm = true, bool drop = true) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector sum(col); //  = no_init_vector // Initialize faster -> Yes !!!
//     if(narm) { 
//       for(int j = col; j--; ) { // Instead Am(j,_) you can use Am.row(j).
//         NumericMatrix::Column column = x( _ , j); 
//         int k = l-1, nj = 1;
//         double sumj = column[k]; 
//         while(std::isnan(sumj) && k!=0) sumj = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(std::isnan(column[i])) continue;
//           sumj += column[i]; 
//           ++nj;
//         }
//         sum[j] = sumj/nj;
//       }
//     } else {
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         double sumj = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(column[i])) {
//             sumj = column[i]; 
//             break;
//           } else { 
//             sumj += column[i];
//           }
//         }
//         sum[j] = sumj/l;
//       }
//     }
//     if(drop) sum.attr("names") = colnames(x); // Slight speed loss 31 to 34 milliseconds on WDIM, but doing it in R not faster !!
//     else {
//       sum.attr("dim") = Dimension(1, col);
//       // sum.attr("dimnames") = List::create(R_NilValue,colnames(x)); 
//       colnames(sum) = colnames(x); // NEW! faster than R ?? -> yes, good !!!
//     }
//     return sum;
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     if(narm) {
//       NumericMatrix sum = no_init_matrix(ng, col);
//       std::fill(sum.begin(), sum.end(), NA_REAL); 
//       IntegerMatrix n = no_init_matrix(ng, col); // Numerically stable ??????? -> Yes !!!
//       for(int j = col; j--; ) { 
//         NumericMatrix::Column column = x( _ , j); 
//         NumericMatrix::Column sumj = sum( _ , j);
//         IntegerMatrix::Column nj = n( _ , j);
//         for(int i = l; i--; ) { 
//           if(!std::isnan(column[i])) { 
//             if(std::isnan(sumj[g[i]-1])) {
//               sumj[g[i]-1] = column[i];
//               nj[g[i]-1] = 1;
//             } else {
//               sumj[g[i]-1] += column[i];
//               ++nj[g[i]-1];
//             }
//           }
//         }
//         for(int i = ng; i--; ) sumj[i] /= nj[i];
//       }
//       return sum;
//     } else {
//       if(gs.size() != ng) stop("Vector of group-sizes must match number of groups"); 
//       NumericMatrix sum(ng, col); // no init numerically unstable !!!
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j); 
//         NumericMatrix::Column sumj = sum( _ , j);
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
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
//       }
//       return sum;
//     }
//   }
// }


// Previous (old) version: before introducing efficient algorithms !!
// // [[Rcpp::export]]
// SEXP fgmeanmCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, 
//                 bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { // No groups !!
//     NumericVector sum(col); 
//     if(narm) { 
//       IntegerVector n(col); 
//         for(int i = l; i--; ) { // Faster way??
//           for(int j = col; j--; ) { // x(i,j) == NA_REAL || R_NaN
//             if(std::isnan(x(i,j))) continue; // isnan(x(i,j))!=0
//             sum[j] += x(i,j); 
//             n[j]++; 
//           } 
//         } 
//         for(int j = col; j--; ) sum[j] /= n[j]; 
//         // for(int j = col; j--; ) { // This is definitely slower !! I guess it copies column!!
//         //   NumericVector column = x(_,j);
//         //   double sumi = 0;
//         //   int ni = 0;
//         //   for(int i = l; i--; ) {
//         //     if(std::isnan(column[i])) continue; // faster way??
//         //     sumi += column[i];
//         //     ni++;
//         //   }
//         //   sum[j] = sumi/ni;
//         // } 
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericMatrix out(l, col);
//           if(fill) for(int j = col; j--; ) out(_,j) = rep(sum[j], l); // Fastest way?? better loop through cols??
//           else for(int i = l; i--; ) { 
//             for(int j = col; j--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//               else out(i,j) = sum[j];
//             }
//           }
//           return out;
//         } else if (ret == 2) {
//           NumericMatrix out(l, col); // Double loop is just as fast!!
//           if(fill) for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j]; 
//           else for(int i = l; i--; ) { 
//             for(int j = col; j--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//               else out(i,j) = x(i,j) - sum[j];
//             }
//           }
//           return out;
//         } else return sum;
//     } else {
//       for(int i = l; i--; ) {
//         for(int j = col; j--; ) sum[j] += x(i,j); 
//       }
//       sum = sum/l; 
//       if(ret == 1) { // possibility of reducing the number of passes??
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = rep(sum[j], l);
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j];
//         return out;
//       } else return sum;
//     }
//   } else { // with groups !!
//     NumericMatrix sum(ng, col);
//     if(narm) {
//       IntegerMatrix n(ng, col);
//         for(int i = l; i--; ) {
//           for(int j = col; j--; ) {
//             if(std::isnan(x(i,j))) continue; // ISNAN
//             sum(g[i]-1,j) += x(i,j); 
//             n(g[i]-1,j)++;
//           }
//         }
//         for(int i = ng; i--; ) { 
//           for(int j = col; j--; ) sum(i,j) /= n(i,j); // fastest??? -> probably!!
//         }
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericMatrix out(l, col);
//           if(fill) for(int i = l; i--; ) { // faster?? -> Yes !!
//             for(int j = col; j--; ) out(i,j) = sum(g[i]-1,j);
//           }
//           else for(int i = l; i--; ) { 
//             for(int j = col; j--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//               else out(i,j) = sum(g[i]-1,j);
//             }
//           }
//           return out;
//         } else if (ret == 2) {
//           NumericMatrix out(l, col);
//           if(fill) for(int i = l; i--; ) { // faster?? -> Yes !!
//             for(int j = col; j--; ) out(i,j) = x(i,j) - sum(g[i]-1,j);
//           }
//           else for(int i = l; i--; ) { 
//             for(int j = col; j--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//               else out(i,j) = x(i,j) - sum(g[i]-1,j);
//             }
//           }
//           return out;
//         } else return sum;
//     } else {
//       if(gs.size() == 1) {
//         IntegerVector n(ng);
//         for(int i = l; i--; ) {
//           n[g[i]-1]++;
//           for(int j = col; j--; ) sum(g[i]-1,j) += x(i,j); 
//         }
//         for(int i = ng; i--; ) { 
//           for(int j = col; j--; ) sum(i,j) /= n[i]; // fastest?? -> Yes!!
//         }
//       } else {
//         for(int i = l; i--; ) {
//           for(int j = col; j--; ) sum(g[i]-1,j) += x(i,j); 
//         }
//         for(int i = ng; i--; ) { 
//           for(int j = col; j--; ) sum(i,j) /= gs[i]; // fastest?? -> Yes!!
//         }
//       }
//     if(ret == 1) { // possibility of reducing the number of passes??
//       NumericMatrix out(l, col);
//       for(int i = l; i--; ) { // faster?? -> yes!!
//        for(int j = col; j--; ) out(i,j) = sum(g[i]-1,j);
//       }
//       return out;
//     } else if (ret == 2) {
//       NumericMatrix out(l, col);
//       for(int i = l; i--; ) { // faster?? -> yes!!
//         for(int j = col; j--; ) out(i,j) = x(i,j) - sum(g[i]-1,j);
//       }
//       return out;
//     } else return sum;
//   }
//  }
// }
 
 
// Older version: Just shrank the code after realizing that two times std::isnan is faster !!
 
//  // // [[Rcpp::plugins(cpp11)]]
//  // #include <R.h>
//  // #include <math.h>
//  // #include <cmath>
//  // //[[Rcpp::depends(math)]]
//  //  // #include <stan/math.hpp> 
//  // #include <numeric>
// #include <Rcpp.h>
//  using namespace Rcpp;
//  // using namespace cmath;
//  
//  // Note: still try to use math.h isnan... but it also sais this is difficult to use in a package !!
//  
//  // [[Rcpp::export]]
//  SEXP fgmeanacpp(NumericMatrix x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, bool narm = true, int ret = 0, bool fill = false) { // void
//    int l = x.nrow(); //sizeof(x)/sizeof(x[0]); 
//    int col = x.ncol(); // sizeof(x)/l;
//    
//    if(ng == 0) { // No groups !!
//      NumericVector sum(col);
//      if(narm) {
//        IntegerVector n(col);
//        if(fill) {
//          for(int i = l; i--; ) { // Faster way??
//            for(int j = col; j--; ) { // x(i,j) == NA_REAL || R_NaN // This does not really work!!!
//              if(std::isnan(x(i,j))) continue; // if(std::isnan(x(i,j))) continue; isnan(x(i,j))!=0
//              sum[j] += x(i,j); // x[i][j];
//              n[j]++;
//            }
//          }
//          for(int j = col; j--; ) sum[j] /= n[j]; 
//          if(ret == 1) { // possibility of reducing the number of passes??
//            NumericMatrix out(l, col);
//            for(int j = col; j--; ) out(_,j) = rep(sum[j], l); // Fastest way?? better loop through cols??
//            return out;
//          } else if (ret == 2) {
//            NumericMatrix out(l, col); // Double loop is just as fast!!
//            for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j]; // for(int i = l; i--; )  out(i,_) = x(i,_) - sum; // vectorized??
//            return out;
//          } else return sum;
//        } else {
//          // LogicalMatrix isnan(l, col); // perhaps it's faster without?? -> Yes!!!
//          for(int i = l; i--; ) {
//            for(int j = col; j--; ) {
//              // isnan(i,j) = std::isnan(x(i,j));
//              // if(isnan(i,j)) continue;
//              if(std::isnan(x(i,j))) continue;
//              sum[j] += x(i,j); // x[i][j];
//              n[j]++;
//            }
//          }
//          for(int j = col; j--; ) sum[j] /= n[j]; // vectorized faster??
//          if(ret == 1) { // possibility of reducing the number of passes??
//            NumericMatrix out(l, col);
//            for(int i = l; i--; ) { // fastest way? Vectorized better??
//              for(int j = col; j--; ) {
//                // if(isnan(i,j)) out(i,j) = x(i,j);
//                if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//                else out(i,j) = sum[j];
//              }
//            }
//            return out;
//          } else if (ret == 2) {
//            NumericMatrix out(l, col);
//            for(int i = l; i--; ) { // fastest way? Vectorized better??
//              for(int j = col; j--; ) {
//                //if(isnan(i,j)) out(i,j) = x(i,j);
//                if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//                else out(i,j) = x(i,j) - sum[j];
//              }
//            }
//            return out;
//          } else return sum;
//        }
//      } else {
//        for(int i = l; i--; ) {
//          for(int j = col; j--; ) sum[j] += x(i,j); // x[i][j];
//        }
//        sum = sum/l; // fastest?? -> Yes, good
//        if(ret == 1) { // possibility of reducing the number of passes??
//          NumericMatrix out(l, col);
//          for(int j = col; j--; ) out(_,j) = rep(sum[j], l); // for(int i = l; i--; ) out(i,_) = sum; // Fastest way?? better loop through cols??
//          return out;
//        } else if (ret == 2) {
//          NumericMatrix out(l, col);
//          for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j]; // for(int i = l; i--; )  out(i,_) = x(i,_) - sum; // vectorized??
//          return out;
//        } else return sum;
//      }
//    } else { // with groups !!
//      NumericMatrix sum(ng, col);
//      if(narm) {
//        IntegerMatrix n(ng, col);
//        if(fill) {
//          for(int i = l; i--; ) {
//            for(int j = col; j--; ) {
//              if(std::isnan(x(i,j))) continue; // ISNAN
//              sum(g[i]-1,j) += x(i,j); // x[i][j];
//              n(g[i]-1,j)++;
//            }
//          }
//          for(int i = ng; i--; ) { 
//            for(int j = col; j--; ) {
//              sum(i,j) /= n(i,j); // fastest??? -> probably!!
//            }
//          }
//          if(ret == 1) { // possibility of reducing the number of passes??
//            NumericMatrix out(l, col);
//            // for(int i = l; i--; ) out(i,_) = sum(g[i]-1,_);
//            for(int i = l; i--; ) { // faster?? -> Yes !!
//              for(int j = col; j--; ) out(i,j) = sum(g[i]-1,j);
//            }
//            return out;
//          } else if (ret == 2) {
//            NumericMatrix out(l, col);
//            // for(int i = l; i--; )  out(i,_) = x(i,_) - sum(g[i]-1,_);
//            for(int i = l; i--; ) { // faster?? -> yes!!
//              for(int j = col; j--; ) out(i,j) = x(i,j) - sum(g[i]-1,j);
//            }
//            return out;
//          } else return sum;
//        } else {
//          // LogicalMatrix isnan(l, col); // With this and ret = 0 318 milliseconds else 224 milliseconds
//          for(int i = l; i--; ) {
//            for(int j = col; j--; ) {
//              // isnan(i,j) = std::isnan(x(i,j));
//              // if(isnan(i,j)) continue;
//              if(std::isnan(x(i,j))) continue;
//              sum(g[i]-1,j) += x(i,j); // x[i][j];
//              n(g[i]-1,j)++;
//            }
//          }
//          for(int i = ng; i--; ) { 
//            for(int j = col; j--; ) {
//              sum(i,j) /= n(i,j); // fastest???
//            }
//          }
//          if(ret == 1) { // possibility of reducing the number of passes??
//            NumericMatrix out(l, col);
//            for(int i = l; i--; ) { 
//              for(int j = col; j--; ) {
//                //if(isnan(i,j)) out(i,j) = x(i,j); // 499 milliseconds
//                if(std::isnan(x(i,j))) out(i,j) = x(i,j); // 407 milliseconds
//                else out(i,j) = sum(g[i]-1,j);
//              }
//            }
//            return out;
//          } else if (ret == 2) {
//            NumericMatrix out(l, col);
//            for(int i = l; i--; ) { 
//              for(int j = col; j--; ) {
//                // if(isnan(i,j)) out(i,j) = x(i,j); // 517 milliseconds
//                if(std::isnan(x(i,j))) out(i,j) = x(i,j); // 406 milliseconds
//                else out(i,j) = x(i,j) - sum(g[i]-1,j);
//              }
//            }
//            return out;
//          } else return sum;
//        }
//      } else {
//        if(gs.size() == 1) {
//          IntegerVector n(ng);
//          for(int i = l; i--; ) {
//            n[g[i]-1]++;
//            for(int j = col; j--; ) {
//              sum(g[i]-1,j) += x(i,j); // x[i][j];
//              // n(g[i]-1,j)++; // YOu actually only need a vector here!!!
//            }
//          }
//          for(int i = ng; i--; ) { 
//            // sum(i,_) /= n[i]; // fastest?? could do this!!
//            for(int j = col; j--; ) {
//              sum(i,j) /= n[i]; // fastest?? -> Yes!!
//            }
//          }
//        } else {
//          for(int i = l; i--; ) {
//            // sum(g[i]-1,_) = sum(g[i]-1,_) + x(i,_);
//            for(int j = col; j--; ) {
//              sum(g[i]-1,j) += x(i,j); // x[i][j];
//            }
//          }
//          for(int i = ng; i--; ) { 
//            // sum(i,_) = sum(i,_) / gs[i]; // fastest?? could do this!! -> Seems to be faster!!
//            for(int j = col; j--; ) {
//              sum(i,j) /= gs[i]; // fastest?? -> Yes!!
//            }
//          }
//        }
//        if(ret == 1) { // possibility of reducing the number of passes??
//          NumericMatrix out(l, col);
//          // for(int i = l; i--; ) out(i,_) = sum(g[i]-1,_); // Fastest way?? better loop through cols??
//          for(int i = l; i--; ) { // faster?? -> yes!!
//            for(int j = col; j--; ) out(i,j) = sum(g[i]-1,j);
//          }
//          return out;
//        } else if (ret == 2) {
//          NumericMatrix out(l, col);
//          //      for(int i = l; i--; )  out(i,_) = x(i,_) - sum(g[i]-1,_);
//          for(int i = l; i--; ) { // faster?? -> yes!!
//            for(int j = col; j--; ) out(i,j) = x(i,j) - sum(g[i]-1,j);
//          }
//          return out;
//        } else return sum;
//      }
//    }
//  }
 
 
 
 
 
 
 
