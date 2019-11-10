#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fminmCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0,  
              bool narm = true, bool drop = true) { 
  int l = x.nrow(), col = x.ncol(); 
  
  if(ng == 0) { 
    NumericVector min = no_init_vector(col); // Initialize faster -> Nope !!!
    if(narm) { 
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j); 
        int k = l-1;
        double minj = column[k]; 
        while(std::isnan(minj) && k!=0) minj = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(minj > column[i]) minj = column[i]; 
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
            if(minj > column[i]) minj = column[i];
          }
        }
        min[j] = minj;
      }
    }
    if(drop) min.attr("names") = colnames(x); 
    else {
      min.attr("dim") = Dimension(1, col);
      colnames(min) = colnames(x); 
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
          if(!std::isnan(column[i])) { // Keeping this is faster !!!!
            if(minj[g[i]-1] > column[i] || std::isnan(minj[g[i]-1])) minj[g[i]-1] = column[i];
          }
        }
      }
    } else {
      std::fill(min.begin(), min.end(), R_PosInf); 
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
            if(minj[g[i]-1] > column[i]) minj[g[i]-1] = column[i];
          }
        }
      }
    }
    colnames(min) = colnames(x);
    return min;
  }
}

// Previous version: wthout column-referencing and with return 
// // [[Rcpp::export]]
// SEXP fgminmCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector min(col); 
//     if(narm) { 
//       int k = 0;
//       for(int j = col; j--; ) {
//         k = l-1;
//         min[j] = x(k,j); // faster using a double instead of min[j]?? -> NOPE, NOT FASTER !!!
//         while(std::isnan(min[j]) && k!=0) min[j] = x(--k,j);
//         if(k != 0) for(int i = k; i--; ) {
//             if(min[j] > x(i,j)) min[j] = x(i,j);
//         }
//       }
//       // int k = 0;
//       // double minj = 0;
//       // for(int j = col; j--; ) {
//       //   k = l-1;
//       //   minj = x(k,j); // faster using a double instead of min[j]?? -> NOPE, NOT FASTER !!!
//       //   while(std::isnan(minj) && k!=0) minj = x(--k,j);
//       //   if(k != 0) for(int i = k; i--; ) {
//       //     if(minj > x(i,j)) minj = x(i,j);
//       //   }
//       //   min[j] = minj;
//       // }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) out(_,j) = rep(min[j], l); 
//           // double minj = 0; // above is faster !!!
//           // for(int j = col; j--; ) {
//           //   minj = min[j];
//           //   for(int i = l; i--; ) out(i,j) = minj;
//           // }
//         } else {
//           for(int j = col; j--; ) { 
//            for(int i = l; i--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//             else out(i,j) = min[j];
//            }
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col); 
//         if(fill) {
//           for(int j = col; j--; ) out(_,j) = x(_,j) - min[j]; 
//         } else {
//           for(int j = col; j--; ) { 
//            for(int i = l; i--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//             else out(i,j) = x(i,j) - min[j];
//            }
//           }
//         }
//         return out;
//       } else return min;
//     } else {
//       for(int j = col; j--; ) {
//         min[j] = x(0,j); // faster using a double instead of min[j]?? -> I Guess not !!!
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(x(i,j))) {
//             min[j] = x(i,j); 
//             break;
//           } else { 
//             if(min[j] > x(i,j)) min[j] = x(i,j);
//           }
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = rep(min[j], l);
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = x(_,j) - min[j];
//         return out;
//       } else return min;
//     }
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericMatrix min(ng, col);
//     if(narm) {
//      std::fill(min.begin(), min.end(), NA_REAL); // NAN -> This is not the Issue !!!
//     for(int j = col; j--; ) { // It was 740 milliseconds for numeric WDI matrix, now I reversed the loop nesting -> 80 milliseconds !!!!!!!!!!
//       for(int i = l; i--; ) { // Intuitively it makes sense to put the loop with the most entities on the inside
//           // Fastest ??
//           if(min(g[i]-1,j) > x(i,j) || std::isnan(min(g[i]-1,j))) min(g[i]-1,j) = x(i,j);
//       }
//     }
//     // Do in Steps !!!
//     // NumericVector xj(l); // -> Nope, a lot slower, as seen before !!!!
//     // for(int j = col; j--; ) { // It was 740 milliseconds for numeric WDI matrix, now I reversed the loop nesting -> 80 milliseconds !!!!!!!!!!
//     //   xj = x(_,j);
//     //   NumericVector minj(l, NA_REAL); 
//     //   for(int i = l; i--; ) { // Intuitively it makes sense to put the loop with the most entities on the inside
//     //     // Fastest ?? 
//     //     if(minj(g[i]-1) > xj(i) || std::isnan(minj(g[i]-1))) minj(g[i]-1) = xj(i); 
//     //   }
//     //   min(_,j) = minj;
//     // }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) out(i,j) = min(g[i]-1,j);
//           }
//         } else { 
//           for(int j = col; j--; ) { 
//            for(int i = l; i--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//             else out(i,j) = min(g[i]-1,j);
//            }
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) out(i,j) = x(i,j) - min(g[i]-1,j);
//           }
//         } else { 
//          for(int j = col; j--; ) { 
//           for(int i = l; i--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//             else out(i,j) = x(i,j) - min(g[i]-1,j);
//           }
//          }
//         }
//         return out;
//       } else return min;
//     } else {
//       std::fill(min.begin(), min.end(), INFINITY);
//       int ngs = 0;
//       for(int j = col; j--; ) {
//         ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(x(i,j))) { // fastest ??
//             if(!std::isnan(min(g[i]-1,j))) {
//               min(g[i]-1,j) = x(i,j); 
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           } else { 
//             if(min(g[i]-1,j) > x(i,j)) min(g[i]-1,j) = x(i,j);
//           }
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) { 
//           for(int i = l; i--; ) out(i,j) = min(g[i]-1,j);
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) { 
//           for(int i = l; i--; ) out(i,j) = x(i,j) - min(g[i]-1,j);
//         }
//         return out;
//       } else return min;
//     }
//   }
// }
