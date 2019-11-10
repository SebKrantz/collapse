#include <Rcpp.h>
using namespace Rcpp;

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

// Previous version: Not efficient algorithm and return !!
// // [[Rcpp::export]]
// SEXP fgmaxmCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector max(col); 
//     if(narm) { 
//       int k = 0;
//       for(int j = col; j--; ) {
//         k = l-1;
//         max[j] = x(k,j);
//         while(std::isnan(max[j]) && k!=0) max[j] = x(--k,j);
//         if(k != 0) for(int i = k; i--; ) {
//           if(max[j] < x(i,j)) max[j] = x(i,j);
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) out(_,j) = rep(max[j], l); 
//         } else {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//               else out(i,j) = max[j];
//             }
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col); 
//         if(fill) {
//           for(int j = col; j--; ) out(_,j) = x(_,j) - max[j]; 
//         } else {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//               else out(i,j) = x(i,j) - max[j];
//             }
//           }
//         }
//         return out;
//       } else return max;
//     } else {
//       for(int j = col; j--; ) {
//         max[j] = x(0,j); //  Necessary because if no NA's and 0 is the largest value -> mistake !!
//         for(int i = 0; i != l; ++i) { // also need to start from 0 beause x(0,j) could be NA !!
//           if(std::isnan(x(i,j))) {
//             max[j] = x(i,j); 
//             break;
//           } else { 
//             if(max[j] < x(i,j)) max[j] = x(i,j);
//           }
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = rep(max[j], l);
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = x(_,j) - max[j];
//         return out;
//       } else return max;
//     }
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericMatrix max(ng, col);
//     if(narm) {
//       std::fill(max.begin(), max.end(), NA_REAL); 
//       for(int j = col; j--; ) { 
//         for(int i = l; i--; ) { 
//           if(max(g[i]-1,j) < x(i,j) || std::isnan(max(g[i]-1,j))) max(g[i]-1,j) = x(i,j);
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) out(i,j) = max(g[i]-1,j);
//           }
//         } else { 
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//               else out(i,j) = max(g[i]-1,j);
//             }
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) out(i,j) = x(i,j) - max(g[i]-1,j);
//           }
//         } else { 
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//               else out(i,j) = x(i,j) - max(g[i]-1,j);
//             }
//           }
//         }
//         return out;
//       } else return max;
//     } else {
//       std::fill(max.begin(), max.end(), INFINITY); // better way ?? 
//       int ngs = 0;
//       for(int j = col; j--; ) {
//         ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(x(i,j))) { // fastest ??
//             if(!std::isnan(max(g[i]-1,j))) {
//               max(g[i]-1,j) = x(i,j); 
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           } else { 
//             if(max(g[i]-1,j) < x(i,j)) max(g[i]-1,j) = x(i,j);
//           }
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) { 
//           for(int i = l; i--; ) out(i,j) = max(g[i]-1,j);
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) { 
//           for(int i = l; i--; ) out(i,j) = x(i,j) - max(g[i]-1,j);
//         }
//         return out;
//       } else return max;
//     }
//   }
// }
