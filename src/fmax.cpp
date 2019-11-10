#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
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

// Previous version: With return !!
// // [[Rcpp::export]]
// NumericVector fgmaxCpp(NumericVector x, int ng = 0, IntegerVector g = 0, 
//                        bool narm = true, int ret = 0, bool fill = false) {
//   int l = x.size();
//   if (ng == 0) {
//     double max = x[l-1]; 
//     if(narm) {
//       int j = l-1; 
//       while(std::isnan(max) && j!=0) max = x[--j]; 
//       if(j != 0) for(int i = j; i--; ) {
//          if(max < x[i]) max = x[i];
//       } 
//       if(ret == 1) { 
//         if(fill) return rep(max, l);
//         else {
//           NumericVector out(l, max); 
//           for(int i = l; i--; ) if(std::isnan(x[i])) out[i] = x[i]; 
//           return out;
//         }
//       } else if (ret == 2) {
//         if(fill) return x - max; 
//         else {
//           NumericVector out(l);
//           for(int i = l; i--; ) { 
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = x[i] - max; 
//           }
//           return out;
//         }
//       } else return NumericVector::create(max);
//     } else {
//       for(int i = 0; i != l; ++i) { 
//         if(std::isnan(x[i])) {
//           max = x[i]; 
//           break;
//         } else { 
//           if(max < x[i]) max = x[i];
//         }
//       }
//       if(ret == 1) { 
//         return rep(max, l);
//       } else if (ret == 2) {
//         return x - max; 
//       } else return NumericVector::create(max);
//     } 
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     if(narm) {
//       NumericVector max(ng, NA_REAL); 
//       for(int i = l; i--; ) {
//         if(max[g[i]-1] < x[i] || std::isnan(max[g[i]-1])) max[g[i]-1] = x[i];
//       }
//       if(ret == 1) { 
//         NumericVector out(l);
//         if(fill) {
//           for(int i = l; i--; ) out[i] = max[g[i]-1];
//           return out;
//         } else {
//           for(int i = l; i--; ) {
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = max[g[i]-1]; 
//           }
//           return out; 
//         }
//       } else if (ret == 2) {
//         NumericVector out(l);
//         if(fill) {
//           for(int i = l; i--; ) out[i] = x[i] - max[g[i]-1];
//           return out;
//         } else {
//           for(int i = l; i--; ) {
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = x[i] - max[g[i]-1]; 
//           }
//           return out; 
//         }
//       } else return max; 
//     } else {
//       NumericVector max(ng, INFINITY);
//       int ngs = 0;
//       for(int i = 0; i != l; ++i) { 
//         if(std::isnan(x[i])) {
//           if(!std::isnan(max[g[i]-1])) {
//             max[g[i]-1] = x[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         } else {
//           if(max[g[i]-1] < x[i]) max[g[i]-1] = x[i];
//         }
//       }
//       if(ret == 1) { 
//         NumericVector out(l);
//         for(int i = l; i--; ) out[i] = max[g[i]-1];
//         return out; 
//       } else if (ret == 2) {
//         NumericVector out(l);
//         for(int i = l; i--; ) out[i] = x[i] - max[g[i]-1]; 
//         return out; 
//       } else return max; 
//     }
//   }
// }
