#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fsumCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0, bool narm = true) {
  int l = x.size();
  
  if(ng == 0) {
    if(narm) {
      int j = l-1;
      long double sum = x[j]; 
      while(std::isnan(sum) && j!=0) sum = x[--j]; 
      if(j != 0) for(int i = j; i--; ) {
        if(!std::isnan(x[i])) sum += x[i]; // Fastest ??
      } 
      return NumericVector::create((double)sum); // Converting long double directly to numeric vector is slow !!!
    } else {
      long double sum = 0;
      for(int i = 0; i != l; ++i) { 
        if(std::isnan(x[i])) {
          sum = x[i]; 
          break;
        } else { 
          sum += x[i];
        }
      }
      return NumericVector::create((double)sum);
    }
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    if(narm) {
      NumericVector sum(ng, NA_REAL); // Other way ??
      for(int i = l; i--; ) { 
        if(!std::isnan(x[i])) { // faster way to code this ??? -> Not Bad at all
          if(std::isnan(sum[g[i]-1])) sum[g[i]-1] = x[i];
          else sum[g[i]-1] += x[i];
        }
      }
      DUPLICATE_ATTRIB(sum, x);
      return sum;
    } else {
      NumericVector sum(ng); // good?? -> yes !! // Not initializing numerically unstable !!!
      int ngs = 0;
      for(int i = 0; i != l; ++i) { 
        if(std::isnan(x[i])) { 
          if(!std::isnan(sum[g[i]-1])) {
            sum[g[i]-1] = x[i];
            ++ngs;
            if(ngs == ng) break;
          }
        } else {
          sum[g[i]-1] += x[i];
        }
      }
      DUPLICATE_ATTRIB(sum, x);
      return sum;
    }
  }
}


// Version with internal return 
// // [[Rcpp::export]]
// NumericVector fgsumCpp(NumericVector x, int ng = 0, IntegerVector g = 0, 
//                        bool narm = true, int ret = 0, bool fill = false) {
//   int l = x.size();
//   if (ng == 0) {
//     if(narm) {
//       int j = l-1;
//       double sum = x[j]; 
//       while(std::isnan(sum) && j!=0) sum = x[--j]; 
//       if(j != 0) for(int i = j; i--; ) {
//         if(!std::isnan(x[i])) sum += x[i]; // Fastest ??
//       } 
//       if(ret == 1) { 
//         if(fill) return rep(sum, l);
//         else {
//           NumericVector out(l, sum); 
//           for(int i = l; i--; ) if(std::isnan(x[i])) out[i] = x[i]; 
//           return out;
//         }
//       } else if (ret == 2) {
//         if(fill) return x - sum; 
//         else {
//           NumericVector out = no_init_vector(l); // any speed loss from this ?? 
//           for(int i = l; i--; ) { 
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = x[i] - sum; 
//           }
//           return out;
//         }
//       } else return NumericVector::create(sum);
//     } else {
//       double sum = 0;
//       for(int i = 0; i != l; ++i) { 
//         if(std::isnan(x[i])) {
//           sum = x[i]; 
//           break;
//         } else { 
//           sum += x[i];
//         }
//       }
//       if(ret == 1) { 
//         return rep(sum, l);
//       } else if (ret == 2) {
//         return x - sum; 
//       } else return NumericVector::create(sum);
//     } 
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     if(narm) {
//       NumericVector sum(ng, NA_REAL); // Other way ??
//       for(int i = l; i--; ) { // Best ??
//         if(!std::isnan(x[i])) { // faster way to code this ??? -> Not Bad at all
//           if(std::isnan(sum[g[i]-1])) sum[g[i]-1] = x[i];
//           else sum[g[i]-1] += x[i];
//         }
//         //if(!std::isnan(x[i])) sum[g[i]-1] += x[i];
//       }
//       if(ret == 1) { 
//         NumericVector out = no_init_vector(l);
//         if(fill) {
//           for(int i = l; i--; ) out[i] = sum[g[i]-1];
//           return out;
//         } else {
//           for(int i = l; i--; ) {
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = sum[g[i]-1]; 
//           }
//           return out; 
//         }
//       } else if (ret == 2) {
//         NumericVector out = no_init_vector(l);
//         if(fill) {
//           for(int i = l; i--; ) out[i] = x[i] - sum[g[i]-1];
//           return out;
//         } else {
//           for(int i = l; i--; ) {
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = x[i] - sum[g[i]-1]; 
//           }
//           return out; 
//         }
//       } else return sum; 
//     } else {
//       NumericVector sum = no_init_vector(ng); // good?? -> yes !!
//       int ngs = 0;
//       for(int i = 0; i != l; ++i) { 
//         if(std::isnan(x[i])) { // maybe the first option is still the fastest i.e. checking both at the same time ?? -> test that !!
//           if(!std::isnan(sum[g[i]-1])) {
//             sum[g[i]-1] = x[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         } else {
//           sum[g[i]-1] += x[i];
//         }
//         // if(std::isnan(x[i]) && !std::isnan(sum[g[i]-1])) { // Doesn't really make any difference !!!!
//         //     sum[g[i]-1] = x[i];
//         //     ++ngs;
//         //     if(ngs == ng) break;
//         // } else {
//         //   sum[g[i]-1] += x[i];
//         // }
//       }
//       if(ret == 1) { 
//         NumericVector out = no_init_vector(l);
//         for(int i = l; i--; ) out[i] = sum[g[i]-1];
//         return out; 
//       } else if (ret == 2) {
//         NumericVector out = no_init_vector(l);
//         for(int i = l; i--; ) out[i] = x[i] - sum[g[i]-1]; 
//         return out; 
//       } else return sum; 
//     }
//   }
// }


// Previous Version
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// NumericVector fgsumCpp(NumericVector x, int ng = 0, IntegerVector g = 0, 
//                        bool narm = true, int ret = 0, bool fill = false) {
//   int l = x.size();
//   if (ng == 0) {
//     double sum = 0;
//     if(narm) {
//        // double inf = std::numeric_limits<double>::infinity();
//        // double v = 0;
//         for(int i = l; i--; ) {
//           // v = x[i];
//           // if(std::isnan(v)) continue; // also 22 millisconds !!
//           // sum += v;
//           if(std::isnan(x[i])) continue; // 22 milliseconds for 10 mio obs with 10% NA's
//           // if(!(x[i] <= INFINITY)) continue; // also 22 milliseconds !!
//           sum += x[i];
//           // if(x[i] <= INFINITY) sum += x[i]; // 21 milliseconds -> slightly faster
//         }
//       if(ret == 1) { 
//         if(fill) return rep(sum, l);
//         else {
//           NumericVector out(l, sum); 
//           for(int i = l; i--; ) if(std::isnan(x[i])) out[i] = x[i]; 
//           return out;
//         }
//       } else if (ret == 2) {
//         if(fill) return x - sum; 
//         else {
//           NumericVector out(l);
//           for(int i = l; i--; ) { 
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = x[i] - sum; 
//           }
//           return out;
//         }
//       } else return NumericVector::create(sum);
//     } else {
//         for(int i = l; i--; ) sum += x[i];
//       if(ret == 1) { 
//         return rep(sum, l);
//       } else if (ret == 2) {
//         return x - sum; 
//       } else return NumericVector::create(sum);
//     }
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericVector sum(ng);
//     if(narm) {
//         for(int i = l; i--; ) {
//           if(std::isnan(x[i])) continue; 
//           sum[g[i]-1] += x[i];
//         }
//       if(ret == 1) { 
//         NumericVector out(l);
//         if(fill) {
//           for(int i = l; i--; ) out[i] = sum[g[i]-1];
//           return out;
//         } else {
//           for(int i = l; i--; ) {
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = sum[g[i]-1]; 
//           }
//           return out; 
//         }
//       } else if (ret == 2) {
//         NumericVector out(l);
//         if(fill) {
//           for(int i = l; i--; ) out[i] = x[i] - sum[g[i]-1];
//           return out;
//         } else {
//           for(int i = l; i--; ) {
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = x[i] - sum[g[i]-1]; 
//           }
//           return out; 
//         }
//       } else return sum; 
//     } else {
//       for(int i = l; i--; ) sum[g[i]-1] += x[i];
//       if(ret == 1) { 
//         NumericVector out(l);
//         for(int i = l; i--; ) out[i] = sum[g[i]-1];
//         return out; 
//       } else if (ret == 2) {
//         NumericVector out(l);
//         for(int i = l; i--; ) out[i] = x[i] - sum[g[i]-1]; 
//         return out; 
//       } else return sum; 
//     }
//   }
// }
