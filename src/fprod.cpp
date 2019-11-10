#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fprodCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0, 
                       bool narm = true) {
  int l = x.size();
  
  if(ng == 0) {
    if(narm) {
      int j = l-1;
      long double prod = x[j]; 
      while(std::isnan(prod) && j!=0) prod = x[--j]; 
      if(j != 0) for(int i = j; i--; ) {
        if(!std::isnan(x[i])) prod *= x[i]; // Fastest ??
      } 
      return NumericVector::create((double)prod);
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
      return NumericVector::create((double)prod);
    }
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    if(narm) {
      NumericVector prod(ng, NA_REAL); // Other way ??
      for(int i = l; i--; ) {
        if(!std::isnan(x[i])) { // faster way to code this ??? -> Not Bad at all
          if(std::isnan(prod[g[i]-1])) prod[g[i]-1] = x[i];
          else prod[g[i]-1] *= x[i];
        }
      }
      DUPLICATE_ATTRIB(prod, x);
      return prod;
    } else {
      NumericVector prod(ng, 1.0); // good?? -> yes
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
      DUPLICATE_ATTRIB(prod, x);
      return prod;
    }
  }
}



// Previous Version: Not efficient lagorithm and with return
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// NumericVector fgprodCpp(NumericVector x, int ng = 0, IntegerVector g = 0, 
//                        bool narm = true, int ret = 0, bool fill = false) {
//   int l = x.size();
//   if (ng == 0) {
//     double prod = 1.0;
//     if(narm) {
//       for(int i = l; i--; ) {
//         if(std::isnan(x[i])) continue; 
//         prod *= x[i];
//       }
//       if(ret == 1) { 
//         if(fill) return rep(prod, l);
//         else {
//           NumericVector out(l, prod); 
//           for(int i = l; i--; ) if(std::isnan(x[i])) out[i] = x[i]; 
//           return out;
//         }
//       } else if (ret == 2) {
//         if(fill) return x - prod; 
//         else {
//           NumericVector out(l);
//           for(int i = l; i--; ) { 
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = x[i] - prod; 
//           }
//           return out;
//         }
//       } else return NumericVector::create(prod);
//     } else {
//       for(int i = l; i--; ) prod *= x[i];
//       if(ret == 1) { 
//         return rep(prod, l);
//       } else if (ret == 2) {
//         return x - prod; 
//       } else return NumericVector::create(prod);
//     }
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericVector prod(ng, 1.0);
//     if(narm) {
//       for(int i = l; i--; ) {
//         if(std::isnan(x[i])) continue; 
//         prod[g[i]-1] *= x[i];
//       }
//       if(ret == 1) { 
//         NumericVector out(l);
//         if(fill) {
//           for(int i = l; i--; ) out[i] = prod[g[i]-1];
//           return out;
//         } else {
//           for(int i = l; i--; ) {
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = prod[g[i]-1]; 
//           }
//           return out; 
//         }
//       } else if (ret == 2) {
//         NumericVector out(l);
//         if(fill) {
//           for(int i = l; i--; ) out[i] = x[i] - prod[g[i]-1];
//           return out;
//         } else {
//           for(int i = l; i--; ) {
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = x[i] - prod[g[i]-1]; 
//           }
//           return out; 
//         }
//       } else return prod; 
//     } else {
//       for(int i = l; i--; ) prod[g[i]-1] *= x[i];
//       if(ret == 1) { 
//         NumericVector out(l);
//         for(int i = l; i--; ) out[i] = prod[g[i]-1];
//         return out; 
//       } else if (ret == 2) {
//         NumericVector out(l);
//         for(int i = l; i--; ) out[i] = x[i] - prod[g[i]-1]; 
//         return out; 
//       } else return prod; 
//     }
//   }
// }
