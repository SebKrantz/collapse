#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fminCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0, 
                      bool narm = true) {
  int l = x.size();
  
  if(ng == 0) {
    if(narm) {
      int j = l-1;
      double min = x[j]; 
      while(std::isnan(min) && j!=0) min = x[--j]; 
      if(j != 0) for(int i = j; i--; ) {
        if(min > x[i]) min = x[i]; 
      } 
      return NumericVector::create(min);
    } else {
      double min = x[0];
      for(int i = 0; i != l; ++i) { 
        if(std::isnan(x[i])) {
          min = x[i]; 
          break;
        } else { 
          if(min > x[i]) min = x[i];
        }
      }
      return NumericVector::create(min);
    }
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    if(narm) {
      NumericVector min(ng, NA_REAL); // Other way ??
      for(int i = l; i--; ) { // adding if isnan(x[i]) before is not faster !!!
          if(min[g[i]-1] > x[i] || std::isnan(min[g[i]-1])) min[g[i]-1] = x[i];  // fastest !!
      }
      DUPLICATE_ATTRIB(min, x);
      return min;
    } else {
      NumericVector min(ng, R_PosInf); // INFINITY // DBL_MAX // good?? -> yes, same value bu better output !!
      int ngs = 0;
      for(int i = 0; i != l; ++i) {
        if(std::isnan(x[i])) {
          if(!std::isnan(min[g[i]-1])) {
            min[g[i]-1] = x[i];
            ++ngs;
            if(ngs == ng) break;
          }
        } else {
          if(min[g[i]-1] > x[i]) min[g[i]-1] = x[i];
        }
      }
      DUPLICATE_ATTRIB(min, x);
      return min;
    }
  }
}


// Previous Version: Lots of comments and attempts of finding the most efficient algorithms implemented above !! with return !!
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// NumericVector fgminCpp(NumericVector x, int ng = 0, IntegerVector g = 0, 
//                         bool narm = true, int ret = 0, bool fill = false) {
//   int l = x.size();
//   if (ng == 0) {
//     double min = x[l-1]; // R_PosInf; // All operations involving NA or NaN evaluate as false!!
//     if(narm) {
//       int j = l-1; // right??
//       while(std::isnan(min) && j!=0) { // fastest so far!! j!=0 faster than j>0
//          // --j; // yes !! --j faster than --j ?? // https://stackoverflow.com/questions/21447332/increment-i-i-and-i-1
//          min = x[--j]; // https://stackoverflow.com/questions/24901/is-there-a-performance-difference-between-i-and-i-in-c?noredirect=1&lq=1
//       }
//       // // bool t = true;
//       // double inf = std::numeric_limits<double>::infinity();
//       if(j != 0) for(int i = j; i--; ) {
//         // if(std::isnan(x[i])) continue; 
//         // // t = min < inf;
//         // if(min > x[i] || !(min > -INFINITY)) min = x[i]; // INFINITY // min > x[i] || !(min < inf)
//         // min > x[i] || !(min > inf) 7 milliseconds on 10 mio obs, but this statement is always true
//         // if(min > x[i] || std::isnan(min)) min = x[i]; // best so far: 28 milliseconds
//         if(min > x[i]) min = x[i];
//         // if(x[i] <= INFINITY) min += x[i]; // 21 milliseconds for the sum !!
//       } // evalCpp("NAN != INFINITY") consider this !! always TRUE
//       if(ret == 1) { 
//         if(fill) return rep(min, l);
//         else {
//           NumericVector out(l, min); 
//           for(int i = l; i--; ) if(std::isnan(x[i])) out[i] = x[i]; 
//           return out;
//         }
//       } else if (ret == 2) {
//         if(fill) return x - min; 
//         else {
//           NumericVector out(l);
//           for(int i = l; i--; ) { 
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = x[i] - min; 
//           }
//           return out;
//         }
//       } else return NumericVector::create(min);
//     } else {
//       for(int i = 0; i != l; ++i) { // better since it seems most series have missing values at the beginning??
//         if(std::isnan(x[i])) {
//           min = x[i]; // Fast ?? -> Yup, great !!
//           break;
//         } else {// else if ?? -> nope, this is faster !!! -> good to know !!
//           if(min > x[i]) min = x[i];
//         }
//       }
//       if(ret == 1) { 
//         return rep(min, l);
//       } else if (ret == 2) {
//         return x - min; 
//       } else return NumericVector::create(min);
//     } 
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     // NumericVector min(ng, NA_REAL); //INFINITY); //R_PosInf); // INFINITY is slightly faster 
//     if(narm) {
//       NumericVector min(ng, NA_REAL); //INFINITY); //R_PosInf); // INFINITY is slightly faster 
//       // int group = 0; // These don't improve speed!!
//       // double val = 0;
//       // int ngs = 0;
//       for(int i = l; i--; ) {
//         // if(std::isnan(x[i])) continue; // Not necessary -> 10 milliseconds gain
//         // group = g[i]-1;
//         // val = x[i];
//         // val = min[g[i]-1];
//         // if(min[group] > val) min[group] = val; 
//         if(min[g[i]-1] > x[i] || std::isnan(min[g[i]-1])) min[g[i]-1] = x[i];  // fastest !!!!!!
//         // if(std::isnan(min[g[i]-1]) && !std::isnan(x[i])) { Neat Idea but takes too long
//         //   min[g[i]-1] = x[i];
//         //   ++ngs;
//         //   if(ngs == ng) {
//         //     for(int j = i; j--; ) {
//         //       if(min[g[j]-1] > x[j]) min[g[j]-1] = x[j];
//         //     }
//         //     break;
//         //   }
//         // }
//         // if(min[g[i]-1] > x[i]) min[g[i]-1] = x[i];
//       }
//       if(ret == 1) { 
//         NumericVector out(l);
//         if(fill) {
//           for(int i = l; i--; ) out[i] = min[g[i]-1];
//           return out;
//         } else {
//           for(int i = l; i--; ) {
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = min[g[i]-1]; 
//           }
//           return out; 
//         }
//       } else if (ret == 2) {
//         NumericVector out(l);
//         if(fill) {
//           for(int i = l; i--; ) out[i] = x[i] - min[g[i]-1];
//           return out;
//         } else {
//           for(int i = l; i--; ) {
//             if(std::isnan(x[i])) out[i] = x[i]; 
//             else out[i] = x[i] - min[g[i]-1]; 
//           }
//           return out; 
//         }
//       } else return min; 
//     } else {
//       NumericVector min(ng, INFINITY); //R_PosInf); // INFINITY is slightly faster 
//       //for(int i = l; i--; ) if(min[g[i]-1] > x[i]) min[g[i]-1] = x[i];
//       int ngs = 0;
//       for(int i = 0; i != l; ++i) { // better since it seems most series have missing values at the beginning??
//         // if(std::isnan(x[i]) && !std::isnan(min[g[i]-1])) { // Basic idea: 135 milliseconds
//         //   min[g[i]-1] = x[i]; // Fast ?? -> Yup, great !!
//         //   ++ngs;
//         //   if(ngs == ng) break;
//         // } else {
//         //     if(min[g[i]-1] > x[i]) min[g[i]-1] = x[i]; // else faster?? 
//         // }
//         
//         // if(!std::isnan(min[g[i]-1])) { // 188 milliseconds -> slower !!
//         //   if(std::isnan(x[i])) {
//         //     min[g[i]-1] = x[i]; // Fast ?? -> Yup, great !!
//         //     ++ngs;
//         //     if(ngs == ng) break;
//         //   } else {
//         //      if(min[g[i]-1] > x[i]) min[g[i]-1] = x[i];
//         //   }
//         // }
//         
//         if(std::isnan(x[i])) { // 135 milliseconds -> same as above!! still compare speed on other data!!
//           if(!std::isnan(min[g[i]-1])) {
//             min[g[i]-1] = x[i]; // Fast ?? -> Yup, great !!
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         } else { // This is the best !!! not else if
//          if(min[g[i]-1] > x[i]) min[g[i]-1] = x[i];
//         }
//       }
//       if(ret == 1) { 
//         NumericVector out(l);
//         for(int i = l; i--; ) out[i] = min[g[i]-1];
//         return out; 
//       } else if (ret == 2) {
//         NumericVector out(l);
//         for(int i = l; i--; ) out[i] = x[i] - min[g[i]-1]; 
//         return out; 
//       } else return min; 
//     }
//   }
// }
