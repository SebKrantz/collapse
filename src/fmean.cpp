#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fmeanCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0, 
                       const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, bool narm = true) {
  int l = x.size();
  
  if (Rf_isNull(w)) { // No weights !!
    if (ng == 0) {
      if(narm) {
        int j = l-1, n = 1; // 1 because for-loop starts from 2
        long double sum = x[j]; 
        while(std::isnan(sum) && j!=0) sum = x[--j]; 
        if(j != 0) for(int i = j; i--; ) {
          if(std::isnan(x[i])) continue;
          sum += x[i]; // Fastest ??
          ++n;
        } 
        sum = sum/n;
        return NumericVector::create((double)sum);
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
        sum = sum/l;
        return NumericVector::create((double)sum);
      }
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      if(narm) {
        NumericVector sum(ng, NA_REAL); // Other way ??
        IntegerVector n(ng, 1); // could also do no_init_vector and then add n[g[i]-1] = 1 in fir if condition... -> Nope, that is slower !!!
        for(int i = l; i--; ) { 
          if(!std::isnan(x[i])) { // faster way to code this ??? -> Not Bad at all -> index for g[i]-1?? -> Nope, no noticeable improvement !!
            if(std::isnan(sum[g[i]-1])) sum[g[i]-1] = x[i];
            else {
              sum[g[i]-1] += x[i]; 
              ++n[g[i]-1];
            }
          }
        }
        for(int i = ng; i--; ) sum[i] /= n[i];
        DUPLICATE_ATTRIB(sum, x);
        return sum;
      } else {
        NumericVector sum(ng); // no_init_vector // good?? -> yes, but not initializing is numerically unstable.. 
        int ngs = 0;
        if(Rf_isNull(gs)) {
          IntegerVector gsv(ng);
          if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
          for(int i = 0; i != l; ++i) { 
            if(std::isnan(x[i])) { 
              if(!std::isnan(sum[g[i]-1])) {
                sum[g[i]-1] = x[i];
                ++ngs;
                if(ngs == ng) break;
              }
            } else {
              sum[g[i]-1] += x[i];
              ++gsv[g[i]-1];
            }
          }
          for(int i = ng; i--; ) sum[i] /= gsv[i]; // Adding n takes twice as long, 
        } else {
          IntegerVector gsv = gs;
          if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
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
          for(int i = ng; i--; ) sum[i] /= gsv[i]; // This is good because adding n takes twice as long, if factor, supply gs = tabulate(f,nlevels(f))
        }
        DUPLICATE_ATTRIB(sum, x);
        return sum;
      }
    }
  } else { // With weights
    NumericVector wg = w; // wg(w) Identical speed
    if(l != wg.size()) stop("length(w) must match length(x)");
    if (ng == 0) {
      if(narm) {
        int j = l-1; // 1 because for-loop starts from 2
        while((std::isnan(x[j]) || std::isnan(wg[j])) && j!=0) --j; // This does not make a difference in performance but is more parsimonious. 
        long double sum = x[j]*wg[j], sumw = wg[j]; 
        if(j != 0) for(int i = j; i--; ) {
          if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
          sum += x[i]*wg[i]; // Fastest ??
          sumw += wg[i];
        } 
        sum = sum/sumw;
        return NumericVector::create((double)sum);
      } else {
        long double sum = 0, sumw = 0;
        for(int i = 0; i != l; ++i) { 
          if(std::isnan(x[i]) || std::isnan(wg[i])) { // good, check both ?? -> yes!!
            sum = x[i]+wg[i]; 
            break;
          } else { 
            sum += x[i]*wg[i];
            sumw += wg[i];
          }
        }
        sum = sum/sumw;
        return NumericVector::create((double)sum);
      }
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      if(narm) {
        NumericVector sum(ng, NA_REAL); // Other way ?? -> Nope, this is as good as it gets !!
        NumericVector sumw = no_init_vector(ng); // what if only NA ?? -> Works !! for some reason no problem !!, and faster !!
        for(int i = l; i--; ) { 
          if(std::isnan(x[i]) || std::isnan(wg[i])) continue; // faster way to code this ??? -> Not Bad at all -> index for g[i]-1?? -> Nope, no noticeable improvement !!
          if(std::isnan(sum[g[i]-1])) {
            sum[g[i]-1] = x[i]*wg[i];
            sumw[g[i]-1] = wg[i];
          } else {
            sum[g[i]-1] += x[i]*wg[i]; 
            sumw[g[i]-1] += wg[i];
          }
        }
        sum = sum/sumw; // good ?? better return sum/sumw?? -> Nope, slower !!
        DUPLICATE_ATTRIB(sum, x);
        return sum;
      } else {
        NumericVector sum(ng), sumw(ng); // good?? -> yes !! //  = no_init_vector// Not initializing numerically unstable !!
        int ngs = 0;
        for(int i = 0; i != l; ++i) { 
          if(std::isnan(x[i]) || std::isnan(wg[i])) { 
            if(!std::isnan(sum[g[i]-1])) {
              sum[g[i]-1] = sumw[g[i]-1] = x[i]+wg[i]; // or NA_REAL ?? -> Nope, good !!
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            sum[g[i]-1] += x[i]*wg[i];
            sumw[g[i]-1] += wg[i];
          }
        }
        sum = sum/sumw; 
        DUPLICATE_ATTRIB(sum, x);
        return sum;
      }
    }
  }
}



// // Final version without weights
// // [[Rcpp::export]]
// NumericVector fgmeanCpp(NumericVector x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, bool narm = true) {
//   int l = x.size();
//   if (ng == 0) {
//     if(narm) {
//       int j = l-1, n = 1; // 1 because for-loop starts from 2
//       double sum = x[j]; 
//       while(std::isnan(sum) && j!=0) sum = x[--j]; 
//       if(j != 0) for(int i = j; i--; ) {
//         if(std::isnan(x[i])) continue;
//         sum += x[i]; // Fastest ??
//         ++n;
//       } 
//       return NumericVector::create(sum/n);
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
//       return NumericVector::create(sum/l);
//     }
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     if(narm) {
//       NumericVector sum(ng, NA_REAL); // Other way ??
//       IntegerVector n(ng, 1); // could also do no_init_vector and then add n[g[i]-1] = 1 in fir if condition... -> Nope, that is slower !!!
//       for(int i = l; i--; ) { 
//         if(!std::isnan(x[i])) { // faster way to code this ??? -> Not Bad at all -> index for g[i]-1?? -> Nope, no noticeable improvement !!
//           if(std::isnan(sum[g[i]-1])) sum[g[i]-1] = x[i];
//           else {
//             sum[g[i]-1] += x[i]; 
//             ++n[g[i]-1];
//           }
//         }
//       }
//       for(int i = ng; i--; ) sum[i] /= n[i];
//       return sum;
//     } else {
//       NumericVector sum(ng); // good?? -> yes !!
//       int ngs = 0;
//       for(int i = 0; i != l; ++i) { 
//         if(std::isnan(x[i])) { 
//           if(!std::isnan(sum[g[i]-1])) {
//             sum[g[i]-1] = x[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         } else {
//           sum[g[i]-1] += x[i];
//         }
//       }
//       for(int i = ng; i--; ) sum[i] /= gs[i]; // This is good because adding n takes twice as long, if factor, supply gs = tabulate(f,nlevels(f))
//       return sum;
//     }
//   }
// }

// // Previous version -> non-efficient algorithms and internal return 
// // [[Rcpp::export]]
// NumericVector fgmeancpp(NumericVector x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, 
//                         bool narm = true, int ret = 0, bool fill = false) {
//   int l = x.size();
//   if (ng == 0) {
//     double sum = 0;
//     if(narm) {
//       int n = 0;
//         for(int i = l; i--; ) {
//           if(std::isnan(x[i])) continue; 
//           sum += x[i];
//           n++;
//         }
//         sum /= n;
//         if(ret == 1) { // possibility of reducing the number of passes??
//           if(fill) return rep(sum, l);
//           else {
//             NumericVector out(l, sum); 
//             for(int i = l; i--; ) if(std::isnan(x[i])) out[i] = x[i]; // 9.546032
//             return out;
//           }
//         } else if (ret == 2) {
//           if(fill) return x - sum; 
//           else {
//             NumericVector out(l);
//             for(int i = l; i--; ) { 
//               if(std::isnan(x[i])) out[i] = x[i]; // 9.266529
//               else out[i] = x[i] - sum; 
//             }
//             return out;
//           }
//         } else return NumericVector::create(sum);
//     } else {
//       for(int i = l; i--; ) sum += x[i];
//       sum /= l;
//       if(ret == 1) { // possibility of reducing the number of passes??
//         return rep(sum, l);
//       } else if (ret == 2) {
//         return x - sum; // x // x[i] -= sum[g[i]-1];
//       } else return NumericVector::create(sum);
//     }
//   } else { // with groups
//     NumericVector sum(gs);
//     if(narm) {
//       IntegerVector n(gs);
//         for(int i = l; i--; ) {
//           if(std::isnan(x[i])) continue; 
//           sum[g[i]-1] += x[i];
//           n[g[i]-1]++;
//         }
//         for(int i = ng; i--; ) sum[i] /= n[i];
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector out(l);
//           if(fill) {
//             for(int i = l; i--; ) out[i] = sum[g[i]-1];
//             return out;
//           } else {
//             for(int i = l; i--; ) {
//               if(std::isnan(x[i])) out[i] = x[i]; // 18.1844
//               else out[i] = sum[g[i]-1]; 
//             }
//             return out; 
//           }
//         } else if (ret == 2) {
//           NumericVector out(l);
//           if(fill) {
//             for(int i = l; i--; ) out[i] = x[i] - sum[g[i]-1];
//             return out;
//           } else {
//             for(int i = l; i--; ) {
//               if(std::isnan(x[i])) out[i] = x[i]; 
//               else out[i] = x[i] - sum[g[i]-1]; 
//             }
//             return out; 
//           }
//         } else return sum; 
//     } else {
//       if(gs.size() == 1) {
//         IntegerVector n(gs);
//         for(int i = l; i--; ) {
//           sum[g[i]-1] += x[i];
//           n[g[i]-1]++;
//         }
//         for(int i = ng; i--; ) sum[i] /= n[i];
//       } else {
//         for(int i = l; i--; ) sum[g[i]-1] += x[i];
//         for(int i = ng; i--; ) sum[i] /= gs[i];
//       }
//       if(ret == 1) { // possibility of reducing the number of passes??
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




// Version before reducing code::
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// NumericVector fgmeancpp(NumericVector x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, 
//                         bool narm = true, int ret = 0, bool fill = false) {
//   int l = x.size();
//   if (ng == 0) {
//     double sum = 0;
//     if(narm) {
//       int n = 0;
//       if(fill) {
//         for(int i = l; i--; ) {
//           if(std::isnan(x[i])) continue; // what about NA's with ret = 1 or 2
//           sum += x[i];
//           n++;
//         }
//         sum /= n;
//         if(ret == 1) { // possibility of reducing the number of passes??
//           return rep(sum, l);
//         } else if (ret == 2) {
//           return x - sum; // x // x[i] -= sum[g[i]-1];
//         } else return NumericVector::create(sum);
//       } else {
//        // LogicalVector isnan(l); // This is taking time -> Allocating space !!!
//         for(int i = l; i--; ) {
//           // isnan[i] = std::isnan(x[i]); // 4.944284
//           // if(isnan[i]) continue; // what about NA's with ret = 1 or 2
//           if(std::isnan(x[i])) continue; // 2.368566
//           sum += x[i];
//           n++;
//         }
//         sum /= n;
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector out(l, sum); 
//           // for(int i = l; i--; ) if(isnan[i]) out[i] = x[i]; // 11.97013  // faster way?? vector subsetting??
//           for(int i = l; i--; ) if(std::isnan(x[i])) out[i] = x[i]; // 9.546032
//           return out;
//         } else if (ret == 2) {
//           NumericVector out(l);
//           for(int i = l; i--; ) { 
//             // if(isnan[i]) out[i] = x[i]; // 11.85642  // faster way?? vector subsetting??
//             if(std::isnan(x[i])) out[i] = x[i]; // 9.266529
//             else out[i] = x[i] - sum; 
//           }
//           return out;
//         } else return NumericVector::create(sum);
//       }
//     } else {
//       for(int i = l; i--; ) sum += x[i];
//       sum /= l;
//       if(ret == 1) { // possibility of reducing the number of passes??
//         return rep(sum, l);
//       } else if (ret == 2) {
//         return x - sum; // x // x[i] -= sum[g[i]-1];
//       } else return NumericVector::create(sum);
//     }
//   } else { // with groups
//     NumericVector sum(gs);
//     if(narm) {
//       IntegerVector n(gs);
//       if(fill) {
//         for(int i = l; i--; ) {
//           if(std::isnan(x[i])) continue; 
//           sum[g[i]-1] += x[i];
//           n[g[i]-1]++;
//         }
//         for(int i = ng; i--; ) sum[i] /= n[i];
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector out(l);
//           for(int i = l; i--; ) out[i] = sum[g[i]-1];
//           return out; // return sum[g-1];
//         } else if (ret == 2) {
//           NumericVector out(l);
//           for(int i = l; i--; )  out[i] = x[i] - sum[g[i]-1]; // x[i] -= sum[g[i]-1];
//           return out; // x
//         } else return sum; 
//       } else {
//         // LogicalVector isnan(l);
//         for(int i = l; i--; ) {
//           // isnan[i] = std::isnan(x[i]); // 11.79592
//           // if(isnan[i]) continue;
//           if(std::isnan(x[i])) continue; // 8.242974
//           sum[g[i]-1] += x[i];
//           n[g[i]-1]++;
//         }
//         for(int i = ng; i--; ) sum[i] /= n[i];
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector out(l);
//           for(int i = l; i--; ) {
//             // if(isnan[i]) out[i] = x[i]; // 21.71893 // faster way?? vector subsetting??
//             if(std::isnan(x[i])) out[i] = x[i]; // 18.1844
//             else out[i] = sum[g[i]-1]; 
//           }
//           return out; // return sum[g-1];
//         } else if (ret == 2) {
//           NumericVector out(l);
//           for(int i = l; i--; ) {
//             // if(isnan[i]) out[i] = x[i]; // 21.61112  // faster way?? vector subsetting??
//             if(std::isnan(x[i])) out[i] = x[i]; // 18.1551
//             else out[i] = x[i] - sum[g[i]-1]; 
//           }
//           return out; // return sum[g-1];
//         } else return sum; 
//       }
//     } else {
//      if(gs.size() == 1) {
//       IntegerVector n(gs);
//       for(int i = l; i--; ) {
//         sum[g[i]-1] += x[i];
//         n[g[i]-1]++;
//       }
//       for(int i = ng; i--; ) sum[i] /= n[i];
//      } else {
//        for(int i = l; i--; ) sum[g[i]-1] += x[i];
//        for(int i = ng; i--; ) sum[i] /= gs[i];
//      }
//      if(ret == 1) { // possibility of reducing the number of passes??
//        NumericVector out(l);
//        for(int i = l; i--; ) out[i] = sum[g[i]-1];
//        return out; // return sum[g-1]; this is slower!!
//      } else if (ret == 2) {
//        NumericVector out(l);
//        for(int i = l; i--; ) out[i] = x[i] - sum[g[i]-1]; // x[i] -= sum[g[i]-1];
//        return out; // x
//      } else return sum; 
//    }
//   }
// }

// // [[Rcpp::plugins(cpp11)]]
// #include <numeric>
// #include <Rcpp.h>
// using namespace Rcpp;

// Old version before introducing some efficient vectorizations:
// // [[Rcpp::export]]
// NumericVector fgmeancpp(NumericVector x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, bool narm = true, int ret = 0, bool fill = false) {
//   int l = x.size();
//   if (ng == 0) {
//     double sum = 0;
//     if(narm) {
//       int n = 0;
//       if(fill) {
//         for(int i = l; i--; ) {
//           if(ISNAN(x[i])) continue; // what about NA's with ret = 1 or 2
//           sum += x[i];
//           n++;
//         }
//         sum /= n;
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector out(l);
//           for(int i = l; i--; ) out[i] = sum;
//           return out;
//         } else if (ret == 2) {
//           NumericVector out(l);
//           for(int i = l; i--; )  out[i] = x[i] - sum; // x[i] -= sum[g[i]-1];
//           return out; // x
//         } else return sum;
//       } else {
//         LogicalVector isnan(l);
//         for(int i = l; i--; ) {
//           isnan[i] = ISNAN(x[i]);
//           if(isnan[i]) continue; // what about NA's with ret = 1 or 2
//           sum += x[i];
//           n++;
//         }
//         sum /= n;
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector out(l);
//           for(int i = l; i--; ) {
//             if(isnan[i]) out[i] = x[i];  // faster way?? vector subsetting??
//             else out[i] = sum; 
//           }
//           return out;
//         } else if (ret == 2) {
//           NumericVector out(l);
//           for(int i = l; i--; ) {
//             if(isnan[i]) out[i] = x[i];  // faster way?? vector subsetting??
//             else out[i] = x[i] - sum; 
//           }
//           return out;
//         } else return sum;
//       }
//     } else {
//       for(int i = l; i--; ) sum += x[i];
//       sum /= l;
//       if(ret == 1) { // possibility of reducing the number of passes??
//         NumericVector out(l);
//         for(int i = l; i--; ) out[i] = sum;
//         return out;
//       } else if (ret == 2) {
//         NumericVector out(l);
//         for(int i = l; i--; )  out[i] = x[i] - sum; // x[i] -= sum[g[i]-1];
//         return out; // x
//       } else return sum;
//     }
//   } else { // with groups
//     NumericVector sum(gs);
//     if(narm) {
//       IntegerVector n(gs);
//       if(fill) {
//         for(int i = l; i--; ) {
//           if(ISNAN(x[i])) continue; 
//           sum[g[i]-1] += x[i];
//           n[g[i]-1]++;
//         }
//         for(int i = ng; i--; ) sum[i] /= n[i];
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector out(l);
//           for(int i = l; i--; ) out[i] = sum[g[i]-1];
//           return out; // return sum[g-1]; // This is slower !!
//         } else if (ret == 2) {
//           NumericVector out(l);
//           for(int i = l; i--; )  out[i] = x[i] - sum[g[i]-1]; // x[i] -= sum[g[i]-1];
//           return out; // x
//         } else return sum; 
//       } else {
//         LogicalVector isnan(l);
//         for(int i = l; i--; ) {
//           isnan[i] = ISNAN(x[i]);
//           if(isnan[i]) continue; 
//           sum[g[i]-1] += x[i];
//           n[g[i]-1]++;
//         }
//         for(int i = ng; i--; ) sum[i] /= n[i];
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector out(l);
//           for(int i = l; i--; ) {
//             if(isnan[i]) out[i] = x[i];  // faster way?? vector subsetting??
//             else out[i] = sum[g[i]-1]; 
//           }
//           return out; // return sum[g-1];
//         } else if (ret == 2) {
//           NumericVector out(l);
//           for(int i = l; i--; ) {
//             if(isnan[i]) out[i] = x[i];  // faster way?? vector subsetting??
//             else out[i] = x[i] - sum[g[i]-1]; 
//           }
//           return out; // return sum[g-1];
//         } else return sum; 
//       }
//     } else {
//       if(gs.size() == 1) {
//         IntegerVector n(gs);
//         for(int i = l; i--; ) {
//           sum[g[i]-1] += x[i];
//           n[g[i]-1]++;
//         }
//         for(int i = ng; i--; ) sum[i] /= n[i];
//       } else {
//         for(int i = l; i--; ) sum[g[i]-1] += x[i];
//         for(int i = ng; i--; ) sum[i] /= gs[i];
//       }
//       if(ret == 1) { // possibility of reducing the number of passes??
//         NumericVector out(l);
//         for(int i = l; i--; ) out[i] = sum[g[i]-1];
//         return out; // return sum[g-1];
//       } else if (ret == 2) {
//         NumericVector out(l);
//         for(int i = l; i--; )  out[i] = x[i] - sum[g[i]-1]; // x[i] -= sum[g[i]-1];
//         return out; // x
//       } else return sum; 
//     }
//   }
// }