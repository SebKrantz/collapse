#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fsumlCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
              bool narm = true, bool drop = true) {
  int l = x.size();
  
  if (ng == 0) {
    NumericVector sum(l); // not initializing not faster WIth NWDI (35 instead of 32 milliseconds)
    if(narm) {
      for(int j = l; j--; ) {
      // for(int j = 0; j != l; ++j) { // Not necessarily faster !!
        NumericVector column = x[j];
        int k = column.size()-1;
        long double sumi = column[k]; // a bit extra speed with double, 31 vs 36 milliseconds on NWDI
        while(std::isnan(sumi) && k!=0) sumi = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(!std::isnan(column[i])) sumi += column[i];
        }
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
        for(int i = gss; i--; ) {
          if(std::isnan(column[i])) continue; // faster way to code this ??? -> Not Bad at all, 54.. millisec for WDIM
            if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
            else sumj[g[i]-1] += column[i];
        }
        SHALLOW_DUPLICATE_ATTRIB(sumj, column);
        sum[j] = sumj;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector sumj(ng); //  = no_init_vector(ng); // Not initializing in loop is numerically unstable !!
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
        SHALLOW_DUPLICATE_ATTRIB(sumj, column);
        sum[j] = sumj;
      }
    }
    DUPLICATE_ATTRIB(sum, x);
    sum.attr("row.names") = NumericVector::create(NA_REAL, -ng);
    return sum;
  }
}


// // Optimal Version witj built in return!!! Formerly fgsumlCppSW2 -> Not all in one loop, a tiny bit slower but more practical, with 3 return and a drop option !!!
// // [[Rcpp::export]]
// SEXP fgsumlCpp(List x, int ng = 0, IntegerVector g = 0,
//                   bool narm = true, int ret = 0, bool drop = true) {
//   int l = x.size(), row = 0;
//   
//   if (ng == 0) {
//     NumericVector sum(l); // not initializing not faster WIth NWDI (35 instead of 32 milliseconds)
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         int k = row-1;
//         double sumi = column[k];
//         while(std::isnan(sumi) && k!=0) sumi = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(!std::isnan(column[i])) sumi += column[i];
//         }
//         sum[j] = sumi;
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
//         sum[j] = sumi;
//       }
//     }
//   switch(ret) {
//     case 1: 
//     {
//       List out(l);
//       for(int i = l; i--; ) out[i] = rep(sum[i], row);
//       return out;
//     }
//     case 2: 
//     {
//       List out(l);
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         NumericVector sgj(row); //  = no_init_vector withoit initialization, above is faster !!. This is definitely faster than clone !!
//         double sumj = sum[j];
//         for(int i = row; i--; ) {
//           if(std::isnan(column[i])) sgj[i] = column[i];
//           else sgj[i] = sumj;
//         }
//         out[j] = sgj; 
//       }
//       return out;
//     }
//     case 3: 
//     {
//       List out(l);
//       for(int j = l; j--; ) { // Could make faster ?? 
//         NumericVector column = x[j];
//         out[j] = column - sum[j];
//       }
//       return out;
//     }
//     default: 
//     {
//       if(drop) {
//         sum.attr("names") = x.attr("names");
//         return sum;
//       } else {
//         List out(l);
//         for(int j = l; j--; ) out[j] = sum[j];
//         return out;
//       }
//     }
//    }
//   } else { // With groups !!
//     List sum(l);
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng, NA_REAL);
//         for(int i = row; i--; ) {
//           if(!std::isnan(column[i])) { // faster way to code this ??? -> Not Bad at all, 54.. millisec for WDIM
//             if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//             else sumj[g[i]-1] += column[i];
//           }
//         }
//         sum[j] = sumj;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj = no_init_vector(ng);
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
//         sum[j] = sumj;
//       }
//     }
//   switch(ret) {
//     case 1: 
//     {
//       for(int j = l; j--; ) {
//       NumericVector sgj = no_init_vector(row);
//       NumericVector sumj = sum[j];
//       for(int i = row; i--; ) sgj[i] = sumj[g[i]-1];
//       sum[j] = sgj; // clone(sgj); // overwriting fastest ?? 
//     }
//       return sum;
//     }
//     case 2: 
//     {
//       for(int j = l; j--; ) {
//       NumericVector sgj = no_init_vector(row);
//       NumericVector column = x[j];
//       NumericVector sumj = sum[j];
//       for(int i = row; i--; ) {
//         if(std::isnan(column[i])) sgj[i] = column[i];
//         else sgj[i] = sumj[g[i]-1];
//       }
//       sum[j] = sgj;
//     }
//       return sum;
//     }
//     case 3: 
//     {
//       for(int j = l; j--; ) {
//       NumericVector sgj = no_init_vector(row);
//       NumericVector column = x[j];
//       NumericVector sumj = sum[j];
//       for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1];
//       sum[j] = sgj; // clone(sgj); Seems slower.. 
//       }
//       return sum;
//     }
//     default: 
//       return sum;
//     }
//   }
// }

// Putting 3 options but with a switch -> A bit faster than previous (second fastest)
// // [[Rcpp::export]]
// List fgsumlCppSW(List x, int ng = 0, IntegerVector g = 0,
//                  bool narm = true, int ret = 0, bool drop = true) {
//   int l = x.size();
//   List sum(l);
//   
//   if (ng == 0) {
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         int k = row-1;
//         double sumi = column[k];
//         while(std::isnan(sumi) && k!=0) sumi = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(!std::isnan(column[i])) sumi += column[i];
//         }
//         switch(ret) {
//         case 1: 
//           sum[j] = rep(sumi, row);
//           break;
//         case 2: 
//           {
//             NumericVector sgj(row, sumi); // fastest??
//             for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i];
//             sum[j] = sgj;
//             break;
//           }
//         case 3: 
//           sum[j] = column - sumi;
//           break;
//         default: 
//           sum[j] = sumi;
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double sumi = 0;
//         int row = column.size();
//         for(int i = 0; i != row; ++i) {
//           if(std::isnan(column[i])) {
//             sumi = column[i];
//             break;
//           } else {
//             sumi += column[i];
//           }
//         }
//         switch(ret) {
//         case 1: 
//           sum[j] = rep(sumi, row);
//           break;
//         case 2: 
//           {
//             NumericVector sgj(row, sumi);
//             for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i];
//             sum[j] = sgj;
//             break;
//           }
//         case 3: 
//           sum[j] = column - sumi;
//           break;
//         default: 
//           sum[j] = sumi;
//         }
//       }
//     }
//     // if(drop) { // Nah, not so brilliant
//     //   NumericVector out(l);
//     //   for(int i = l; i--; ) out[i] = sum[i];
//     //   return out;
//     // }
//   } else { // With groups !!
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng, NA_REAL);
//         for(int i = row; i--; ) {
//           if(!std::isnan(column[i])) { // faster way to code this ??? -> Not Bad at all, 54.. millisec for WDIM
//             if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//             else sumj[g[i]-1] += column[i];
//           }
//         }
//         switch(ret) {
//         case 1: 
//         {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = sumj[g[i]-1];
//           sum[j] = sgj;
//           break;
//         }
//         case 2: 
//         {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i];
//             else sgj[i] = sumj[g[i]-1];
//           }
//           sum[j] = sgj;
//           break;
//         }
//         case 3: 
//         {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1];
//           sum[j] = sgj;
//           break;
//         }
//         default: 
//           sum[j] = sumj;
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng);
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
//         switch(ret) {
//         case 1: 
//         {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = sumj[g[i]-1];
//           sum[j] = sgj;
//           break;
//         }
//         case 2: 
//         {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i];
//             else sgj[i] = sumj[g[i]-1];
//           }
//           sum[j] = sgj;
//           break;
//         }
//         case 3: 
//         {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1];
//           sum[j] = sgj;
//           break;
//         }
//         default: 
//           sum[j] = sumj;
//         }
//       }
//     }
//   }
//   return sum;
// }


// Giving 3 return options -> Not quite as fast as prvious !!
// // [[Rcpp::export]]
// List fgsumlCppv2(List x, int ng = 0, IntegerVector g = 0, 
//                  bool narm = true, int ret = 0) { 
//   int l = x.size(); 
//   List sum(l); 
//   
//   if (ng == 0) { 
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size(); 
//         int k = row-1;
//         double sumi = column[k];
//         while(std::isnan(sumi) && k!=0) sumi = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(!std::isnan(column[i])) sumi += column[i];
//         }
//         if(ret == 1) {
//           sum[j] = rep(sumi, row);
//         } else if(ret == 2) {
//           NumericVector sgj(row, sumi);
//           for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i];
//           sum[j] = sgj;
//         } else if (ret == 3) {
//           sum[j] = column - sumi; 
//         } else {
//           sum[j] = sumi;
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         // sumi = column[0]; This would give a mistake !!!
//         double sumi = 0; // This is necessary !!
//         int row = column.size();
//         for(int i = 0; i != row; ++i) { 
//           if(std::isnan(column[i])) {
//             sumi = column[i]; 
//             break;
//           } else { 
//             sumi += column[i];
//           }
//         }
//         if(ret == 1) {
//           sum[j] = rep(sumi, row);
//         } else if(ret == 2) {
//           NumericVector sgj(row, sumi);
//           for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i];
//           sum[j] = sgj;
//         } else if (ret == 3) {
//           sum[j] = column - sumi; 
//         } else {
//           sum[j] = sumi;
//         }
//       }
//     }
//   } else { // With groups !!
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng, NA_REAL);
//         for(int i = row; i--; ) {
//           if(!std::isnan(column[i])) { // faster way to code this ??? -> Not Bad at all, 54.. millisec for WDIM
//             if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//             else sumj[g[i]-1] += column[i];
//           }
//         }
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; 
//           sum[j] = sgj;
//         } else if(ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i]; 
//             else sgj[i] = sumj[g[i]-1];
//           }
//           sum[j] = sgj;
//         } else if (ret == 3) {
//           NumericVector sgj(row); 
//           for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1]; 
//           sum[j] = sgj;
//         } else {
//           sum[j] = sumj;
//         }
//       }
//     } else { 
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng);
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
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; 
//           sum[j] = sgj;
//         } else if(ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i]; 
//             else sgj[i] = sumj[g[i]-1];
//           }
//           sum[j] = sgj;
//         } else if (ret == 3) {
//           NumericVector sgj(row); 
//           for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1]; 
//           sum[j] = sgj;
//         } else {
//           sum[j] = sumj;
//         }
//       }
//     }
//   }
//   return sum;  
// }


// Fastest Version !!! Adjusted Algorithms! But above is just as fast and more practical..
// // [[Rcpp::export]]
// List fgsumlCpp(List x, int ng = 0, IntegerVector g = 0, 
//                bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.size(); 
//   List sum(l); 
//   
//   if (ng == 0) { 
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size(); 
//         int k = row-1;
//         double sumi = column[k];
//         while(std::isnan(sumi) && k!=0) sumi = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(!std::isnan(column[i])) sumi += column[i];
//         }
//         if(ret == 1) {
//           if(fill) sum[j] = rep(sumi, row);
//           else {
//             NumericVector sgj(row, sumi);
//             for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i];
//             sum[j] = sgj;
//           }
//         } else if (ret == 2) {
//           if(fill) sum[j] = column - sumi; 
//           else {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i];
//               else sgj[i] = column[i] - sumi; 
//             }
//             sum[j] = sgj;
//           }
//         } else sum[j] = sumi;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         // sumi = column[0]; This would give a mistake !!!
//         double sumi = 0; // This is necessary !!
//         int row = column.size();
//         for(int i = 0; i != row; ++i) { 
//           if(std::isnan(column[i])) {
//             sumi = column[i]; 
//             break;
//           } else { 
//             sumi += column[i];
//           }
//         }
//         if(ret == 1) { 
//           sum[j] = rep(sumi, row);
//         } else if (ret == 2) {
//           sum[j] = column - sumi;
//         } else sum[j] = sumi;
//       }
//     }
//   } else { // With groups !!
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng, NA_REAL);
//         for(int i = row; i--; ) {
//           if(!std::isnan(column[i])) { // faster way to code this ??? -> Not Bad at all, 54.. millisec for WDIM
//             if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//             else sumj[g[i]-1] += column[i];
//           }
//         }
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           if(fill) { 
//             for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; 
//           } else {
//             for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i]; 
//               else sgj[i] = sumj[g[i]-1];
//             }
//           }
//           sum[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row); 
//           if(fill) { 
//             for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1]; 
//           } else {
//             for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i]; 
//               else sgj[i] = column[i] - sumj[g[i]-1];
//             }
//           }
//           sum[j] = sgj;
//         } else sum[j] = sumj;
//       }
//     } else { 
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng);
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
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; 
//           sum[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1];
//           sum[j] = sgj;
//         } else sum[j] = sumj;
//       }
//     }
//   }
//   return sum;  
// }


// Minimal Example !!!
// // [[Rcpp::export]]
// List fgsumlCppmin(List x, int ng = 0, IntegerVector g = 0, bool narm = true) {
//   int l = x.size();
//   List sum(l);
//   
//   if (ng == 0) {
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         int k = row-1;
//         double sumi = column[k];
//         while(std::isnan(sumi) && k!=0) sumi = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(!std::isnan(column[i])) sumi += column[i];
//         }
//         sum[j] = sumi;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double sumi = 0; // This is necessary !!
//         int row = column.size();
//         for(int i = 0; i != row; ++i) {
//           if(std::isnan(column[i])) {
//             sumi = column[i];
//             break;
//           } else {
//             sumi += column[i];
//           }
//         }
//         sum[j] = sumi;
//       }
//     }
//   } else { // With groups !!
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng, NA_REAL);
//         for(int i = row; i--; ) {
//           if(!std::isnan(column[i])) { // faster way to code this ??? -> Not Bad at all, 54.. millisec for WDIM
//             if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//             else sumj[g[i]-1] += column[i];
//           }
//         }
//         sum[j] = sumj;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng);
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
//         sum[j] = sumj;
//       }
//     }
//   }
//   return sum;
// }


// Previous Version (before adjusting algorithms)
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// List fgsumlCpp(List x, int ng = 0, IntegerVector g = 0, 
//                bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.size(); 
//   List sum(l); 
//   
//   if (ng == 0) { 
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double sumi = 0;
//         int row = column.size();
//         for(int i = row; i--; ) {
//           if(std::isnan(column[i])) continue; 
//           sumi += column[i];
//         }
//         if(ret == 1) {
//           if(fill) sum[j] = rep(sumi, row);
//           else {
//             NumericVector sgj(row, sumi);
//             for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i];
//             sum[j] = sgj;
//           }
//         } else if (ret == 2) {
//           if(fill) sum[j] = column - sumi; 
//           else {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i];
//               else sgj[i] = column[i] - sumi; 
//             }
//             sum[j] = sgj;
//           }
//         } else sum[j] = sumi;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double sumi = 0;
//         int row = column.size();
//         for(int i = row; i--; ) sumi += column[i];
//         if(ret == 1) { 
//           sum[j] = rep(sumi, row);
//         } else if (ret == 2) {
//           sum[j] = column - sumi;
//         } else sum[j] = sumi;
//       }
//     }
//   } else { // With groups !!
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng);
//         for(int i = row; i--; ) {
//           if(std::isnan(column[i])) continue; 
//           sumj[g[i]-1] += column[i];  
//         }
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           if(fill) for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; 
//           else for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i]; 
//             else sgj[i] = sumj[g[i]-1];
//           }
//           sum[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row); 
//           if(fill) for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1]; 
//           else for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i]; 
//             else sgj[i] = column[i] - sumj[g[i]-1];
//           }
//           sum[j] = sgj;
//         } else sum[j] = sumj;
//       }
//     } else { 
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector sumj(ng);
//         for(int i = row; i--; ) sumj[g[i]-1] += column[i]; 
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = sumj[g[i]-1]; 
//           sum[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1];
//           sum[j] = sgj;
//         } else sum[j] = sumj;
//       }
//     }
//   }
//   return sum;  
// }