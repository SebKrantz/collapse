#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fminlCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
              bool narm = true, bool drop = true) {
  int l = x.size();
  
  if (ng == 0) {
    NumericVector min = no_init_vector(l); // good and fast here !!
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        int k = column.size()-1;
        double mini = column[k];
        while(std::isnan(mini) && k!=0) mini = column[--k];
        if(k != 0) for(int i = k; i--; ) {
            if(mini > column[i]) mini = column[i];
        }
        min[j] = mini;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        double mini = column[0];
        int row = column.size();
        for(int i = 0; i != row; ++i) {
          if(std::isnan(column[i])) {
            mini = column[i];
            break;
          } else {
            if(mini > column[i]) mini = column[i];
          }
        }
        min[j] = mini;
      }
    }
    if(drop) {
      min.attr("names") = x.attr("names");
      return min;
    } else {
      List out(l);
      for(int j = l; j--; ) {
        out[j] = min[j];
        SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
      }
      DUPLICATE_ATTRIB(out, x);
      out.attr("row.names") = 1;
      return out;
    }
  } else { // With groups !!
    List min(l);
    int gss = g.size();
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector minj(ng, NA_REAL);
        for(int i = gss; i--; ) {
          if(!std::isnan(column[i])) { // Keeping this is faster !!!!
            if(minj[g[i]-1] > column[i] || std::isnan(minj[g[i]-1])) minj[g[i]-1] = column[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(minj, column);
        min[j] = minj;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector minj(ng, R_PosInf); 
        int ngs = 0;
        for(int i = 0; i != gss; ++i) {
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
        SHALLOW_DUPLICATE_ATTRIB(minj, column);
        min[j] = minj;
      }
    }
    DUPLICATE_ATTRIB(min, x);
    min.attr("row.names") = NumericVector::create(NA_REAL, -ng);
    return min;
  }
}

// Previous Version: With return, but fully efficient !!
// // [[Rcpp::export]]
// List fgminlCpp(List x, int ng = 0, IntegerVector g = 0, 
//                 bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.size(); 
//   List min(l); 
//   
//   if (ng == 0) { 
//     if(narm) {
//       int row = 0; // extra speed ??, not defining again ?? 
//       int k = 0;
//       double mini = 0;
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size(); // necessary here, could put only if ret != 0
//         k = row-1;
//         mini = column[k];
//         while(std::isnan(mini) && k!=0) mini = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(mini > column[i]) mini = column[i];
//         }
//         if(ret == 1) {
//           if(fill) min[j] = rep(mini, row);
//           else {
//             NumericVector sgj(row, mini);
//             for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i];
//             min[j] = sgj;
//           }
//         } else if (ret == 2) {
//           if(fill) min[j] = column - mini; 
//           else {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i];
//               else sgj[i] = column[i] - mini; 
//             }
//             min[j] = sgj;
//           }
//         } else min[j] = mini;
//       }
//     } else {
//       double mini = 0;
//       int row = 0;
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         mini = column[0];
//         row = column.size();
//         for(int i = 0; i != row; ++i) { 
//           if(std::isnan(column[i])) {
//             mini = column[i]; 
//             break;
//           } else { 
//             if(mini > column[i]) mini = column[i];
//           }
//         }
//         if(ret == 1) { 
//           min[j] = rep(mini, row);
//         } else if (ret == 2) {
//           min[j] = column - mini;
//         } else min[j] = mini;
//       }
//     }
//   } else { // With groups !!
//     int gss = g.size();
//     if(narm) {
//       int row = 0;
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector minj(ng, NA_REAL);
//         for(int i = row; i--; ) {
//           if(minj[g[i]-1] > column[i] || std::isnan(minj[g[i]-1])) minj[g[i]-1] = column[i]; 
//         }
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           if(fill) for(int i = row; i--; ) sgj[i] = minj[g[i]-1]; 
//           else for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i]; 
//             else sgj[i] = minj[g[i]-1];
//           }
//           min[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row); 
//           if(fill) for(int i = row; i--; ) sgj[i] = column[i] - minj[g[i]-1]; 
//           else for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i]; 
//             else sgj[i] = column[i] - minj[g[i]-1];
//           }
//           min[j] = sgj;
//         } else min[j] = minj;
//       }
//     } else { 
//       int row = 0;
//       int ngs = 0;
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector minj(ng, INFINITY);
//         ngs = 0;
//         for(int i = 0; i != row; ++i) {
//           if(std::isnan(column[i])) {
//             if(!std::isnan(minj[g[i]-1])) {
//               minj[g[i]-1] = column[i]; 
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           } else { 
//             if(minj[g[i]-1] > column[i]) minj[g[i]-1] = column[i];
//           }
//         }
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = minj[g[i]-1]; 
//           min[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = column[i] - minj[g[i]-1];
//           min[j] = sgj;
//         } else min[j] = minj;
//       }
//     }
//   }
//   return min;  
// }


// // [[Rcpp::plugins(cpp11)]]
// #include <numeric>
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // Note: it seems that changing the g[i]-1 doesn't make things faster!!!
// 
// // gs IS OF NO USE HERE!!
// 
// // [[Rcpp::export]]
// List fgminlcpp(List x, IntegerVector g = 0, int ng = 0, IntegerVector gs = 0, bool narm = true, int ret = 0, bool fill = false) { // void
//   int l = x.size();
//   List min(l);
// 
//   if (ng == 0) { // No groups
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double mini = 0;
//           int row = column.size();
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             mini *= column[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = mini; // faster with loop??: for(int i = row; i--; ) sgj[i] = minj[g[i]-1];
//             min[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - mini;
//             min[j] = sgj;
//           } else min[j] = mini;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double mini = 0;
//           int row = column.size();
//           LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             isnan[i] = ISNAN(column[i]);
//             if(isnan[i]) continue; // faster way??
//             mini *= column[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = mini; 
//             }
//             min[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = column[i] - mini; 
//             }
//             min[j] = sgj;
//           } else min[j] = mini;
//         }        
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double mini = 0;
//         int row = column.size();
//         for(int i = row; i--; ) mini *= column[i];
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = mini; // faster with loop??: for(int i = row; i--; ) sgj[i] = minj[g[i]-1];
//           min[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = column[i] - mini;
//           min[j] = sgj;
//         } else min[j] = mini;
//       }
//     }
//   } else { // With groups
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector minj(ng);
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             minj[g[i]-1] *= column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = minj[g[i]-1]; // Fastest way?? better loop through cols??
//             min[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - minj[g[i]-1];
//             min[j] = sgj;
//           } else min[j] = minj;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector minj(ng);
//           LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             isnan[i] = ISNAN(column[i]);
//             if(isnan[i]) continue; // faster way??
//             minj[g[i]-1] *= column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = minj[g[i]-1]; 
//             }
//             min[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = column[i] - minj[g[i]-1]; 
//             }
//             min[j] = sgj;
//           } else min[j] = minj;
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         NumericVector minj(ng);
//         for(int i = row; i--; ) minj[g[i]-1] *= column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = minj[g[i]-1]; // Fastest way?? better loop through cols??
//           min[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = column[i] - minj[g[i]-1];
//           min[j] = sgj;
//         } else min[j] = minj;
//       }
//     }
//   }
//   return min;
// }



// // [[Rcpp::plugins(cpp11)]]
// #include <numeric>
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // Note: it seems that changing the g[i]-1 doesn't make things faster!!!
// 
// // gs IS OF NO USE HERE!!
// 
// // [[Rcpp::export]]
// List fgminlcpp(List x, IntegerVector g = 0, int ng = 0, IntegerVector gs = 0, bool narm = true, int ret = 0, bool fill = false) { // void
//   int l = x.size();
//   List min(l);
//   
//   if (ng == 0) { // No groups
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double mini = R_PosInf;
//           int row = column.size();
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             if(mini>column[i]) mini = column[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = mini; // faster with loop??: for(int i = row; i--; ) sgj[i] = minj[g[i]-1];
//             min[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - mini;
//             min[j] = sgj;
//           } else min[j] = mini;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double mini = R_PosInf;
//           int row = column.size();
//           LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             isnan[i] = ISNAN(column[i]);
//             if(isnan[i]) continue; // faster way??
//             if(mini>column[i]) mini = column[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = mini; 
//             }
//             min[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = column[i] - mini; 
//             }
//             min[j] = sgj;
//           } else min[j] = mini;
//         }        
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double mini = R_PosInf;
//         int row = column.size();
//         for(int i = row; i--; ) if(mini>column[i]) mini = column[i];
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = mini; // faster with loop??: for(int i = row; i--; ) sgj[i] = minj[g[i]-1];
//           min[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = column[i] - mini;
//           min[j] = sgj;
//         } else min[j] = mini;
//       }
//     }
//   } else { // With groups
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector minj(ng, R_PosInf); // Best way?? 
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             if(minj[g[i]-1]>column[i]) minj[g[i]-1] = column[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = minj[g[i]-1]; // Fastest way?? better loop through cols??
//             min[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - minj[g[i]-1];
//             min[j] = sgj;
//           } else min[j] = minj;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector minj(ng, R_PosInf); // Best way??
//           LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             isnan[i] = ISNAN(column[i]);
//             if(isnan[i]) continue; // faster way??
//             if(minj[g[i]-1]>column[i]) minj[g[i]-1] = column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = minj[g[i]-1]; 
//             }
//             min[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = column[i] - minj[g[i]-1]; 
//             }
//             min[j] = sgj;
//           } else min[j] = minj;
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         NumericVector minj(ng, R_PosInf); // Best way??
//         for(int i = row; i--; ) if(minj[g[i]-1]>column[i]) minj[g[i]-1] = column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = minj[g[i]-1]; // Fastest way?? better loop through cols??
//           min[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = column[i] - minj[g[i]-1];
//           min[j] = sgj;
//         } else min[j] = minj;
//       }
//     }
//   }
//   return min;
// }