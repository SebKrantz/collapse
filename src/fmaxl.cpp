#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fmaxlCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
              bool narm = true, bool drop = true) {
  int l = x.size();
  
  if (ng == 0) {
    NumericVector max = no_init_vector(l); // Good and faster !!
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        int k = column.size()-1;
        double maxi = column[k];
        while(std::isnan(maxi) && k!=0) maxi = column[--k];
        if(k != 0) for(int i = k; i--; ) { // This is fastest, checking isnan(column[i]) does not give extra speed!!
          if(maxi < column[i]) maxi = column[i];
        }
        max[j] = maxi;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        double maxi = column[0];
        int row = column.size();
        for(int i = 0; i != row; ++i) {
          if(std::isnan(column[i])) {
            maxi = column[i];
            break;
          } else {
            if(maxi < column[i]) maxi = column[i];
          }
        }
        max[j] = maxi;
      }
    }
    if(drop) {
      max.attr("names") = x.attr("names");
      return max;
    } else {
      List out(l);
      for(int j = l; j--; ) {
        out[j] = max[j];
        SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
      }
      DUPLICATE_ATTRIB(out, x);
      out.attr("row.names") = 1;
      return out;
    }
  } else { // With groups !!
    List max(l);
    int gss = g.size();
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector maxj(ng, NA_REAL);
        for(int i = gss; i--; ) {
          if(!std::isnan(column[i])) { // Keeping this is faster !!!!
            if(maxj[g[i]-1] < column[i] || std::isnan(maxj[g[i]-1])) maxj[g[i]-1] = column[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(maxj, column);
        max[j] = maxj;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector maxj(ng, R_NegInf); // -INFINITY
        int ngs = 0;
        for(int i = 0; i != gss; ++i) {
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
        SHALLOW_DUPLICATE_ATTRIB(maxj, column);
        max[j] = maxj;
      }
    }
    DUPLICATE_ATTRIB(max, x);
    max.attr("row.names") = NumericVector::create(NA_REAL, -ng);
    return max;
  }
}


// Previous Version with return !!
// // [[Rcpp::export]]
// List fgmaxlCpp(List x, int ng = 0, IntegerVector g = 0, 
//                bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.size(); 
//   List max(l); 
//   
//   if (ng == 0) { 
//     if(narm) {
//       int row = 0; 
//       int k = 0;
//       double maxi = 0;
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size(); 
//         k = row-1;
//         maxi = column[k];
//         while(std::isnan(maxi) && k!=0) maxi = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(maxi < column[i]) maxi = column[i];
//         }
//         if(ret == 1) {
//           if(fill) max[j] = rep(maxi, row);
//           else {
//             NumericVector sgj(row, maxi);
//             for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i];
//             max[j] = sgj;
//           }
//         } else if (ret == 2) {
//           if(fill) max[j] = column - maxi; 
//           else {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i];
//               else sgj[i] = column[i] - maxi; 
//             }
//             max[j] = sgj;
//           }
//         } else max[j] = maxi;
//       }
//     } else {
//       double maxi = 0;
//       int row = 0;
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         maxi = column[0];
//         row = column.size();
//         for(int i = 0; i != row; ++i) { 
//           if(std::isnan(column[i])) {
//             maxi = column[i]; 
//             break;
//           } else { 
//             if(maxi < column[i]) maxi = column[i];
//           }
//         }
//         if(ret == 1) { 
//           max[j] = rep(maxi, row);
//         } else if (ret == 2) {
//           max[j] = column - maxi;
//         } else max[j] = maxi;
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
//         NumericVector maxj(ng, NA_REAL);
//         for(int i = row; i--; ) {
//           if(maxj[g[i]-1] < column[i] || std::isnan(maxj[g[i]-1])) maxj[g[i]-1] = column[i]; 
//         }
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           if(fill) { 
//             for(int i = row; i--; ) sgj[i] = maxj[g[i]-1]; 
//           } else {
//             for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i]; 
//               else sgj[i] = maxj[g[i]-1];
//             }
//           }
//           max[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row); 
//           if(fill) { 
//             for(int i = row; i--; ) sgj[i] = column[i] - maxj[g[i]-1]; 
//           } else {
//             for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i]; 
//               else sgj[i] = column[i] - maxj[g[i]-1];
//             }
//           }
//           max[j] = sgj;
//         } else max[j] = maxj;
//       }
//     } else { 
//       int row = 0;
//       int ngs = 0;
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector maxj(ng, INFINITY);
//         ngs = 0;
//         for(int i = 0; i != row; ++i) {
//           if(std::isnan(column[i])) {
//             if(!std::isnan(maxj[g[i]-1])) {
//               maxj[g[i]-1] = column[i]; 
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           } else { 
//             if(maxj[g[i]-1] < column[i]) maxj[g[i]-1] = column[i];
//           }
//         }
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = maxj[g[i]-1]; 
//           max[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = column[i] - maxj[g[i]-1];
//           max[j] = sgj;
//         } else max[j] = maxj;
//       }
//     }
//   }
//   return max;  
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
// List fgmaxlcpp(List x, IntegerVector g = 0, int ng = 0, IntegerVector gs = 0, bool narm = true, int ret = 0, bool fill = false) { // void
//   int l = x.size();
//   List max(l);
//   
//   if (ng == 0) { // No groups
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double maxi = R_NegInf;
//           int row = column.size();
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             if(maxi<column[i]) maxi = column[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = maxi; // faster with loop??: for(int i = row; i--; ) sgj[i] = maxj[g[i]-1];
//             max[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - maxi;
//             max[j] = sgj;
//           } else max[j] = maxi;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double maxi = R_NegInf;
//           int row = column.size();
//           LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             isnan[i] = ISNAN(column[i]);
//             if(isnan[i]) continue; // faster way??
//             if(maxi<column[i]) maxi = column[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = maxi; 
//             }
//             max[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = column[i] - maxi; 
//             }
//             max[j] = sgj;
//           } else max[j] = maxi;
//         }        
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double maxi = R_NegInf;
//         int row = column.size();
//         for(int i = row; i--; ) if(maxi<column[i]) maxi = column[i];
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = maxi; // faster with loop??: for(int i = row; i--; ) sgj[i] = maxj[g[i]-1];
//           max[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = column[i] - maxi;
//           max[j] = sgj;
//         } else max[j] = maxi;
//       }
//     }
//   } else { // With groups
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector maxj(ng, R_NegInf);
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             if(maxj[g[i]-1]<column[i]) maxj[g[i]-1] = column[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = maxj[g[i]-1]; // Fastest way?? better loop through cols??
//             max[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - maxj[g[i]-1];
//             max[j] = sgj;
//           } else max[j] = maxj;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector maxj(ng, R_NegInf);
//           LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             isnan[i] = ISNAN(column[i]);
//             if(isnan[i]) continue; // faster way??
//             if(maxj[g[i]-1]<column[i]) maxj[g[i]-1] = column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = maxj[g[i]-1]; 
//             }
//             max[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = column[i] - maxj[g[i]-1]; 
//             }
//             max[j] = sgj;
//           } else max[j] = maxj;
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         NumericVector maxj(ng, R_NegInf);
//         for(int i = row; i--; ) if(maxj[g[i]-1]<column[i]) maxj[g[i]-1] = column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = maxj[g[i]-1]; // Fastest way?? better loop through cols??
//           max[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = column[i] - maxj[g[i]-1];
//           max[j] = sgj;
//         } else max[j] = maxj;
//       }
//     }
//   }
//   return max;
// }