#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fprodlCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
               bool narm = true, bool drop = true) {
  int l = x.size();
  
  if(ng == 0) {
    NumericVector prod(l); 
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        int k = column.size()-1;
        long double prodi = column[k];
        while(std::isnan(prodi) && k!=0) prodi = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(!std::isnan(column[i])) prodi *= column[i];
        }
        prod[j] = (double)prodi;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        long double prodi = 1;
        int row = column.size();
        for(int i = 0; i != row; ++i) {
          if(std::isnan(column[i])) {
            prodi = column[i];
            break;
          } else {
            prodi *= column[i];
          }
        }
        prod[j] = (double)prodi;
      }
    }
    if(drop) {
      prod.attr("names") = x.attr("names");
      return prod;
    } else {
      List out(l);
      for(int j = l; j--; ) {
        out[j] = prod[j];
        SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
      }
      DUPLICATE_ATTRIB(out, x);
      out.attr("row.names") = 1;
      return out;
    }
  } else { // With groups 
    List prod(l);
    int gss = g.size();
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector prodj(ng, NA_REAL);
        for(int i = gss; i--; ) {
          if(!std::isnan(column[i])) { 
            if(std::isnan(prodj[g[i]-1])) prodj[g[i]-1] = column[i];
            else prodj[g[i]-1] *= column[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(prodj, column);
        prod[j] = prodj;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector prodj(ng, 1.0); 
        int ngs = 0;
        for(int i = 0; i != gss; ++i) {
          if(std::isnan(column[i])) {
            if(!std::isnan(prodj[g[i]-1])) {
              prodj[g[i]-1] = column[i];
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            prodj[g[i]-1] *= column[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(prodj, column);
        prod[j] = prodj;
      }
    }
    DUPLICATE_ATTRIB(prod, x);
    prod.attr("row.names") = NumericVector::create(NA_REAL, -ng);
    return prod;
  }
}



// Previous Version: No efficient algorithm and with return
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// List fgprodlCpp(List x, int ng = 0, IntegerVector g = 0, 
//                bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.size(); 
//   List prod(l); 
//   
//   if (ng == 0) { 
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double prodi = 1.0;
//         int row = column.size();
//         for(int i = row; i--; ) {
//           if(std::isnan(column[i])) continue; 
//           prodi *= column[i];
//         }
//         if(ret == 1) {
//           if(fill) prod[j] = rep(prodi, row);
//           else {
//             NumericVector sgj(row, prodi);
//             for(int i = row; i--; ) if(std::isnan(column[i])) sgj[i] = column[i];
//             prod[j] = sgj;
//           }
//         } else if (ret == 2) {
//           if(fill) prod[j] = column - prodi; 
//           else {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i];
//               else sgj[i] = column[i] - prodi; 
//             }
//             prod[j] = sgj;
//           }
//         } else prod[j] = prodi;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double prodi = 1.0;
//         int row = column.size();
//         for(int i = row; i--; ) prodi *= column[i];
//         if(ret == 1) { 
//           prod[j] = rep(prodi, row);
//         } else if (ret == 2) {
//           prod[j] = column - prodi;
//         } else prod[j] = prodi;
//       }
//     }
//   } else { // With groups !!
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector prodj(ng, 1.0);
//         for(int i = row; i--; ) {
//           if(std::isnan(column[i])) continue; 
//           prodj[g[i]-1] *= column[i];  
//         }
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           if(fill) for(int i = row; i--; ) sgj[i] = prodj[g[i]-1]; 
//           else for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i]; 
//             else sgj[i] = prodj[g[i]-1];
//           }
//           prod[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row); 
//           if(fill) for(int i = row; i--; ) sgj[i] = column[i] - prodj[g[i]-1]; 
//           else for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i]; 
//             else sgj[i] = column[i] - prodj[g[i]-1];
//           }
//           prod[j] = sgj;
//         } else prod[j] = prodj;
//       }
//     } else { 
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector prodj(ng, 1.0);
//         for(int i = row; i--; ) prodj[g[i]-1] *= column[i]; 
//         if(ret == 1) { 
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = prodj[g[i]-1]; 
//           prod[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = column[i] - prodj[g[i]-1];
//           prod[j] = sgj;
//         } else prod[j] = prodj;
//       }
//     }
//   }
//   return prod;  
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
// List fgprodlcpp(List x, IntegerVector g = 0, int ng = 0, IntegerVector gs = 0, bool narm = true, int ret = 0, bool fill = false) { // void
//   int l = x.size();
//   List prod(l);
// 
//   if (ng == 0) { // No groups
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double prodi = 0;
//           int row = column.size();
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             prodi *= column[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = prodi; // faster with loop??: for(int i = row; i--; ) sgj[i] = prodj[g[i]-1];
//             prod[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - prodi;
//             prod[j] = sgj;
//           } else prod[j] = prodi;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           double prodi = 0;
//           int row = column.size();
//           LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             isnan[i] = ISNAN(column[i]);
//             if(isnan[i]) continue; // faster way??
//             prodi *= column[i];
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = prodi; 
//             }
//             prod[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = column[i] - prodi; 
//             }
//             prod[j] = sgj;
//           } else prod[j] = prodi;
//         }        
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         double prodi = 0;
//         int row = column.size();
//         for(int i = row; i--; ) prodi *= column[i];
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = prodi; // faster with loop??: for(int i = row; i--; ) sgj[i] = prodj[g[i]-1];
//           prod[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = column[i] - prodi;
//           prod[j] = sgj;
//         } else prod[j] = prodi;
//       }
//     }
//   } else { // With groups
//     if(narm) {
//       if(fill) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector prodj(ng);
//           for(int i = row; i--; ) {
//             if(ISNAN(column[i])) continue; // faster way??
//             prodj[g[i]-1] *= column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) sgj[i] = prodj[g[i]-1]; // Fastest way?? better loop through cols??
//             prod[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; )  sgj[i] = column[i] - prodj[g[i]-1];
//             prod[j] = sgj;
//           } else prod[j] = prodj;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           NumericVector prodj(ng);
//           LogicalVector isnan(row);
//           for(int i = row; i--; ) {
//             isnan[i] = ISNAN(column[i]);
//             if(isnan[i]) continue; // faster way??
//             prodj[g[i]-1] *= column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//           }
//           if(ret == 1) { // possibility of reducing the number of passes??
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = prodj[g[i]-1]; 
//             }
//             prod[j] = sgj;
//           } else if (ret == 2) {
//             NumericVector sgj(row);
//             for(int i = row; i--; ) {
//               if(isnan[i]) sgj[i] = column[i];  // faster way?? vector subsetting??
//               else sgj[i] = column[i] - prodj[g[i]-1]; 
//             }
//             prod[j] = sgj;
//           } else prod[j] = prodj;
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         NumericVector prodj(ng);
//         for(int i = row; i--; ) prodj[g[i]-1] *= column[i]; // x[i][j]; // could think of array implementation, but I think it is good. 
//         if(ret == 1) { // possibility of reducing the number of passes??
//           NumericVector sgj(row);
//           for(int i = row; i--; ) sgj[i] = prodj[g[i]-1]; // Fastest way?? better loop through cols??
//           prod[j] = sgj;
//         } else if (ret == 2) {
//           NumericVector sgj(row);
//           for(int i = row; i--; )  sgj[i] = column[i] - prodj[g[i]-1];
//           prod[j] = sgj;
//         } else prod[j] = prodj;
//       }
//     }
//   }
//   return prod;
// }