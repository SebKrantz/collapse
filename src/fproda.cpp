#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fprodmCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0,  
               bool narm = true, bool drop = true) { 
  int l = x.nrow(), col = x.ncol(); 
  
  if(ng == 0) { 
    NumericVector prod = no_init_vector(col); // Initialize faster -> Nope !!!
    if(narm) { 
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j); 
        int k = l-1;
        long double prodj = column[k]; 
        while(std::isnan(prodj) && k!=0) prodj = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(!std::isnan(column[i])) prodj *= column[i]; 
        }
        prod[j] = (double)prodj;
      }
    } else {
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        long double prodj = 1;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(column[i])) {
            prodj = column[i]; 
            break;
          } else { 
            prodj *= column[i];
          }
        }
        prod[j] = (double)prodj;
      }
    }
    if(drop) prod.attr("names") = colnames(x); 
    else {
      prod.attr("dim") = Dimension(1, col);
      colnames(prod) = colnames(x); 
    }
    return prod;
  } else { // with groups 
    if(g.size() != l) stop("length(g) must match nrow(X)");
    if(narm) {
      NumericMatrix prod = no_init_matrix(ng, col);
      std::fill(prod.begin(), prod.end(), NA_REAL); 
      for(int j = col; j--; ) { 
        NumericMatrix::ConstColumn column = x( _ , j); 
        NumericMatrix::Column prodj = prod( _ , j); 
        for(int i = l; i--; ) { 
          if(!std::isnan(column[i])) { 
            if(std::isnan(prodj[g[i]-1])) prodj[g[i]-1] = column[i];
            else prodj[g[i]-1] *= column[i];
          }
        }
      }
      colnames(prod) = colnames(x);
      return prod;
    } else {
      NumericMatrix prod = no_init_matrix(ng, col); // no init numerically unstable !!!
      std::fill(prod.begin(), prod.end(), 1.0); 
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j); 
        NumericMatrix::Column prodj = prod( _ , j);
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
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
      }
      colnames(prod) = colnames(x);
      return prod;
    }
  }
}



// Previous Version: No efficient algorithm and with return
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// SEXP fgprodmCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector prod(col, 1.0); 
//     if(narm) { 
//       for(int i = l; i--; ) { 
//         for(int j = col; j--; ) { 
//           if(std::isnan(x(i,j))) continue; 
//           prod[j] *= x(i,j); 
//         } 
//       } 
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) for(int j = col; j--; ) out(_,j) = rep(prod[j], l); 
//         else for(int i = l; i--; ) { 
//           for(int j = col; j--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//             else out(i,j) = prod[j];
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col); 
//         if(fill) for(int j = col; j--; ) out(_,j) = x(_,j) - prod[j]; 
//         else for(int i = l; i--; ) { 
//           for(int j = col; j--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//             else out(i,j) = x(i,j) - prod[j];
//           }
//         }
//         return out;
//       } else return prod;
//     } else {
//       for(int i = l; i--; ) {
//         for(int j = col; j--; ) prod[j] *= x(i,j); 
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = rep(prod[j], l);
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = x(_,j) - prod[j];
//         return out;
//       } else return prod;
//     }
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericMatrix prod(ng, col);
//     std::fill(prod.begin(), prod.end(), 1.0); // https://stackoverflow.com/questions/23748572/initializing-a-matrix-to-na-in-rcpp
//     if(narm) {
//       for(int i = l; i--; ) {
//         for(int j = col; j--; ) {
//           if(std::isnan(x(i,j))) continue; 
//           prod(g[i]-1,j) *= x(i,j); 
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) for(int i = l; i--; ) { 
//           for(int j = col; j--; ) out(i,j) = prod(g[i]-1,j);
//         }
//         else for(int i = l; i--; ) { 
//           for(int j = col; j--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//             else out(i,j) = prod(g[i]-1,j);
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         if(fill) for(int i = l; i--; ) { 
//           for(int j = col; j--; ) out(i,j) = x(i,j) - prod(g[i]-1,j);
//         }
//         else for(int i = l; i--; ) { 
//           for(int j = col; j--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//             else out(i,j) = x(i,j) - prod(g[i]-1,j);
//           }
//         }
//         return out;
//       } else return prod;
//     } else {
//       for(int i = l; i--; ) {
//         for(int j = col; j--; ) prod(g[i]-1,j) *= x(i,j); 
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int i = l; i--; ) { 
//           for(int j = col; j--; ) out(i,j) = prod(g[i]-1,j);
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int i = l; i--; ) { 
//           for(int j = col; j--; ) out(i,j) = x(i,j) - prod(g[i]-1,j);
//         }
//         return out;
//       } else return prod;
//     }
//   }
// }
