#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fsummCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0, // No speed loss by putting SEXP !! 
              bool narm = true, bool drop = true) { 
  int l = x.nrow(), col = x.ncol(); 
  
  if(ng == 0) { 
    NumericVector sum = no_init_vector(col); // Initialize faster -> Nope !!!
    if(narm) { 
      for(int j = col; j--; ) { // Instead Am(j,_) you can use Am.row(j).
        NumericMatrix::ConstColumn column = x( _ , j); 
        int k = l-1;
        long double sumj = column[k]; // Slight speed loss 38 vs 32 milliseconds on WDIM
        while(std::isnan(sumj) && k!=0) sumj = column[--k];
        if(k != 0) for(int i = k; i--; ) {
          if(!std::isnan(column[i])) sumj += column[i]; 
        }
        sum[j] = (double)sumj; // No speed loss, but more secure
      }
    } else {
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        long double sumj = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(column[i])) {
            sumj = column[i]; 
            break;
          } else { 
            sumj += column[i];
          }
        }
        sum[j] = (double)sumj;
      }
    }
    if(drop) sum.attr("names") = colnames(x); // Slight speed loss 31 to 34 milliseconds on WDIM, but doing it in R not faster !!
    else {
      // SHALLOW_DUPLICATE_ATTRIB(sum, x);
      sum.attr("dim") = Dimension(1, col);
      colnames(sum) = colnames(x); 
    }
    return sum;
  } else { // with groups 
    if(g.size() != l) stop("length(g) must match nrow(X)");
    if(narm) {
      NumericMatrix sum = no_init_matrix(ng, col);
      std::fill(sum.begin(), sum.end(), NA_REAL); 
      for(int j = col; j--; ) { 
        NumericMatrix::ConstColumn column = x( _ , j); 
        NumericMatrix::Column sumj = sum( _ , j); 
        for(int i = l; i--; ) { 
          if(!std::isnan(column[i])) { 
            if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
            else sumj[g[i]-1] += column[i];
          }
        }
      }
      colnames(sum) = colnames(x);  // extremely efficient !!
      return sum;
    } else {
      NumericMatrix sum(ng, col); // no init numerically unstable !!!
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j); 
        NumericMatrix::Column sumj = sum( _ , j);
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
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
      }
      colnames(sum) = colnames(x); 
      return sum;
    }
    // sum.attr("dimnames") = List::create(R_NilValue,colnames(x));
    // colnames(sum) = colnames(x); // NEW! faster than R ?? -> nope, slower !! -> Do as above !!
  }
}



// // FINAL CODE with return built - in 
// // [[Rcpp::export]]
// SEXP fgsummCpp2(NumericMatrix x, int ng = 0, IntegerVector g = 0,
//                bool narm = true, int ret = 0, bool drop = true) {
//   int l = x.nrow();
//   int col = x.ncol();
// 
//   if(ng == 0) {
//     NumericVector sum = no_init_vector(col);
//     if(narm) {
//       for(int j = col; j--; ) { //Instead Am(j,_) you can use Am.row(j).
//         NumericMatrix::Column column = x( _ , j);
//         int k = l-1;
//         double sumj = column[k];
//         while(std::isnan(sumj) && k!=0) sumj = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(!std::isnan(column[i])) sumj += column[i];
//         }
//         sum[j] = sumj;
//       }
//     } else {
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         double sumj = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(column[i])) {
//             sumj = column[i];
//             break;
//           } else {
//             sumj += column[i];
//           }
//         }
//         sum[j] = sumj;
//       }
//     }
//     switch(ret) {
//     case 1:
//     {
//       NumericMatrix out = no_init_matrix(l, col);
//       for(int j = col; j--; ) out(_,j) = rep(sum[j], l);
//       return out;
//     }
//     case 2:
//     {
//       NumericMatrix out = no_init_matrix(l, col);
//       for(int j = col; j--; ) {
//         NumericMatrix::Column colo = out( _ , j);
//         NumericMatrix::Column column = x( _ , j);
//         double sumj = sum[j];
//         for(int i = l; i--; ) {
//           if(std::isnan(column[i])) colo[i] = column[i];
//           else colo[i] = sumj;
//         }
//       }
//       return out;
//     }
//     case 3:
//     {
//       NumericMatrix out = no_init_matrix(l, col);
//       for(int j = col; j--; ) {
//         NumericMatrix::Column colo = out( _ , j);
//         NumericMatrix::Column column = x( _ , j);
//         colo = column - sum[j];
//       }
//       return out;
//     }
//     default:
//       if(drop) sum.attr("names") = colnames(x); // Slight speed loss 31 to 34 milliseconds on WDIM, but doing it in R not faster !!
//       else sum.attr("dim") = Dimension(1, col);
//       return sum;
//     }
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericMatrix sum = no_init_matrix(ng, col);
//     if(narm) {
//       std::fill(sum.begin(), sum.end(), NA_REAL);
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         NumericMatrix::Column sumj = sum( _ , j);
//         for(int i = l; i--; ) {
//           if(!std::isnan(column[i])) {
//             if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//             else sumj[g[i]-1] += column[i];
//           }
//         }
//       }
//     } else {
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         NumericMatrix::Column sumj = sum( _ , j);
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
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
//       }
//     }
//     switch(ret) {
//     case 1:
//     {
//       NumericMatrix out = no_init_matrix(l, col);
//       for(int j = col; j--; ) {
//         NumericMatrix::Column colo = out( _ , j);
//         NumericMatrix::Column sumj = sum( _ , j);
//         for(int i = l; i--; ) colo[i] = sumj[g[i]-1];
//       }
//       return out;
//     }
//     case 2:
//     {
//       NumericMatrix out = no_init_matrix(l, col);
//       for(int j = col; j--; ) {
//         NumericMatrix::Column colo = out( _ , j);
//         NumericMatrix::Column column = x( _ , j);
//         NumericMatrix::Column sumj = sum( _ , j);
//         for(int i = l; i--; ) {
//           if(std::isnan(column[i])) colo[i] = column[i];
//           else colo[i] = sumj[g[i]-1];
//         }
//       }
//       return out;
//     }
//     case 3:
//     {
//       NumericMatrix out = no_init_matrix(l, col);
//       for(int j = col; j--; ) {
//         NumericMatrix::Column colo = out( _ , j);
//         NumericMatrix::Column column = x( _ , j);
//         NumericMatrix::Column sumj = sum( _ , j);
//         for(int i = l; i--; ) colo[i] = column[i] - sumj[g[i]-1];
//       }
//       return out;
//     }
//     default:
//       return sum;
//     }
//   }
// }



// FINAL CODE: No COMMENTS, BUT EXTERIMENTING A BIT WITH MATRIX MEMORY OPTIMIZATION. OPTIMIZED CODE ABOVE
// // [[Rcpp::export]]
// SEXP fgsummCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                bool narm = true, int ret = 0, bool drop = true) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector sum(col); 
//     if(narm) { 
//       for(int j = col; j--; ) { 
//         NumericMatrix::Column column = x( _ , j); 
//         int k = l-1;
//         double sumj = column[k]; 
//         while(std::isnan(sumj) && k!=0) sumj = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(!std::isnan(column[i])) sumj += column[i]; 
//         }
//         sum[j] = sumj;
//       }
//     } else {
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         double sumj = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(column[i])) {
//             sumj = column[i]; 
//             break;
//           } else { 
//             sumj += column[i];
//           }
//         }
//         sum[j] = sumj;
//       }
//     }
//     switch(ret) {
//     case 1: 
//     {
//       NumericMatrix out(l, col);
//       for(int j = col; j--; ) out(_,j) = rep(sum[j], l);  
//       return out;
//     }
//     case 2: 
//     {
//       NumericMatrix out(l, col); 
//       for(int j = col; j--; ) {
//         NumericMatrix::Column colo = out( _ , j);
//         NumericMatrix::Column column = x( _ , j);
//         double sumj = sum[j];
//         for(int i = l; i--; ) {
//           if(std::isnan(column[i])) colo[i] = column[i];
//           else colo[i] = sumj;
//         }
//       }
//       return out;
//       
//       // NumericMatrix out;
//       // for(int j = col; j--; ) {
//       //   NumericVector colo(l);
//       //   NumericMatrix::Column column = x( _ , j);
//       //   double sumj = sum[j];
//       //   for(int i = l; i--; ) {
//       //     if(std::isnan(column[i])) colo[i] = column[i];
//       //     else colo[i] = sumj;
//       //   }
//       //   out = cbind(out,colo);
//       // }
//       // return out;
//       
//       // // NumericMatrix out(l, col);
//       // ListMatrix out(l, col); // list is a lot faster !! 100 vs 140 milliseconds on WDIM !!
//       // for(int j = col; j--; ) {
//       //   // NumericMatrix::Column colo = out( _ , j);
//       //   NumericVector colo(l);
//       //   NumericMatrix::Column column = x( _ , j);
//       //   double sumj = sum[j];
//       //   for(int i = l; i--; ) {
//       //     if(std::isnan(column[i])) colo[i] = column[i];
//       //     else colo[i] = sumj;
//       //   }
//       //   out[j] = colo;
//       // }
//       // return out;
//       
//       
//       // NumericVector out(l*col); // Also slow !! It is this memory allocation that takes time !!
//       // for(int j = col; j--; ) { 
//       //   NumericMatrix::Column column = x( _ , j);
//       //   double sumj = sum[j]; 
//       //   int k = j*col;
//       //   for(int i = l; i--; ) {
//       //     if(std::isnan(column[i])) out[k+i] = column[i];
//       //     else out[k+i] = sumj;
//       //   }
//       // }
//       // out.attr("dim") = Dimension(l, col);
//       // return out;
//       // 
//     }
//     case 3: 
//     {
//       NumericMatrix out(l, col); 
//       for(int j = col; j--; ) { 
//         NumericMatrix::Column colo = out( _ , j);
//         NumericMatrix::Column column = x( _ , j);
//         colo = column - sum[j]; 
//       }
//       return out;
//     }
//     default: 
//       if(drop) sum.attr("names") = colnames(x); // Slight speed loss 31 to 34 milliseconds on WDIM, but doing it in R not faster !!
//       else sum.attr("dim") = Dimension(1, col);
//       return sum;
//     } 
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericMatrix sum(ng, col);
//     if(narm) {
//       std::fill(sum.begin(), sum.end(), NA_REAL); 
//       for(int j = col; j--; ) { 
//         NumericMatrix::Column column = x( _ , j); 
//         NumericMatrix::Column sumj = sum( _ , j); 
//         for(int i = l; i--; ) { 
//           if(!std::isnan(column[i])) { 
//             if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//             else sumj[g[i]-1] += column[i];
//           }
//         }
//       }
//     } else {
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j); 
//         NumericMatrix::Column sumj = sum( _ , j);
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
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
//       }
//     }
//     switch(ret) {
//     case 1: 
//     {
//       NumericMatrix out(l, col);
//       for(int j = col; j--; ) {
//         NumericMatrix::Column colo = out( _ , j);
//         NumericMatrix::Column sumj = sum( _ , j); 
//         for(int i = l; i--; ) colo[i] = sumj[g[i]-1];
//       }
//       return out;
//     }
//     case 2: 
//     {
//       NumericMatrix out(l, col);
//       for(int j = col; j--; ) {
//         NumericMatrix::Column colo = out( _ , j); 
//         NumericMatrix::Column column = x( _ , j);
//         NumericMatrix::Column sumj = sum( _ , j);
//         for(int i = l; i--; ) {
//           if(std::isnan(column[i])) colo[i] = column[i];
//           else colo[i] = sumj[g[i]-1];
//         }
//       }
//       return out;
//     }
//     case 3: 
//     {
//       NumericMatrix out(l, col); 
//       for(int j = col; j--; ) {
//         NumericMatrix::Column colo = out( _ , j); 
//         NumericMatrix::Column column = x( _ , j);
//         NumericMatrix::Column sumj = sum( _ , j);
//         for(int i = l; i--; ) colo[i] = column[i] - sumj[g[i]-1];
//       }
//       return out;
//     }
//     default: 
//       return sum;
//     } 
//   }
// }

// Another Version of the Code: One loop (like list-version): ABout same speed, maybe a tiny bit faster, but gives stragge error on some data. 
// See also commented version of above code below !!
// // [[Rcpp::export]]
// SEXP fgsummCpp2(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                 bool narm = true, int ret = 0) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     if(ret != 0) {
//       NumericMatrix out(l, col);
//       if(narm) { 
//         for(int j = col; j--; ) { 
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column outj = out( _ , j);
//           int k = l-1;
//           double sumj = column[k]; 
//           while(std::isnan(sumj) && k!=0) sumj = column[--k];
//           if(k != 0) for(int i = k; i--; ) {
//             if(!std::isnan(column[i])) sumj += column[i]; 
//           }
//           switch(ret) {
//           case 1: 
//             std::fill(outj.begin(), outj.end(), sumj); // or std::fill_n ?? 
//             break;
//           case 2: 
//             for(int i = l; i--; ) {
//               if(std::isnan(column[i])) outj[i] = column[i];
//               else outj[i] = sumj;
//             }
//             break;
//           case 3: 
//             outj = column - sumj;
//             break;
//           }  // End of switch
//         }
//       } else {
//         for(int j = col; j--; ) { // Without seems faster than referencing column!! -> But I keep it for symmetry !!!
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column outj = out( _ , j);
//           double sumj = 0;
//           for(int i = 0; i != l; ++i) {
//             if(std::isnan(column[i])) {
//               sumj = column[i]; 
//               break;
//             } else { 
//               sumj += column[i];
//             }
//           }
//           switch(ret) {
//           case 1: 
//             std::fill(outj.begin(), outj.end(), sumj); // or std::fill_n ?? 
//             break;
//           case 2: 
//             for(int i = l; i--; ) {
//               if(std::isnan(column[i])) outj[i] = column[i];
//               else outj[i] = sumj;
//             }
//             break;
//           case 3: 
//             outj = column - sumj;
//             break;
//           }  // End of switch
//         }
//       }
//       return out;
//     } else {
//       NumericVector sum(col);
//       if(narm) { 
//         for(int j = col; j--; ) { 
//           NumericMatrix::Column column = x( _ , j);
//           int k = l-1;
//           double sumj = column[k]; 
//           while(std::isnan(sumj) && k!=0) sumj = column[--k];
//           if(k != 0) for(int i = k; i--; ) {
//             if(!std::isnan(column[i])) sumj += column[i]; 
//           }
//           sum[j] = sumj;
//         }
//       } else {
//         for(int j = col; j--; ) { // Without seems faster than referencing column!! -> But I keep it for symmetry !!!
//           NumericMatrix::Column column = x( _ , j);
//           double sumj = 0;
//           for(int i = 0; i != l; ++i) {
//             if(std::isnan(column[i])) {
//               sumj = column[i]; 
//               break;
//             } else { 
//               sumj += column[i];
//             }
//           }
//           sum[j] = sumj;
//         }
//       }
//       return sum;
//     }
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     if(ret != 0) {
//       NumericMatrix out(l, col);
//       if(narm) { 
//         for(int j = col; j--; ) { 
//           NumericMatrix::Column column = x( _ , j); 
//           NumericMatrix::Column outj = out( _ , j);
//           NumericVector sumj(col, NA_REAL); 
//           for(int i = l; i--; ) { 
//             if(!std::isnan(column[i])) { 
//               if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//               else sumj[g[i]-1] += column[i];
//             }
//           }
//           switch(ret) {
//           case 1: 
//             for(int i = l; i--; ) outj[i] = sumj[g[i]-1];
//             break;
//           case 2: 
//             for(int i = l; i--; ) {
//               if(std::isnan(column[i])) outj[i] = column[i];
//               else outj[i] = sumj[g[i]-1];
//             }
//             break;
//           case 3: 
//             for(int i = l; i--; ) outj[i] = column[i] - sumj[g[i]-1];
//             break;
//           }  // End of switch
//         }
//       } else {
//         for(int j = col; j--; ) {
//           NumericMatrix::Column column = x( _ , j); 
//           NumericMatrix::Column outj = out( _ , j);
//           NumericVector sumj(col);
//           int ngs = 0;
//           for(int i = 0; i != l; ++i) {
//             if(std::isnan(column[i])) { // fastest ?? -> Yes, with WDI matrix the alternative is slower !!!
//               if(!std::isnan(sumj[g[i]-1])) {
//                 sumj[g[i]-1] = column[i]; 
//                 ++ngs;
//                 if(ngs == ng) break;
//               }
//             } else { 
//               sumj[g[i]-1] += column[i];
//             }
//           }
//           switch(ret) {
//           case 1: 
//             for(int i = l; i--; ) outj[i] = sumj[g[i]-1];
//             break;
//           case 2: 
//             for(int i = l; i--; ) {
//               if(std::isnan(column[i])) outj[i] = column[i];
//               else outj[i] = sumj[g[i]-1];
//             }
//             break;
//           case 3: 
//             for(int i = l; i--; ) outj[i] = column[i] - sumj[g[i]-1];
//             break;
//           }  // End of switch
//         }
//       }
//       return out;
//     } else {
//       NumericMatrix sum(ng, col);
//       if(narm) { 
//         for(int j = col; j--; ) { 
//           NumericMatrix::Column column = x( _ , j); 
//           NumericMatrix::Column sumj = sum( _ , j); 
//           for(int i = l; i--; ) { 
//             if(!std::isnan(column[i])) { 
//               if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//               else sumj[g[i]-1] += column[i];
//             }
//           }
//         }
//       } else {
//         for(int j = col; j--; ) {
//           NumericMatrix::Column column = x( _ , j); 
//           NumericMatrix::Column sumj = sum( _ , j); 
//           int ngs = 0;
//           for(int i = 0; i != l; ++i) {
//             if(std::isnan(column[i])) { // fastest ?? -> Yes, with WDI matrix the alternative is slower !!!
//               if(!std::isnan(sumj[g[i]-1])) {
//                 sumj[g[i]-1] = column[i]; 
//                 ++ngs;
//                 if(ngs == ng) break;
//               }
//             } else { 
//               sumj[g[i]-1] += column[i];
//             }
//           }
//         }
//       }
//       return sum;
//     }
//   }
// }


// FINAL CODE, with Comments
// // [[Rcpp::export]]
// SEXP fgsummCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                bool narm = true, int ret = 0) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector sum(col); 
//     if(narm) { 
//       // This also takes about 35 milliseconds -> Take it, better algorithm !! Check sources of speed gains again !!
//       for(int j = col; j--; ) { // could do simpler by simple setting NA_REAL as the default value, but this code is more elegant and faster !!
//         NumericMatrix::Column column = x( _ , j); // This gives 3 milliseconds: 51 to 48 on WDIM !!!
//         int k = l-1;
//         double sumj = column[k]; // Using a double gives you 10 milliseconds: 48 to 35 on WDIM !!!
//         while(std::isnan(sumj) && k!=0) sumj = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           if(!std::isnan(column[i])) sumj += column[i]; // Still slightly faster !!
//         }
//         sum[j] = sumj;
//       }
//     } else {
//       for(int j = col; j--; ) { // Without seems faster than referencing column!! -> But I keep it for symmetry !!!
//         NumericMatrix::Column column = x( _ , j);
//         double sumj = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(column[i])) {
//             sumj = column[i]; 
//             break;
//           } else { 
//             sumj += column[i];
//           }
//         }
//         sum[j] = sumj;
//       }
//     }
//     switch(ret) {
//       case 1: 
//       {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = rep(sum[j], l); // 161 milliseconds on WDIM -> about the same.. 
//            // NumericVector out(l*col); // Idea of manually indexing !!!
//            // NumericVector ref(l);
//            // ref = seq(0,l-1);
//            // for(int j = 1; j != col+1; ++j) out[ref+l*j] = sum[j-1];
//            // out.attr("dim") = Dimension(l, col);
//
//         return out;
//       }
//       case 2: 
//       {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) { // 178 milliseconds on WDIM !!
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           double sumj = sum[j]; // Does double really bring extra speed here?? - Yes !!
//           for(int i = l; i--; ) {
//             if(std::isnan(column[i])) colo[i] = column[i];
//             else colo[i] = sumj;
//           }
//         }
//         return out;
//       }
//       case 3: 
//       {
//         NumericMatrix out(l, col); 
//         // for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j]; // 173-175 milliseconds on WDIM
//         for(int j = col; j--; ) { // Slightly faster !!
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           colo = column - sum[j]; // doesn't matter here, because no loop over i !!
//         }
//         return out;
//       }
//       default: 
//         return sum;
//     } // End of switch
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericMatrix sum(ng, col);
//     if(narm) {
//       std::fill(sum.begin(), sum.end(), NA_REAL); 
//       for(int j = col; j--; ) { 
//         NumericMatrix::Column column = x( _ , j); // 10 milliseconds gain -> down to 46 milliseconds on WDIM
//         NumericMatrix::Column sumj = sum( _ , j); // Not really big gain, maybe 1 millisecond !!
//         for(int i = l; i--; ) { 
//           if(!std::isnan(column[i])) { // faster way to code this ??? -> Not Bad at all, 54 millisec for WDIM
//             if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
//             else sumj[g[i]-1] += column[i];
//           }
//         }
//       }
//     } else {
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j); // faster !!!
//         NumericMatrix::Column sumj = sum( _ , j);
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(column[i])) { // fastest ?? -> Yes, with WDI matrix the alternative is slower !!!
//             if(!std::isnan(sumj[g[i]-1])) {
//               sumj[g[i]-1] = column[i]; 
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           } else { 
//             sumj[g[i]-1] += column[i];
//           }
//         }
//       }
//     }
//     switch(ret) {
//       case 1: 
//       {
//         NumericMatrix out(l, col);
//         // for(int i = l; i--; ) out(i, _ ) = sum(g[i]-1, _ ); // A lot slower !! -> Bad idea!!
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j); // No tangible improvement
//           NumericMatrix::Column sumj = sum( _ , j); // Also no tangible improvement -> But for consistency !!
//           for(int i = l; i--; ) colo[i] = sumj[g[i]-1];
//         }
//         return out;
//       }
//       case 2: 
//       {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j); // Gives a slight improvement 
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column sumj = sum( _ , j);
//           for(int i = l; i--; ) {
//             if(std::isnan(column[i])) colo[i] = column[i];
//             else colo[i] = sumj[g[i]-1];
//           }
//         }
//         return out;
//       }
//       case 3: 
//       {
//         NumericMatrix out(l, col); 
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j); // Also brings a slight improvement !!
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column sumj = sum( _ , j);
//           for(int i = l; i--; ) colo[i] = column[i] - sumj[g[i]-1];
//         }
//         return out;
//       }
//       default: 
//         return sum;
//     } // End of switch
//   }
// }



// Previous Version: testing and improving (introducing column reference), but not complete, only without groups !!
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// SEXP fgsummCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                bool narm = true, int ret = 0, bool fill = false, bool set = false) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector sum(col); 
//     if(narm) { 
//       // This takes about 50 milliseconds
//       // for(int j = col; j--; ) { // could do simpler by simple setting NA_REAL as the default value, but this code is more elegant and faster !!
//       //   int k = l-1;
//       //   sum[j] = x(k,j);
//       //   while(std::isnan(sum[j]) && k!=0) sum[j] = x(--k,j);
//       //   if(k != 0) for(int i = k; i--; ) {
//       //     if(!std::isnan(x(i,j))) sum[j] += x(i,j);
//       //   }
//       // }
//       
//       // This takes about 35 Milliseconds !!
//       // for(int j = col; j--; ) { // could do simpler by simple setting NA_REAL as the default value, but this code is more elegant and faster !!
//       //   double sumj = NA_REAL; // 54 to 45 milliseconds on WDIM
//       //   NumericMatrix::Column column = x( _ , j); // 45 to 34.48 !!!
//       //   for(int i = l; i--; ) {
//       //     if(!std::isnan(column[i])) { // faster way to code this ??? -> Not Bad at all, 54 millisec for WDIM
//       //       if(std::isnan(sumj)) sumj = column[i];
//       //       else sumj += column[i];
//       //     }
//       //   }
//       //   sum[j] = sumj;
//       // }
// 
//       // This also takes about 35 milliseconds -> Take it, better algorithm !! Check sources of speed gains again !!
//       for(int j = col; j--; ) { // could do simpler by simple setting NA_REAL as the default value, but this code is more elegant and faster !!
//         NumericMatrix::Column column = x( _ , j); // This gives 3 milliseconds: 51 to 48 on WDIM !!!
//         int k = l-1;
//         double sumj = column[k]; // Using a double gives you 10 milliseconds: 48 to 35 on WDIM !!!
//         while(std::isnan(sumj) && k!=0) sumj = column[--k];
//         if(k != 0) for(int i = k; i--; ) {
//           // if(std::isnan(column[i])) continue;
//           // sumj += column[i];
//           if(!std::isnan(column[i])) sumj += column[i]; // Still slightly faster !!
//         }
//         sum[j] = sumj;
//       }
//       
//       switch(ret) {
//       case 1: 
//        if(set) {
//          if(fill) {
//            for(int j = col; j--; ) {
//              NumericMatrix::Column column = x( _ , j); 
//              std::fill(column.begin(), column.end(), sum[j]); // ... milliseconds on WDIM !!!
//            }
//          } else {
//            for(int j = col; j--; ) { // 85 milliseconds on WDIM !!
//              NumericMatrix::Column column = x( _ , j);
//              double sumj = sum[j]; 
//              for(int i = l; i--; ) {
//                if(!std::isnan(column[i])) column[i] = sumj;
//              }
//            }
//          }
//          return R_NilValue;
//        } else {
//         NumericMatrix out(l, col);
//         if(fill) {
//           // for(int j = col; j--; ) out(_,j) = rep(sum[j], l); // 161 milliseconds on WDIM -> about the same.. 
//            for(int j = col; j--; ) {
//              NumericMatrix::Column colo = out( _ , j); 
//           //   colo = rep(sum[j], l); // 161 milliseconds on WDIM !! -> slightly slower 
//              std::fill(colo.begin(), colo.end(), sum[j]); // 160 milliseconds on WDIM !!!
//           // 
//           }
//         } else {
//           // for(int j = col; j--; ) { // 200 milliseconds on WDIM !!
//           //   for(int i = l; i--; ) {
//           //     if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//           //     else out(i,j) = sum[j];
//           //   }
//           // }
//           for(int j = col; j--; ) { // 178 milliseconds on WDIM !!
//             NumericMatrix::Column colo = out( _ , j);
//             NumericMatrix::Column column = x( _ , j);
//             double sumj = sum[j]; // Does double really bring extra speed here??
//             for(int i = l; i--; ) {
//               if(std::isnan(column[i])) colo[i] = column[i];
//               else colo[i] = sumj;
//             }
//           }
//         }
//         return out;
//        }
//        // break; // Break needed after return statement ?? -> Nope!! https://stackoverflow.com/questions/6330114/do-you-need-break-in-switch-when-return-is-used
//       case 2: 
//         if(set) {
//         if(fill) {
//           for(int j = col; j--; ) {
//             NumericMatrix::Column column = x( _ , j); 
//             column = column - sum[j]; // ... milliseconds on WDIM !!!
//           }
//         } else {
//           for(int j = col; j--; ) { // ... milliseconds on WDIM !!
//             NumericMatrix::Column column = x( _ , j);
//             double sumj = sum[j]; 
//             for(int i = l; i--; ) {
//               if(!std::isnan(column[i])) column[i] -= sumj;
//             }
//           }
//         }
//         return R_NilValue;
//         } else {
//         NumericMatrix out(l, col); 
//         if(fill) {
//           // for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j]; // 173-175 milliseconds on WDIM
//           for(int j = col; j--; ) {
//             NumericMatrix::Column colo = out( _ , j);
//             colo = colo - sum[j]; // 164 milliseconds on WDIM !!
//           }
//         } else {
//           // for(int j = col; j--; ) { // 203 milliseconds on WDIM !!
//           //   for(int i = l; i--; ) {
//           //     if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//           //     else out(i,j) = x(i,j) - sum[j];
//           //   }
//           // }
//           for(int j = col; j--; ) { // 180 milliseconds on WDIM !!
//             NumericMatrix::Column colo = out( _ , j);
//             NumericMatrix::Column column = x( _ , j);
//             double sumj = sum[j];
//             for(int i = l; i--; ) {
//               if(std::isnan(column[i])) colo[i] = column[i];
//               else colo[i] = column[i] - sumj;
//             }
//           }
//         }
//         return out;
//        }
//         // break; // Break needed after return statement ?? -> Nope!! https://stackoverflow.com/questions/6330114/do-you-need-break-in-switch-when-return-is-used
//       default: return sum;
//      } // End of switch
//     } else {
//       for(int j = col; j--; ) {
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(x(i,j))) {
//             sum[j] = x(i,j); 
//             break;
//           } else { 
//             sum[j] += x(i,j);
//           }
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = rep(sum[j], l);
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j];
//         return out;
//       } else return sum;
//     }
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericMatrix sum(ng, col);
//     if(narm) {
//       std::fill(sum.begin(), sum.end(), NA_REAL); 
//       for(int j = col; j--; ) { 
//         for(int i = l; i--; ) { 
//           if(!std::isnan(x(i,j))) { // faster way to code this ??? -> Not Bad at all, 54 millisec for WDIM
//             if(std::isnan(sum(g[i]-1,j))) sum(g[i]-1,j) = x(i,j);
//             else sum(g[i]-1,j) += x(i,j);
//           }
//           // if(std::isnan(x(i,j))) {
//           //   sum(g[i]-1,j) += x(i,j);
//           // } else {
//           //   sum(g[i]-1,j) += x(i,j);
//           // }
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) out(i,j) = sum(g[i]-1,j);
//           }
//         } else { 
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//               else out(i,j) = sum(g[i]-1,j);
//             }
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) out(i,j) = x(i,j) - sum(g[i]-1,j);
//           }
//         } else { 
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//               else out(i,j) = x(i,j) - sum(g[i]-1,j);
//             }
//           }
//         }
//         return out;
//       } else return sum;
//     } else {
//       for(int j = col; j--; ) {
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(x(i,j))) { // fastest ?? -> Yes, with WDI matrix the alternative is slower !!!
//             if(!std::isnan(sum(g[i]-1,j))) {
//               sum(g[i]-1,j) = x(i,j); 
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           } else { 
//             sum(g[i]-1,j) += x(i,j);
//           }
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) { 
//           for(int i = l; i--; ) out(i,j) = sum(g[i]-1,j);
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) { 
//           for(int i = l; i--; ) out(i,j) = x(i,j) - sum(g[i]-1,j);
//         }
//         return out;
//       } else return sum;
//     }
//   }
// }


// Previous Version: fast adjusted algorithms exchanged loop order !!
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// SEXP fgsummCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                bool narm = true, int ret = 0, bool fill = false, bool set = false) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector sum(col); 
//     if(narm) { 
//       int k = 0;
//       for(int j = col; j--; ) {
//         k = l-1;
//         sum[j] = x(k,j);
//         while(std::isnan(sum[j]) && k!=0) sum[j] = x(--k,j);
//         if(k != 0) for(int i = k; i--; ) {
//           if(!std::isnan(x(i,j))) sum[j] += x(i,j);
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) out(_,j) = rep(sum[j], l); 
//         } else {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//               else out(i,j) = sum[j];
//             }
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col); 
//         if(fill) {
//           for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j]; 
//         } else {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//               else out(i,j) = x(i,j) - sum[j];
//             }
//           }
//         }
//         return out;
//       } else return sum;
//     } else {
//       for(int j = col; j--; ) {
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(x(i,j))) {
//             sum[j] = x(i,j); 
//             break;
//           } else { 
//             sum[j] += x(i,j);
//           }
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = rep(sum[j], l);
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j];
//         return out;
//       } else return sum;
//     }
//       } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericMatrix sum(ng, col);
//     if(narm) {
//       std::fill(sum.begin(), sum.end(), NA_REAL); 
//       for(int j = col; j--; ) { 
//         for(int i = l; i--; ) { 
//           if(!std::isnan(x(i,j))) { // faster way to code this ??? -> Not Bad at all, 54 millisec for WDIM
//             if(std::isnan(sum(g[i]-1,j))) sum(g[i]-1,j) = x(i,j);
//             else sum(g[i]-1,j) += x(i,j);
//           }
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) out(i,j) = sum(g[i]-1,j);
//           }
//         } else { 
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//               else out(i,j) = sum(g[i]-1,j);
//             }
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         if(fill) {
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) out(i,j) = x(i,j) - sum(g[i]-1,j);
//           }
//         } else { 
//           for(int j = col; j--; ) { 
//             for(int i = l; i--; ) {
//               if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//               else out(i,j) = x(i,j) - sum(g[i]-1,j);
//             }
//           }
//         }
//         return out;
//       } else return sum;
//     } else {
//       for(int j = col; j--; ) {
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(x(i,j))) { // fastest ?? -> Yes, with WDI matrix the alternative is slower !!!
//             if(!std::isnan(sum(g[i]-1,j))) {
//               sum(g[i]-1,j) = x(i,j); 
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           } else { 
//             sum(g[i]-1,j) += x(i,j);
//           }
//         }
//       }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) { 
//           for(int i = l; i--; ) out(i,j) = sum(g[i]-1,j);
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) { 
//           for(int i = l; i--; ) out(i,j) = x(i,j) - sum(g[i]-1,j);
//         }
//         return out;
//       } else return sum;
//     }
//   }
// }


// Previous Version: Before adjusting the algorithm and exchanging loop order !!!
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// SEXP fgsummCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                  bool narm = true, int ret = 0, bool fill = false) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector sum(col); 
//     if(narm) { 
//         for(int i = l; i--; ) { 
//           for(int j = col; j--; ) { 
//             if(std::isnan(x(i,j))) continue; 
//             sum[j] += x(i,j); 
//           } 
//         } 
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) for(int j = col; j--; ) out(_,j) = rep(sum[j], l); 
//         else for(int i = l; i--; ) { 
//           for(int j = col; j--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//             else out(i,j) = sum[j];
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col); 
//         if(fill) for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j]; 
//         else for(int i = l; i--; ) { 
//           for(int j = col; j--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j);
//             else out(i,j) = x(i,j) - sum[j];
//           }
//         }
//         return out;
//       } else return sum;
//     } else {
//         for(int i = l; i--; ) {
//           for(int j = col; j--; ) sum[j] += x(i,j); 
//         }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = rep(sum[j], l);
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int j = col; j--; ) out(_,j) = x(_,j) - sum[j];
//         return out;
//       } else return sum;
//     }
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     NumericMatrix sum(ng, col);
//     if(narm) {
//         for(int i = l; i--; ) {
//           for(int j = col; j--; ) {
//             if(std::isnan(x(i,j))) continue; 
//             sum(g[i]-1,j) += x(i,j); 
//           }
//         }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         if(fill) for(int i = l; i--; ) { 
//           for(int j = col; j--; ) out(i,j) = sum(g[i]-1,j);
//         }
//         else for(int i = l; i--; ) { 
//           for(int j = col; j--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//             else out(i,j) = sum(g[i]-1,j);
//           }
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         if(fill) for(int i = l; i--; ) { 
//           for(int j = col; j--; ) out(i,j) = x(i,j) - sum(g[i]-1,j);
//         }
//         else for(int i = l; i--; ) { 
//           for(int j = col; j--; ) {
//             if(std::isnan(x(i,j))) out(i,j) = x(i,j); 
//             else out(i,j) = x(i,j) - sum(g[i]-1,j);
//           }
//         }
//         return out;
//       } else return sum;
//     } else {
//           for(int i = l; i--; ) {
//             for(int j = col; j--; ) sum(g[i]-1,j) += x(i,j); 
//           }
//       if(ret == 1) { 
//         NumericMatrix out(l, col);
//         for(int i = l; i--; ) { 
//           for(int j = col; j--; ) out(i,j) = sum(g[i]-1,j);
//         }
//         return out;
//       } else if (ret == 2) {
//         NumericMatrix out(l, col);
//         for(int i = l; i--; ) { 
//           for(int j = col; j--; ) out(i,j) = x(i,j) - sum(g[i]-1,j);
//         }
//         return out;
//       } else return sum;
//     }
//   }
// }
