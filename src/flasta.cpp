// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// with better handling row.names 

template <int RTYPE>
inline bool isnaNUM(typename Rcpp::traits::storage_type<RTYPE>::type x) { 
  return x != x;
}

template <int RTYPE>
inline bool isnaOTH(typename Rcpp::traits::storage_type<RTYPE>::type x) { 
  return x == Vector<RTYPE>::get_na();
}

template <int RTYPE>
SEXP flastmCppImpl(const Matrix<RTYPE>& x, int ng, const IntegerVector& g, bool narm, bool drop) { 
  int l = x.nrow(), col = x.ncol(); 
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;
  
  if(ng == 0) { 
    Vector<RTYPE> last = no_init_vector(col); 
    if(narm) { 
      for(int j = col; j--; ) {
        ConstMatrixColumn<RTYPE> column = x( _ , j); 
        int k = l-1;
        auto lastj = column[k];
        while(isnanT(lastj) && k!=0) lastj = column[--k];
        last[j] = lastj;
      }
    } else {
      last = x(l-1, _); 
    }
    if(drop) last.attr("names") = colnames(x); 
    else {
      last.attr("dim") = Dimension(1, col);
      colnames(last) = colnames(x); 
    }
    return last;
  } else { // with groups 
    if(g.size() != l) stop("length(g) must match nrow(X)");
    Matrix<RTYPE> last = no_init_matrix(ng, col);
    if(narm) {
      std::fill(last.begin(), last.end(), Vector<RTYPE>::get_na()); 
      for(int j = col; j--; ) { 
        ConstMatrixColumn<RTYPE> column = x( _ , j); 
        MatrixColumn<RTYPE> lastj = last( _ , j); 
        int ngs = 0;
        for(int i = l; i--; ) {
          if(!isnanT(column[i])) { 
            if(isnanT(lastj[g[i]-1])) { 
              lastj[g[i]-1] = column[i]; 
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
      }
      colnames(last) = colnames(x);
    } else {
      List dn = x.attr("dimnames");
      if(dn[0] != R_NilValue) {
        CharacterVector rn = dn[0];
        CharacterVector newrn = no_init_vector(ng);
        LogicalVector glj(ng, true); // using std::vector<bool> here is more memory efficient but not faster !!
        int ngs = 0;
        for(int i = l; i--; ) {
          if(glj[g[i]-1]) {
            glj[g[i]-1] = false;
            last(g[i]-1, _) = x(i, _);
            newrn[g[i]-1] = rn[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
        last.attr("dimnames") = List::create(newrn, dn[1]); // best way !!
      } else {
        LogicalVector glj(ng, true); // using std::vector<bool> here is more memory efficient but not faster !!
        int ngs = 0;
        for(int i = l; i--; ) {
          if(glj[g[i]-1]) {
            glj[g[i]-1] = false;
            last(g[i]-1, _) = x(i, _);
            ++ngs;
            if(ngs == ng) break;
          }
        }
        colnames(last) = colnames(x);
      }
    }
    return last;
  }
}

template <>
SEXP flastmCppImpl(const Matrix<CPLXSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP flastmCppImpl(const Matrix<VECSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP flastmCppImpl(const Matrix<RAWSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP flastmCppImpl(const Matrix<EXPRSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP flastmCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true, bool drop = true){
  RCPP_RETURN_MATRIX(flastmCppImpl, x, ng, g, narm, drop);
}





// Old code: Solving the na.rm problem using a macro and an iterator
// #define isnan2(x) ((x) != (x)) // http://www.ebyte.it/library/codesnippets/WritingCppMacros.html
// 
// template <int RTYPE>
// SEXP flastmCppImpl(Matrix<RTYPE> x, int ng, IntegerVector g, bool narm, bool drop) { 
//   int l = x.nrow(), col = x.ncol(); 
//   
//   if(ng == 0) { 
//     Vector<RTYPE> last = no_init_vector(col);
//     if(narm) { 
//       if(TYPEOF(x) == REALSXP) {
//         for(int j = col; j--; ) {
//           MatrixColumn<RTYPE> column = x( _ , j); 
//           int k = l-1;
//           auto lastj = column[k];
//           while(lastj != lastj && k!=0) lastj = column[--k];
//           last[j] = lastj;
//         }
//       } else {
//         for(int j = col; j--; ) {
//           MatrixColumn<RTYPE> column = x( _ , j); 
//           int k = l-1;
//           auto lastj = column[k];
//           while(lastj == Vector<RTYPE>::get_na() && k!=0) lastj = column[--k];
//           last[j] = lastj;
//         }
//       }
//     } else {
//       for(int j = col; j--; ) {
//         MatrixColumn<RTYPE> column = x( _ , j);
//         last[j] = column[l-1]; // Good?? What about character?? 
//       }
//     }
//     if(drop) last.attr("names") = colnames(x); 
//     else {
//       last.attr("dim") = Dimension(1, col);
//       colnames(last) = colnames(x); 
//     }
//     return last;
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     Matrix<RTYPE> last = no_init_matrix(ng, col);
//     if(narm) {
//       std::fill(last.begin(), last.end(), Vector<RTYPE>::get_na()); 
//       if(TYPEOF(x) == REALSXP) {
//         for(int j = col; j--; ) { 
//           MatrixColumn<RTYPE> column = x( _ , j); 
//           MatrixColumn<RTYPE> lastj = last( _ , j); 
//           auto fjg = lastj.begin();
//           int ngs = 0, i = l; 
//           for(auto it = column.end(); --it, --i; ) { 
//             if(!isnan2(*it)) { // !isnan
//               if(isnan2(*(fjg+g[i]-1))) { // isnan 
//                 *(fjg+g[i]-1) = *it; 
//                 ++ngs;
//                 if(ngs == ng) break;
//               }
//             }
//           }
//         }
//       } else {
//         for(int j = col; j--; ) { 
//           MatrixColumn<RTYPE> column = x( _ , j); 
//           MatrixColumn<RTYPE> lastj = last( _ , j); 
//           int ngs = 0;
//           for(int i = l; i--; ) {
//             if(column[i] != Vector<RTYPE>::get_na()) { // !isnan
//               if(lastj[g[i]-1] == Vector<RTYPE>::get_na()) { // isnan
//                 lastj[g[i]-1] = column[i]; 
//                 ++ngs;
//                 if(ngs == ng) break;
//               }
//             }
//           }
//         }
//       }
//     } else {
//       for(int j = col; j--; ) {
//         MatrixColumn<RTYPE> column = x( _ , j); 
//         MatrixColumn<RTYPE> lastj = last( _ , j); 
//         LogicalVector glj(ng, true); 
//         int ngs = 0;
//         for(int i = l; i--; ) {
//           if(glj[g[i]-1]) {
//             glj[g[i]-1] = false;
//             lastj[g[i]-1] = column[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         }
//       }
//     }
//     colnames(last) = colnames(x);
//     return last;
//   }
// }
// 
// template <>
// SEXP flastmCppImpl(Matrix<CPLXSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// SEXP flastmCppImpl(Matrix<VECSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// SEXP flastmCppImpl(Matrix<RAWSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// SEXP flastmCppImpl(Matrix<EXPRSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// // [[Rcpp::export]]
// SEXP flastmCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true, bool drop = true){
//   RCPP_RETURN_MATRIX(flastmCppImpl, x, ng, g, narm, drop);
// }




// Only Numeric Version
// // [[Rcpp::export]]
// SEXP flastmCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                 bool narm = true, bool drop = true) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector last = no_init_vector(col); // Initialize faster -> Nope !!!
//     if(narm) { 
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j); 
//         int k = l-1;
//         double lastj = column[k];
//         while(std::isnan(lastj) && k!=0) lastj = column[--k];
//         last[j] = lastj;
//       }
//     } else {
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         last[j] = column[l-1];
//       }
//     }
//     if(drop) last.attr("names") = colnames(x); 
//     else {
//       last.attr("dim") = Dimension(1, col);
//       // last.attr("dimnames") = List::create(R_NilValue,colnames(x)); 
//       colnames(last) = colnames(x); 
//     }
//     return last;
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     if(narm) {
//       NumericMatrix last = no_init_matrix(ng, col);
//       std::fill(last.begin(), last.end(), NA_REAL); 
//       for(int j = col; j--; ) { 
//         NumericMatrix::Column column = x( _ , j); 
//         NumericMatrix::Column lastj = last( _ , j); 
//         int ngs = 0;
//         for(int i = l; i--; ) {
//           if(!std::isnan(column[i])) { 
//             if(std::isnan(lastj[g[i]-1])) {
//               lastj[g[i]-1] = column[i]; 
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }
//         }
//       }
//       return last;
//     } else {
//       NumericMatrix last = no_init_matrix(ng, col); 
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j); 
//         NumericMatrix::Column lastj = last( _ , j);
//         LogicalVector glj(ng, true); // Other way around ?? -> Nope, this is faster !!
//         int ngs = 0;
//         for(int i = l; i--; ) {
//           if(glj[g[i]-1]) {
//             glj[g[i]-1] = false;
//             lastj[g[i]-1] = column[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         }
//       }
//       return last;
//     }
//   }
// }
