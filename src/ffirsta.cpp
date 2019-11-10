// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// with better handling row.names !!

template <int RTYPE>
inline bool isnaNUM(typename Rcpp::traits::storage_type<RTYPE>::type x) { 
  return x != x;
}

template <int RTYPE>
inline bool isnaOTH(typename Rcpp::traits::storage_type<RTYPE>::type x) { 
  return x == Vector<RTYPE>::get_na();
}


template <int RTYPE>
SEXP ffirstmCppImpl(const Matrix<RTYPE>& x, int ng, const IntegerVector& g, bool narm, bool drop) { 
  int l = x.nrow(), col = x.ncol(); 
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;
  
  if(ng == 0) { 
    Vector<RTYPE> first = no_init_vector(col); 
    if(narm) { 
      int end = l-1;
      for(int j = col; j--; ) {
        ConstMatrixColumn<RTYPE> column = x( _ , j); 
        int k = 0;
        auto firstj = column[k];
        while(isnanT(firstj) && k!=end) firstj = column[++k];
        first[j] = firstj;
      }
    } else {
      first = x(0, _); 
    }
    if(drop) first.attr("names") = colnames(x); 
    else {
      first.attr("dim") = Dimension(1, col);
      colnames(first) = colnames(x); 
    }
    return first;
  } else { // with groups 
    if(g.size() != l) stop("length(g) must match nrow(X)");
    Matrix<RTYPE> first = no_init_matrix(ng, col);
    if(narm) {
      std::fill(first.begin(), first.end(), Vector<RTYPE>::get_na()); 
      for(int j = col; j--; ) { 
        ConstMatrixColumn<RTYPE> column = x( _ , j); 
        MatrixColumn<RTYPE> firstj = first( _ , j); 
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(!isnanT(column[i])) { 
            if(isnanT(firstj[g[i]-1])) { 
              firstj[g[i]-1] = column[i]; 
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
      }
      colnames(first) = colnames(x);
    } else {
      List dn = x.attr("dimnames");
      if(dn[0] != R_NilValue) {
        CharacterVector rn = dn[0];
        CharacterVector newrn = no_init_vector(ng);
        LogicalVector glj(ng, true); // using std::vector<bool> here is more memory efficient but not faster !!
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(glj[g[i]-1]) {
            glj[g[i]-1] = false;
            first(g[i]-1, _) = x(i, _);
            newrn[g[i]-1] = rn[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
        first.attr("dimnames") = List::create(newrn, dn[1]); // best way !!
      } else {
        LogicalVector glj(ng, true); // using std::vector<bool> here is more memory efficient but not faster !!
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(glj[g[i]-1]) {
            glj[g[i]-1] = false;
            first(g[i]-1, _) = x(i, _);
            ++ngs;
            if(ngs == ng) break;
          }
        }
        colnames(first) = colnames(x);
      }
    }
    return first;
  }
}

template <>
SEXP ffirstmCppImpl(const Matrix<CPLXSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP ffirstmCppImpl(const Matrix<VECSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP ffirstmCppImpl(const Matrix<RAWSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP ffirstmCppImpl(const Matrix<EXPRSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP ffirstmCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true, bool drop = true){
  RCPP_RETURN_MATRIX(ffirstmCppImpl, x, ng, g, narm, drop);
}


// Old code: Solving the na.rm problem using a macro and an iterator (for commented code see below)
// #define isnan2(x) ((x) != (x)) // http://www.ebyte.it/library/codesnippets/WritingCppMacros.html
// 
// template <int RTYPE>
// SEXP ffirstmCppImpl(Matrix<RTYPE> x, int ng, IntegerVector g, bool narm, bool drop) { 
//   int l = x.nrow(), col = x.ncol(); 
//   
//   
//   if(ng == 0) { 
//     Vector<RTYPE> first = no_init_vector(col); 
//     if(narm) { 
//       if(TYPEOF(x) == REALSXP) {
//         for(int j = col; j--; ) {
//           MatrixColumn<RTYPE> column = x( _ , j); 
//           int k = 0;
//           auto firstj = column[k];
//           while(firstj != firstj && k!=l-1) firstj = column[++k];
//           first[j] = firstj;
//         }
//       } else {
//         for(int j = col; j--; ) {
//           MatrixColumn<RTYPE> column = x( _ , j); 
//           int k = 0;
//           auto firstj = column[k];
//           while(firstj == Vector<RTYPE>::get_na() && k!=l-1) firstj = column[++k];
//           first[j] = firstj;
//         }
//       }
//     } else {
//       for(int j = col; j--; ) {
//         MatrixColumn<RTYPE> column = x( _ , j);
//         first[j] = column[0]; // Good?? What about character?? 
//       }
//     }
//     if(drop) first.attr("names") = colnames(x); 
//     else {
//       first.attr("dim") = Dimension(1, col);
//       colnames(first) = colnames(x); 
//     }
//     return first;
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     Matrix<RTYPE> first = no_init_matrix(ng, col);
//     if(narm) {
//       std::fill(first.begin(), first.end(), Vector<RTYPE>::get_na()); 
//       if(TYPEOF(x) == REALSXP) {
//         for(int j = col; j--; ) { 
//           MatrixColumn<RTYPE> column = x( _ , j); 
//           MatrixColumn<RTYPE> firstj = first( _ , j); 
//           auto fjg = firstj.begin();
//           int ngs = 0, i = 0; 
//           for(auto it = column.begin(); it != column.end(); ++it, ++i) { 
//             if(!isnan2(*it)) { 
//               if(isnan2(*(fjg+g[i]-1))) { 
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
//           MatrixColumn<RTYPE> firstj = first( _ , j); 
//           int ngs = 0;
//           for(int i = 0; i != l; ++i) {
//             if(column[i] != Vector<RTYPE>::get_na()) { 
//               if(firstj[g[i]-1] == Vector<RTYPE>::get_na()) { 
//                 firstj[g[i]-1] = column[i]; 
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
//         MatrixColumn<RTYPE> firstj = first( _ , j); 
//         LogicalVector glj(ng, true); // using std::vector<bool> here is more memory efficient but not faster !!
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(glj[g[i]-1]) {
//             glj[g[i]-1] = false;
//             firstj[g[i]-1] = column[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         }
//       }
//     }
//     colnames(first) = colnames(x);
//     return first;
//   }
// }
// 
// template <>
// SEXP ffirstmCppImpl(Matrix<CPLXSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// SEXP ffirstmCppImpl(Matrix<VECSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// SEXP ffirstmCppImpl(Matrix<RAWSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// SEXP ffirstmCppImpl(Matrix<EXPRSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// // [[Rcpp::export]]
// SEXP ffirstmCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true, bool drop = true){
//   RCPP_RETURN_MATRIX(ffirstmCppImpl, x, ng, g, narm, drop);
// }






// Same Old code, with lots of comments !!!
// template <int RTYPE>
// SEXP ffirstmCppImpl(Matrix<RTYPE> x, int ng, IntegerVector g, bool narm, bool drop) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     Vector<RTYPE> first = no_init_vector(col); // Initialize faster -> Nope !!!
//     if(narm) { 
//       if(TYPEOF(x) == REALSXP) {
//         for(int j = col; j--; ) {
//           MatrixColumn<RTYPE> column = x( _ , j); 
//           // Vector<RTYPE> column = x.column(j); // alternatively !! faster ?? 
//           int k = 0;
//           auto firstj = column[k];
//           while(firstj != firstj && k!=l-1) firstj = column[++k];
//           first[j] = firstj;
//         }
//       } else {
//         for(int j = col; j--; ) {
//           MatrixColumn<RTYPE> column = x( _ , j); 
//           // Vector<RTYPE> column = x.column(j); // alternatively !! faster ?? 
//           int k = 0;
//           auto firstj = column[k];
//           while(firstj == Vector<RTYPE>::get_na() && k!=l-1) firstj = column[++k];
//           first[j] = firstj;
//         }
//       }
//     } else {
//       for(int j = col; j--; ) {
//         MatrixColumn<RTYPE> column = x( _ , j);
//         // Vector<RTYPE> column = x.column(j); // alternatively !! faster ?? 
//         first[j] = column[0]; // Good?? What about character?? 
//       }
//     }
//     if(drop) first.attr("names") = colnames(x); 
//     else {
//       first.attr("dim") = Dimension(1, col);
//       // first.attr("dimnames") = List::create(R_NilValue,colnames(x)); 
//       colnames(first) = colnames(x); 
//     }
//     return first;
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     Matrix<RTYPE> first = no_init_matrix(ng, col);
//     // typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t; Not faster than auto !!
//     if(narm) {
//       std::fill(first.begin(), first.end(), Vector<RTYPE>::get_na()); 
//       if(TYPEOF(x) == REALSXP) {
//         for(int j = col; j--; ) { 
//           MatrixColumn<RTYPE> column = x( _ , j); 
//           // Vector<RTYPE> column = x.column(j); -> This is not faster!! and not identical !!
//           MatrixColumn<RTYPE> firstj = first( _ , j); 
//           // Vector<RTYPE> firstj  = first.column(j);
//           auto fjg = firstj.begin();
//           int ngs = 0;
//           int i = 0; // https://stackoverflow.com/questions/17016376/how-to-write-a-for-loop-that-uses-both-an-iterator-and-an-index-counter
//           for(auto it = column.begin(); it != column.end(); ++it, ++i) { //  int ind = std::distance(column.begin(), i); (get index from iterator) more: https://stackoverflow.com/questions/2152986/what-is-the-most-effective-way-to-get-the-index-of-an-iterator-of-an-stdvector
//             if(*it == *it) { // !isnan
//               if(*(fjg+g[i]-1) != *(fjg+g[i]-1)) { // isnan firstj[g[i]-1]
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
//           // Vector<RTYPE> column = x.column(j);
//           MatrixColumn<RTYPE> firstj = first( _ , j); 
//           // Vector<RTYPE> firstj  = first.column(j);
//           int ngs = 0;
//           for(int i = 0; i != l; ++i) {
//             if(column[i] != Vector<RTYPE>::get_na()) { // !isnan
//               if(firstj[g[i]-1] == Vector<RTYPE>::get_na()) { // isnan
//                 firstj[g[i]-1] = column[i]; 
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
//         // Vector<RTYPE> column = x.column(j);
//         MatrixColumn<RTYPE> firstj = first( _ , j); 
//         // Vector<RTYPE> firstj  = first.column(j);
//         LogicalVector glj(ng, true); // Other way around ?? -> Nope, this is faster !!
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(glj[g[i]-1]) {
//             glj[g[i]-1] = false;
//             firstj[g[i]-1] = column[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         }
//       }
//     }
//     return first;
//   }
// }
// 
// template <>
// SEXP ffirstmCppImpl(Matrix<CPLXSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// SEXP ffirstmCppImpl(Matrix<VECSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// SEXP ffirstmCppImpl(Matrix<RAWSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// SEXP ffirstmCppImpl(Matrix<EXPRSXP> x, int ng, IntegerVector g, bool narm, bool drop) {
//   stop("Not supported SEXP type!");
// }
// 
// // [[Rcpp::export]]
// SEXP ffirstmCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true, bool drop = true){
//   RCPP_RETURN_MATRIX(ffirstmCppImpl, x, ng, g, narm, drop);
// }


// Previous Version: Only for numeric variables !!
// // [[Rcpp::export]]
// SEXP ffirstmCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0,  
//                bool narm = true, bool drop = true) { 
//   int l = x.nrow(); 
//   int col = x.ncol(); 
//   
//   if(ng == 0) { 
//     NumericVector first = no_init_vector(col); // Initialize faster -> Nope !!!
//     if(narm) { 
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j); 
//         int k = 0;
//         double firstj = column[k];
//         while(std::isnan(firstj) && k!=l-1) firstj = column[++k];
//         first[j] = firstj;
//       }
//     } else {
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         first[j] = column[0];
//       }
//     }
//     if(drop) first.attr("names") = colnames(x); 
//     else {
//       first.attr("dim") = Dimension(1, col);
//       // first.attr("dimnames") = List::create(R_NilValue,colnames(x)); 
//       colnames(first) = colnames(x); 
//     }
//     return first;
//   } else { // with groups 
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     if(narm) {
//       NumericMatrix first = no_init_matrix(ng, col);
//       std::fill(first.begin(), first.end(), NA_REAL); 
//       for(int j = col; j--; ) { 
//         NumericMatrix::Column column = x( _ , j); 
//         NumericMatrix::Column firstj = first( _ , j); 
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(!std::isnan(column[i])) { 
//             if(std::isnan(firstj[g[i]-1])) {
//               firstj[g[i]-1] = column[i]; 
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }
//         }
//       }
//       return first;
//     } else {
//       NumericMatrix first = no_init_matrix(ng, col); 
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j); 
//         NumericMatrix::Column firstj = first( _ , j);
//         LogicalVector glj(ng, true); // Other way around ?? -> Nope, this is faster !!
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(glj[g[i]-1]) {
//             glj[g[i]-1] = false;
//             firstj[g[i]-1] = column[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         }
//       }
//       return first;
//     }
//   }
// }
