// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>
IntegerMatrix pwNobsmCppImpl(const Matrix<RTYPE>& x) {
  int l = x.nrow(), col = x.ncol();
  auto isnnanT = (RTYPE == REALSXP) ? [](typename Rcpp::traits::storage_type<RTYPE>::type x) { return x == x; } :
    [](typename Rcpp::traits::storage_type<RTYPE>::type x) { return x != Vector<RTYPE>::get_na(); };
    IntegerMatrix out = no_init_matrix(col, col);
    for(int j = 0; j != col; ++j) {
      ConstMatrixColumn<RTYPE> colj = x( _ , j);
      int nj = std::count_if(colj.begin(), colj.end(), isnnanT);
      out(j, j) = nj;
      for(int k = j+1; k != col; ++k) {
        ConstMatrixColumn<RTYPE> colk = x( _ , k);
        int njk = 0;
        for(int i = l; i--; ) if(isnnanT(colj[i]) && isnnanT(colk[i])) ++njk; // fastest? or save logical vector with colj Non-missing?
        out(j, k) = out(k, j) = njk;
      }
    }
    out.attr("dimnames") = List::create(colnames(x), colnames(x));
    return out;
}

template <>
IntegerMatrix pwNobsmCppImpl(const Matrix<CPLXSXP>& x) {
  stop("Not supported SEXP type!");
}

template <>
IntegerMatrix pwNobsmCppImpl(const Matrix<VECSXP>& x) {
  stop("Not supported SEXP type!");
}

template <>
IntegerMatrix pwNobsmCppImpl(const Matrix<RAWSXP>& x) {
  stop("Not supported SEXP type!");
}

template <>
IntegerMatrix pwNobsmCppImpl(const Matrix<EXPRSXP>& x) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
IntegerMatrix pwNobsmCpp(SEXP x){
  RCPP_RETURN_MATRIX(pwNobsmCppImpl, x);
}



// Old / Experimental:
//
// inline bool nisnan(double x) {
//   return x == x;
// }
//
// Not fast !!! :
// // [[Rcpp::export]]
// IntegerMatrix pwNobslCpp(const List& x) {
//   int l = x.size();
//   IntegerMatrix out = no_init_matrix(l, l);
//   for(int j = 0; j != l; ++j) {
//     switch(TYPEOF(x[j])) {
//     case REALSXP: {
//       NumericVector colj = x[j];
//       int nj = std::count_if(colj.begin(), colj.end(), nisnan);
//       int rowj = colj.size();
//       out(j, j) = nj;
//       for(int k = j+1; k != l; ++k) {
//         switch(TYPEOF(x[k])) {
//         case REALSXP: {
//           NumericVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(nisnan(colj[i]) && nisnan(colk[i])) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case INTSXP: {
//           IntegerVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(nisnan(colj[i]) && colk[i] != NA_INTEGER) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case STRSXP: {
//           CharacterVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(nisnan(colj[i]) && colk[i] != NA_STRING) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case LGLSXP: {
//           CharacterVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(nisnan(colj[i]) && colk[i] != NA_LOGICAL) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         default: stop("incompatible SEXP encountered;");
//         }
//       }
//       break;
//     }
//     case INTSXP: {
//       IntegerVector colj = x[j];
//       int rowj = colj.size();
//       int nj = rowj - std::count(colj.begin(), colj.end(), NA_INTEGER);
//       out(j, j) = nj;
//       for(int k = j+1; k != l; ++k) {
//         switch(TYPEOF(x[k])) {
//         case REALSXP: {
//           NumericVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_INTEGER && nisnan(colk[i])) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case INTSXP: {
//           IntegerVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_INTEGER && colk[i] != NA_INTEGER) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case STRSXP: {
//           CharacterVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_INTEGER && colk[i] != NA_STRING) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case LGLSXP: {
//           CharacterVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_INTEGER && colk[i] != NA_LOGICAL) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         default: stop("incompatible SEXP encountered;");
//         }
//       }
//       break;
//     }
//     case STRSXP: {
//       CharacterVector colj = x[j];
//       int rowj = colj.size();
//       int nj = rowj - std::count(colj.begin(), colj.end(), NA_STRING);
//       out(j, j) = nj;
//       for(int k = j+1; k != l; ++k) {
//         switch(TYPEOF(x[k])) {
//         case REALSXP: {
//           NumericVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_STRING && nisnan(colk[i])) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case INTSXP: {
//           IntegerVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_STRING && colk[i] != NA_INTEGER) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case STRSXP: {
//           CharacterVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_STRING && colk[i] != NA_STRING) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case LGLSXP: {
//           CharacterVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_STRING && colk[i] != NA_LOGICAL) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         default: stop("incompatible SEXP encountered;");
//         }
//       }
//       break;
//     }
//     case LGLSXP: {
//       LogicalVector colj = x[j];
//       int rowj = colj.size();
//       int nj = rowj - std::count(colj.begin(), colj.end(), NA_LOGICAL);
//       out(j, j) = nj;
//       for(int k = j+1; k != l; ++k) {
//         switch(TYPEOF(x[k])) {
//         case REALSXP: {
//           NumericVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_LOGICAL && nisnan(colk[i])) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case INTSXP: {
//           IntegerVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_LOGICAL && colk[i] != NA_INTEGER) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case STRSXP: {
//           CharacterVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_LOGICAL && colk[i] != NA_STRING) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         case LGLSXP: {
//           CharacterVector colk = x[k];
//           if(colk.size() != rowj) stop("All columns need to have the same length!");
//           int njk = 0;
//           for(int i = rowj; i--; ) if(colj[i] != NA_LOGICAL && colk[i] != NA_LOGICAL) ++njk; // fastest? or save logical vector with colj Non-missing?
//           out(j, k) = out(k, j) = njk;
//           break;
//         }
//         default: stop("incompatible SEXP encountered;");
//         }
//       }
//       break;
//     }
//     default:
//       stop("incompatible SEXP encountered;");
//     }
//   }
//   out.attr("dimnames") = List::create(x.attr("names"), x.attr("names"));
//   return out;
// }



//
// // [[Rcpp::export]]
// IntegerMatrix pwNobslCpp(const List& x) {
//   int l = x.size();
//   IntegerMatrix out = no_init_matrix(l, l);
//   for(int j = 0; j != l; ++j) {
//     int RTYPEj = TYPEOF(x[j]);
//     auto isnnanTj = (RTYPEj == REALSXP) ? [](typename Rcpp::traits::storage_type<RTYPEj>::type x) { return x == x; } :
//       [](typename Rcpp::traits::storage_type<RTYPEj>::type x) { return x != Vector<RTYPEj>::get_na(); };
//     Vector<RTYPEj> colj = x[j];
//     int nj = std::count_if(colj.begin(), colj.end(), isnnanTj);
//     int rowj = colj.size();
//     out(j, j) = nj;
//     for(int k = j+1; k != col; ++k) {
//       int RTYPEk = TYPEOF(x[k]);
//       auto isnnanTk = (RTYPEk == REALSXP) ? [](typename Rcpp::traits::storage_type<RTYPEk>::type x) { return x == x; } :
//         [](typename Rcpp::traits::storage_type<RTYPEk>::type x) { return x != Vector<RTYPEk>::get_na(); };
//       Vector<RTYPEk> colk = x[k];
//       if(colk.size() != rowj) stop("All columns need to have the same length!");
//       int njk = 0;
//       for(int i = rowj; i--; ) if(isnnanTj(colj[i]) && isnnanTk(colk[i])) ++njk; // fastest? or save logical vector with colj Non-missing?
//       out(j, k) = out(k, j) = njk;
//     }
//   }
//   out.attr("dimnames") = List::create(names(x), names(x));
//   return out;
// }
