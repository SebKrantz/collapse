#include <Rcpp/Lightest>
using namespace Rcpp;

template <int RTYPE>
List mrtlImpl(Matrix<RTYPE> X, bool names, int ret) {
  int l = X.nrow();
  List out(l);

  for(int i = l; i--; ) out[i] = X(i, _);

  if(names) {
    SEXP dn = Rf_getAttrib(X, R_DimNamesSymbol);
    if(dn == R_NilValue) dn = List::create(R_NilValue, R_NilValue); // should also work for plain matrices !
      if(Rf_isNull(VECTOR_ELT(dn, 0))) {
        CharacterVector rn(l);
        std::string VS = std::string("V"); // faster !
        for (int i = l; i--; ) rn[i] = VS + std::to_string(i+1);
        Rf_namesgets(out, rn);
      } else Rf_namesgets(out, VECTOR_ELT(dn, 0));
      if(ret != 0) {
        if(Rf_isNull(VECTOR_ELT(dn, 1)) || ret == 2) {
          Rf_setAttrib(out, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -X.ncol()));
        } else Rf_setAttrib(out, R_RowNamesSymbol, VECTOR_ELT(dn, 1));
        if(ret == 1) {
          Rf_classgets(out, Rf_mkString("data.frame"));
        } else {
          Rf_classgets(out, CharacterVector::create("data.table","data.frame"));
        }
      }
  } else if (ret != 0) {
    CharacterVector rn(l);
    std::string VS = std::string("V"); // faster !
    for (int i = l; i--; ) rn[i] = VS + std::to_string(i+1);
    Rf_namesgets(out, rn);
    Rf_setAttrib(out, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -X.ncol()));
    if(ret == 1) {
      Rf_classgets(out, Rf_mkString("data.frame"));
    } else {
      Rf_classgets(out, CharacterVector::create("data.table","data.frame"));
    }
  }
  return out;
}


template <>
List mrtlImpl(Matrix<RAWSXP> X, bool names, int ret) {
  stop("Not supported SEXP type!");
}

template <>
List mrtlImpl(Matrix<EXPRSXP> X, bool names, int ret) {
  stop("Not supported SEXP type!");
}


// [[Rcpp::export]]
SEXP mrtl(const SEXP& X, bool names = false, int ret = 0){
  RCPP_RETURN_MATRIX(mrtlImpl, X, names, ret);
}


template <int RTYPE>
List mctlImpl(Matrix<RTYPE> X, bool names, int ret) {
    int l = X.ncol();
    List out(l);

    for(int i = l; i--; ) out[i] = X(_, i);

    if(names) {
      SEXP dn = Rf_getAttrib(X, R_DimNamesSymbol);
      if(dn == R_NilValue) dn = List::create(R_NilValue, R_NilValue); // should also work for plain matrices !
        if(Rf_isNull(VECTOR_ELT(dn, 1))) {
          CharacterVector rn(l);
          std::string VS = std::string("V"); // faster !
          for (int i = l; i--; ) rn[i] = VS + std::to_string(i+1);
          Rf_namesgets(out, rn);
        } else Rf_namesgets(out, VECTOR_ELT(dn, 1));
        if(ret != 0) {
          if(Rf_isNull(VECTOR_ELT(dn, 0)) || ret == 2) {
            Rf_setAttrib(out, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -X.nrow()));
          } else Rf_setAttrib(out, R_RowNamesSymbol, VECTOR_ELT(dn, 0));
          if(ret == 1) {
            Rf_classgets(out, Rf_mkString("data.frame"));
          } else {
            Rf_classgets(out, CharacterVector::create("data.table","data.frame"));
          }
        }
    } else if (ret != 0) {
      CharacterVector rn(l);
      std::string VS = std::string("V"); // faster !
      for (int i = l; i--; ) rn[i] = VS + std::to_string(i+1);
      Rf_namesgets(out, rn);
      Rf_setAttrib(out, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -X.nrow()));
      if(ret == 1) {
        Rf_classgets(out, Rf_mkString("data.frame"));
      } else {
        Rf_classgets(out, CharacterVector::create("data.table","data.frame"));
      }
    }
    return out;
}



template <>
List mctlImpl(Matrix<RAWSXP> X, bool names, int ret) {
  stop("Not supported SEXP type!");
}

template <>
List mctlImpl(Matrix<EXPRSXP> X, bool names, int ret) {
  stop("Not supported SEXP type!");
}



// [[Rcpp::export]]
SEXP mctl(const SEXP& X, bool names = false, int ret = 0){
  RCPP_RETURN_MATRIX(mctlImpl, X, names, ret);
}



// Experimental Matrix apply functions -> Need to make faster, see Hmisc::mApply
// template <int RTYPE> // Slower than lapply(mctl...)
//  List mrtlapplyImpl(Matrix<RTYPE> X, Function FUN, bool names, int ret) {
//   int l = X.nrow();
//   List out(l);
//   for(int i = l; i--; ) {
//     MatrixRow<RTYPE> Xi = X(i,_);
//     out[i] = FUN(Xi);
//   }
//   if(names && X.hasAttribute("dimnames")) {
//     List dn(2);
//     dn = X.attr("dimnames");
//     if (Rf_isNull(dn[0])) {
//       CharacterVector rn(l);
//       for (int i = l; i--; ) {
//         rn[i] = std::string("V") + std::to_string(i+1);
//       }
//       out.attr("names") = rn;
//     } else out.attr("names") = dn[0];
//     if (ret != 0) {
//       if (Rf_isNull(dn[1])) {
//         out.attr("row.names") = NumericVector::create(NA_REAL,-X.ncol());
//       } else out.attr("row.names") = dn[1];
//       if(ret == 1) {
//         out.attr("class") = "data.frame";
//       } else {
//         out.attr("class") = CharacterVector::create("data.table","data.frame");
//       }
//     }
//   } else if (ret != 0) {
//     CharacterVector rn(l);
//     for (int i = l; i--; ) {
//       rn[i] = std::string("V") + std::to_string(i+1);
//     }
//     out.attr("names") = rn;
//     out.attr("row.names") = NumericVector::create(NA_REAL,-X.ncol());
//     if (ret == 1) {
//       out.attr("class") = "data.frame";
//     } else {
//       out.attr("class") = CharacterVector::create("data.table","data.frame");
//     }
//   }
//   return out;
// }

// template <int RTYPE>
//  Matrix<RTYPE> mrtmapplyImpl(Matrix<RTYPE> X, Function FUN) {
//   int l = X.nrow();
//   Vector<RTYPE> out0 = FUN(X(0,_)); // What if not same type ??
//   int col = out0.size();
//   Matrix<RTYPE> out = no_init_matrix(l, col);
//   for(int i = 1; i != l; ++i) {
//     out(i,_) = FUN(X(i,_));
//   }
//   if(X.ncol() == col) SHALLOW_DUPLICATE_ATTRIB(out, X);
//   else rownames(out) = rownames(X);
//   return out;
// }

// template <int RTYPE> // Slower than lapply(mctl...)
//  List mctlapplyImpl(Matrix<RTYPE> X, Function FUN, bool names, int ret) {
//   int l = X.ncol();
//   List out(l);
//   for(int i = l; i--; ) {
//     MatrixColumn<RTYPE> Xi = X(_,i);
//     out[i] = FUN(Xi);
//   }
//   if(names && X.hasAttribute("dimnames")) {
//     List dn(2);
//     dn = X.attr("dimnames");
//     if (Rf_isNull(dn[1])) {
//       CharacterVector cn(l);
//       for (int i = l; i--; ) {
//         cn[i] = std::string("V") + std::to_string(i+1);
//       }
//       out.attr("names") = cn;
//     } else out.attr("names") = dn[1];
//     if (ret != 0) {
//       if (Rf_isNull(dn[0])) {
//         out.attr("row.names") = NumericVector::create(NA_REAL,-X.nrow());
//       } else out.attr("row.names") = dn[0];
//       if(ret == 1) {
//         out.attr("class") =  "data.frame";
//       } else {
//         out.attr("class") = CharacterVector::create("data.table","data.frame");
//       }
//     }
//   } else if (ret != 0) {
//     CharacterVector cn(l);
//     for (int i = l; i--; ) {
//       cn[i] = std::string("V") + std::to_string(i+1);
//     }
//     out.attr("names") = cn;
//     out.attr("row.names") = NumericVector::create(NA_REAL,-X.nrow());
//     if (ret == 1) {
//       out.attr("class") = "data.frame";
//     } else {
//       out.attr("class") = CharacterVector::create("data.table","data.frame");
//     }
//   }
//   return out;
// }

// template <int RTYPE>
//  Matrix<RTYPE> mctmapplyImpl(Matrix<RTYPE> X, Function FUN) {
//   int l = X.ncol();
//   Vector<RTYPE> out0 = FUN(X(_,0)); // What if not same type ??
//   int row = out0.size();
//   Matrix<RTYPE> out = no_init_matrix(row, l);
//   for(int i = 1; i != l; ++i) {
//     NumericMatrix::Column outi = out(_,i);
//      outi = FUN(X(_,i));
//   }
//   if(X.nrow() == row) SHALLOW_DUPLICATE_ATTRIB(out, X);
//   else colnames(out) = colnames(X);
//   return out;
// }


// // [[Rcpp::export]]
// SEXP mrtlapply(SEXP X, Function FUN, bool names = false, int ret = 0){
//   RCPP_RETURN_MATRIX(mrtlapplyImpl, X, FUN, names, ret);
// }

// // [[Rcpp::export]]
// SEXP mrtmapply(SEXP X, Function FUN){
//   RCPP_RETURN_MATRIX(mrtmapplyImpl, X, FUN);
// }


// // [[Rcpp::export]]
// SEXP mctlapply(SEXP X, Function FUN, bool names = false, int ret = 0){
//   RCPP_RETURN_MATRIX(mctlapplyImpl, X, FUN, names, ret);
// }

// // [[Rcpp::export]]
// SEXP mctmapply(SEXP X, Function FUN){
//   RCPP_RETURN_MATRIX(mctmapplyImpl, X, FUN);
// }
