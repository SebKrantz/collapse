#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP setAttributes(SEXP x, SEXP a) {
  SET_ATTRIB(x, Rf_coerceVector(a, LISTSXP));
  SET_OBJECT(x, TYPEOF(x)); // or SET_OBJECT(x, OBJECT(x)) ?? see keepattr function in data.table_assign.c // if(OBJECT(a))
  return x; // wrap(x) // wrap better ?? -> error setting attributes of matrix !!-> Nope !! Still error for dapply(NGGDC, log, return = "matrix")
}

// [[Rcpp::export]]
void setattributes(SEXP x, SEXP a) {
  SET_ATTRIB(x, Rf_coerceVector(a, LISTSXP));
  SET_OBJECT(x, TYPEOF(x)); // if(OBJECT(a))
}

// [[Rcpp::export]]
SEXP setAttr(SEXP x, SEXP a, SEXP v) {
  Rf_setAttrib(x, a, v);
  return x;
}

// [[Rcpp::export]]
void setattr_clp(SEXP x, SEXP a, SEXP v) {
  Rf_setAttrib(x, a, v);
}

// [[Rcpp::export]]
SEXP duplAttributes(SEXP x, SEXP y) {
  DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
  return x;
}

// [[Rcpp::export]]
void duplattributes(SEXP x, SEXP y) {
  DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
}

// [[Rcpp::export]]
SEXP cond_duplAttributes(SEXP x, SEXP y) {
  if(TYPEOF(x) == TYPEOF(y)) DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
  return x;
}

// [[Rcpp::export]]
void cond_duplattributes(SEXP x, SEXP y) {
  if(TYPEOF(x) == TYPEOF(y)) DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
}

// // [[Rcpp::export]] // not faster than base method on large data (PRIO)!!
// List LSubset(const List& x, const IntegerVector& rows = 0, const IntegerVector& cols = 0) {
//   int rs = rows.size(), cs = cols.size(), l = x.size();
//   int outs = (cs == 1 && cols[0] == 0) ? l : cs;
//   List out(outs);
//   if(outs == l && rows[0] != 0) {
//     for(int j = l; j--; ) {
//       switch(TYPEOF(x[j])) {
//       case REALSXP: {
//         NumericVector column = x[j];
//         out[j] = column[rows];
//         break;
//       }
//       case INTSXP: {
//         IntegerVector column = x[j];
//         out[j] = column[rows];
//         break;
//       }
//       case STRSXP: {
//         CharacterVector column = x[j];
//         out[j] = column[rows];
//         break;
//       }
//       case LGLSXP: {
//         LogicalVector column = x[j];
//         out[j] = column[rows];
//         break;
//       }
//       case CPLXSXP: {
//         ComplexVector column = x[j];
//         out[j] = column[rows];
//         break;
//       }
//       case VECSXP: {
//         List column = x[j];
//         out[j] = column[rows];
//         break;
//       }
//       default: stop("Not supported SEXP type!");
//       }
//       SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
//     }
//     // DUPLICATE_ATTRIB(out, x);
//     // if(Rf_isFrame(x)) out.attr("row.names") = wrap(x.attr("row.names"))[rows];
//   } else if(outs != l && rows[0] != 0) {
//     for(int j = cs; j--; ) {
//       int col = cols[j]-1;
//       switch(TYPEOF(x[col])) {
//       case REALSXP: {
//         NumericVector column = x[col];
//         out[j] = column[rows];
//         break;
//       }
//       case INTSXP: {
//         IntegerVector column = x[col];
//         out[j] = column[rows];
//         break;
//       }
//       case STRSXP: {
//         CharacterVector column = x[col];
//         out[j] = column[rows];
//         break;
//       }
//       case LGLSXP: {
//         LogicalVector column = x[col];
//         out[j] = column[rows];
//         break;
//       }
//       case CPLXSXP: {
//         ComplexVector column = x[col];
//         out[j] = column[rows];
//         break;
//       }
//       case VECSXP: {
//         List column = x[col];
//         out[j] = column[rows];
//         break;
//       }
//       default: stop("Not supported SEXP type!");
//       }
//       SHALLOW_DUPLICATE_ATTRIB(out[j], x[col]);
//     }
//   } else if(rs == 1 && rows[0] == 0) {
//    out = x[cols-1];
//   } else stop("Must subset either rows or columns or both!");
//   return out;
// }


// // [[Rcpp::export]] // needed?? for what ??
// void sduplattributes(SEXP x, SEXP y) {
//   SHALLOW_DUPLICATE_ATTRIB(x, y); // DUPLICATE_ATTRIB(x, y);
// }

// // [[Rcpp::export]] // Not necessary !!
// SEXP getattributes(SEXP x) {
//    return ATTRIB(x); // DUPLICATE_ATTRIB(x, y);
// }





// // [[Rcpp::export]]
// CharacterVector lagnamesCpp(CharacterVector x, IntegerVector n) {
//  int xs = x.size(), ns = n.size();
//   CharacterVector out = no_init_vector(xs*ns);
//   for(int j = xs; j--; ) {
//     for(int i = ns; i--; ) {
//       if(n[i]>0) out[ns*j+i] = collapse(CharacterVector::create("L",std::to_string(n[i]),".",x[j])); // "L" + std::to_string(n[i]) + "." + x[j]; // or collapse(CharacterVector::Create(...))
//       else if(n[i] == 0) out[ns*j+i] = x[j];
//       else out[ns*j+i] = collapse(CharacterVector::create("F",std::to_string(abs(n[i])),".",x[j])); // "F" + std::to_string(abs(n[i])) + "." + x[j];
//     }
//   }
//   return out;
// }
