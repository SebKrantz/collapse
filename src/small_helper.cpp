#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector groups2GRPCpp(const List& x, int lx, const IntegerVector& gs) {
  int l = x.size(); // , sum = 0;
  // for(int i = 0; i != l; ++i) sum += gs[i]; // could also use input, just length of data...
  IntegerVector out = no_init_vector(lx); // no_init_vector(sum);
  for(int j = l; j--; ) { // This can go in any direction..
    IntegerVector column = x[j]; // const ??
    int jp = j+1;
    for(int i = gs[j]; i--; ) out[column[i]-1] = jp; // This can go in any direction...
  }
  return out;
}

// [[Rcpp::export]] // not faster than base method on large data (PRIO)
List lassignCpp(const List& x, int s, const SEXP& rows = R_NilValue, double fill = NA_REAL) {
  int l = x.size(), tr = TYPEOF(rows);
  List out(l);
  if(tr == INTSXP) {
    IntegerVector rowsv = rows;
    int rs = rowsv.size();
    for(int j = l; j--; ) {
      NumericVector column = x[j];
      if(column.size() != rs) stop("length(rows) must match nrow(x)");
      NumericVector outj(s, fill);
      for(int i = 0; i != rs; ++i) outj[rowsv[i]-1] = column[i];
      SHALLOW_DUPLICATE_ATTRIB(outj, column);
      out[j] = outj;
    }
  } else if(tr == LGLSXP) {
    LogicalVector rowsl = rows;
    int rs = rowsl.size();
    if(s != rs) stop("length(rows) must match length(s) if rows is a logical vector");
    for(int j = l; j--; ) {
      NumericVector column = x[j];
      NumericVector outj = no_init_vector(s);
      int k = 0;
      for(int i = 0; i != rs; ++i) {
        outj[i] = rowsl[i] ? column[k++] : fill;
      }
      SHALLOW_DUPLICATE_ATTRIB(outj, column);
      out[j] = outj;
    }
  } else stop("rows must be positive integers or a logical vector");
  DUPLICATE_ATTRIB(out, x);
  return out;
}


// // [[Rcpp::export]]
// SEXP setAttributes(SEXP x, SEXP a) {
//   SET_ATTRIB(x, Rf_coerceVector(a, LISTSXP));
//   Rf_classgets(x, Rf_getAttrib(x, R_ClassSymbol));
//   return x;
// }
//
// // [[Rcpp::export]]
// void setattributes(SEXP x, SEXP a) {
//   SET_ATTRIB(x, Rf_coerceVector(a, LISTSXP));
//   // SET_OBJECT(x, TYPEOF(x)); // if(OBJECT(a))  // This does not work with ts-matrices! could also make compatible with S4 objects !
//   Rf_classgets(x, Rf_getAttrib(x, R_ClassSymbol));
// }
//
// // [[Rcpp::export]]
// SEXP setAttr(SEXP x, SEXP a, SEXP v) {
//   Rf_setAttrib(x, a, v);
//   return x;
// }
//
// // [[Rcpp::export]]
// void setattr(SEXP x, SEXP a, SEXP v) {
//   Rf_setAttrib(x, a, v);
// }
//
// // [[Rcpp::export]]
// SEXP duplAttributes(SEXP x, SEXP y) { // also look at data.table's keepattributes ...
//   DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
//   return x;
// }
//
// // [[Rcpp::export]]
// void duplattributes(SEXP x, SEXP y) {
//   DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
// }
//
// // [[Rcpp::export]]
// SEXP cond_duplAttributes(SEXP x, SEXP y) {
//   if(TYPEOF(x) == TYPEOF(y)) DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
//   return x;
// }
//
// // [[Rcpp::export]]
// void cond_duplattributes(SEXP x, SEXP y) {
//   if(TYPEOF(x) == TYPEOF(y)) DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
// }




// Old / Experimental:
// // [[Rcpp::export]] // https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/main/attrib.c
// SEXP setAttributes(SEXP x, SEXP a) { // https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/main/duplicate.c
//   // bool S4l = IS_S4_OBJECT(x);
//   // int obj = OBJECT(x);
//   SET_ATTRIB(x, Rf_coerceVector(a, LISTSXP)); // Rf_shallow_duplicate(
//   // SET_OBJECT(x, obj); // or SET_OBJECT(x, OBJECT(x)) -> just OBJECT does not work !!   ?? see keepattr function in data.table_assign.c // if(OBJECT(a))
//   // SET_OBJECT(x, TYPEOF(x)); // This does not work with ts-matrices!!
//   // (S4l) ?  SET_S4_OBJECT(x) : UNSET_S4_OBJECT(x);
//   // Rf_setAttrib(x, R_ClassSymbol, Rf_getAttrib(x, R_ClassSymbol));
//   Rf_classgets(x, Rf_getAttrib(x, R_ClassSymbol)); // fast on large data (PRIO!! also works for unclassed objects)
//   return x; // wrap(x) // wrap better ?? -> error setting attributes of matrix !!-> Nope !! Still error for dapply(NGGDC, log, return = "matrix")
// }

// // [[Rcpp::export]]
// NumericVector narmCpp(NumericVector x) {
//   NumericVector y = no_init_vector(x.size());
//   std::remove_copy_if(x.begin(), x.end(), y.begin(), isnan2);
//   return y;
// }

// // [[Rcpp::export]]
// bool fanyNAint(IntegerVector x) {
//   for(int i = x.size(); i--; ) if(x[i] == NA_INTEGER) return true;
//   return false;
// }

// // [[Rcpp::export]] // not faster than base method on large data (PRIO)!!
// List NUMlassignCpp(const List& x, const IntegerVector& rows = 0, const IntegerVector& cols = 0, double fill = NA_REAL) {
//   int rs = rows.size(), cs = cols.size(), l = x.size();
//   int outs = (cs == 1 && cols[0] == 0) ? l : cs;
//   List out(outs);
//   if(outs == l && rows[0] != 0) {
//     for(int j = l; j--; ) {
//       NumericVector column = x[j];
//       int row = column.size();
//       NumericVector outj(row, fill);
//       // outj[rows] = column[rows]; // best ?? -> nope, below is better !!
//       for(int i = 0; i != rs; ++i) {
//         int ri = rows[i]-1;
//         outj[ri] = column[ri];
//       }
//       SHALLOW_DUPLICATE_ATTRIB(outj, column);
//       out[j] = outj;
//     }
//     DUPLICATE_ATTRIB(out, x);
//     return out;
//   } else if(outs != l && rows[0] != 0) {
//     for(int j = cs; j--; ) {
//       int col = cols[j]-1;
//       NumericVector column = x[col];
//       int row = column.size();
//       NumericVector outj(row, fill);
//       // outj[rows] = column[rows]; // best ?? -> nope, below is better !!
//       for(int i = 0; i != rs; ++i) {
//         int ri = rows[i]-1;
//         outj[ri] = column[ri];
//       }
//       SHALLOW_DUPLICATE_ATTRIB(outj, column);
//       out[j] = outj;
//     }
//   } else if(rs == 1 && rows[0] == 0) {
//     out = x[cols-1];
//   } else stop("Must subset either rows or columns or both!");
//   DUPLICATE_ATTRIB(out, x);
//   SEXP nam = Rf_getAttrib(x, R_NamesSymbol);
//   if(nam != R_NilValue) {
//     CharacterVector names = nam;
//     out.attr("names") = names[cols-1];
//   }
//   return out;
// }



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
