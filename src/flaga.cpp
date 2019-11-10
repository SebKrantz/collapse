// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

inline SEXP coln_check(SEXP x) { 
    if(Rf_isNull(x)) return NA_STRING;     
    else return x; // Rf_coerceVector(x, STRSXP);             
} 

template <int RTYPE>
Matrix<RTYPE> flagleadmCppImpl(const Matrix<RTYPE>& x, const IntegerVector& n, const SEXP& fill,
                               int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  
  Vector<RTYPE> fil(1);
  if(Rf_isNull(fill)) { //  || fill != fill not necessary !!
    fil = Vector<RTYPE>::get_na();
  } else {
    fil = as<Vector<RTYPE> >(fill); 
  }
  auto ff = fil[0];
  
  int l = x.nrow(), col = x.ncol(), ns = n.size(), pos = 0;
  IntegerVector absn = no_init_vector(ns);
  for(int i = 0; i != ns; ++i) {
    if(n[i]<0) absn[i] = -n[i];
    else absn[i] = n[i];
  }
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING; // NumericVector(abs(n))
  CharacterVector colnam = names ? no_init_vector(col*ns) : no_init_vector(1); // what if no names ?? 
  CharacterVector coln = names ? coln_check(colnames(x)) : NA_STRING;
  if(names && coln[0] == NA_STRING) names = false;

  Matrix<RTYPE> out = no_init_matrix(l, col*ns);

  if(ng == 0) { // No groups 
    if(Rf_isNull(t)) { // Ordered data
      for(int j = 0; j != col; ++j) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p] > l) stop("lag-length exceeds length of vector");
          MatrixColumn<RTYPE> outj = out( _ , pos);
          if(np>0) {
            if(names) colnam[pos] = "L" + nc[p] + "." + coln[j];
            int i = 0;
            while(i != np) outj[i++] = ff;
            for( ; i != l; ++i) outj[i] = column[i - np];
          } else if(np<0) {
            if(names) colnam[pos] = "F" + nc[p] + "." + coln[j];
            int i = l, st = l+np;
            while(i != st) outj[--i] = ff;
            for( ; i--; ) outj[i] = column[i - np];
          } else {
            if(names) colnam[pos] = coln[j];
            outj = column; 
          }
          ++pos;
        }
      }
    } else { // Unordered data: Timevar Provided
      IntegerVector ord = t;
      if(l != ord.size()) stop("length(x) must match length(t)");
      LogicalVector ocheck(l, true); 
      IntegerVector omap = no_init_vector(l); 
      for(int i = 0; i != l; ++i) { 
        if(ord[i] > l) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i; // Note: omap is the same as order(ord) !!
        } else {
          stop("Repeated values in timevar");
        }
      }
      for(int j = 0; j != col; ++j) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p] > l) stop("lag-length exceeds length of vector");
          MatrixColumn<RTYPE> outj = out( _ , pos);
          if(np>0) {
            if(names) colnam[pos] = "L" + nc[p] + "." + coln[j]; 
            int i = 0;
            while(i != np) outj[omap[i++]] = ff; 
            for( ; i != l; ++i) outj[omap[i]] = column[omap[i - np]]; 
          } else if(np<0) {
            if(names) colnam[pos] = "F" + nc[p] + "." + coln[j];
            int st = l+np, i = l; 
            while(i != st) outj[omap[--i]] = ff;
            for( ; i--; ) outj[omap[i]] = column[omap[i - np]]; 
          } else {
            if(names) colnam[pos] = coln[j];
            outj = column; 
          }
          ++pos;
        }
      }
    }
  } else { // With groups
    if(l != g.size()) stop("length(x) must match length(g)");
    int ags = l/ng, ngp = ng+1;
    if(Rf_isNull(t)) { // Ordered data
      int seen[ngp], memsize = sizeof(int)*ngp;
      for(int j = 0; j != col; ++j) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
          MatrixColumn<RTYPE> outj = out( _ , pos);
          if(np>0) {
            if(names) colnam[pos] = "L" + nc[p] + "." + coln[j]; 
            memset(seen, 0, memsize); 
            for(int i = 0; i != l; ++i) {  
              if(seen[g[i]] == np) {
                outj[i] = column[i-np];
              } else {
                outj[i] = ff;
                ++seen[g[i]];
              }
            }
          } else if(np<0) {
            if(names) colnam[pos] = "F" + nc[p] + "." + coln[j]; 
            memset(seen, 0, memsize);
            for(int i = l; i--; ) { // good?? 
              if(seen[g[i]] == np) {
                outj[i] = column[i-np];
              } else {
                outj[i] = ff;
                --seen[g[i]];
              }
            }
          } else {
            if(names) colnam[pos] = coln[j];
            outj = column; 
          }
          ++pos;
        }
      }
    } else { // Unordered data: Timevar provided
      IntegerVector ord = t;
      if(l != ord.size()) stop("length(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX); 
      IntegerVector gsv = NULL; 
      IntegerVector ord2 = no_init_vector(l); // See flag.cpp for any improvements on this code !!
      if(Rf_isNull(gs)) {
        gsv = IntegerVector(ng);
        for(int i = 0; i != l; ++i) {
          ++gsv[g[i]-1];
          if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
        }
      } else {
        gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
      }
      IntegerVector omap(l); 
      int cgs[ngp], index[l]; // See flag.cpp for any improvements on this code !!
      cgs[1] = 0;
      for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; // or get "starts from forderv"
      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]]; 
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        index[i] = cgs[g[i]]+ord2[i];
        if(omap[index[i]] == 0) omap[index[i]] = i; 
        else stop("Repeated values of timevar within one or more groups"); 
      }
      for(int j = 0; j != col; ++j) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
          MatrixColumn<RTYPE> outj = out( _ , pos);
          if(np>0) {
            if(names) colnam[pos] = "L" + nc[p] + "." + coln[j];
            for(int i = 0; i != l; ++i) {
              if(ord2[i] >= np) {
                outj[i] = column[omap[index[i]-np]];
              } else {
                outj[i] = ff;
              }
            }
          } else if(np<0) {
            if(names) colnam[pos] = "F" + nc[p] + "." + coln[j];
            for(int i = 0; i != l; ++i) { // best loop ??
              if(ord2[i] < gsv[g[i]-1]+np) {
                outj[i] = column[omap[index[i]-np]];
              } else {
                outj[i] = ff;
              }
            }
          } else {
            if(names) colnam[pos] = coln[j];
            outj = column;
          }
          ++pos;
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  if(ns != 1) out.attr("dim") = Dimension(l, col*ns);
  if(names) {
    out.attr("dimnames") = List::create(rownames(x), colnam); // colnames(out) = colnam deletes row names !!!
  } else if(ns != 1) {
    out.attr("dimnames") = R_NilValue;
  }
  return out;
}


// template <>
// Vector<CPLXSXP> flagleadCppImpl(Vector<CPLXSXP> x, IntegerVector n, SEXP fill,
//                               int ng, IntegerVector g, SEXP gs, SEXP t) {
//   stop("Not supported SEXP type!");
// }

template <>
Matrix<VECSXP> flagleadmCppImpl(const Matrix<VECSXP>& x, const IntegerVector& n, const SEXP& fill,
                                int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Matrix<RAWSXP> flagleadmCppImpl(const Matrix<RAWSXP>& x, const IntegerVector& n, const SEXP& fill,
                                int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Matrix<EXPRSXP> flagleadmCppImpl(const Matrix<EXPRSXP>& x, const IntegerVector& n, const SEXP& fill,
                                 int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP flagleadmCpp(SEXP x, IntegerVector n = 1, SEXP fill = R_NilValue, 
                 int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue, bool names = true){
  RCPP_RETURN_MATRIX(flagleadmCppImpl, x, n, fill, ng, g, gs, t, names);
}


// Second Numeric Version: Still 2D omap, not working unorderd time-series (non-panel)...
// // [[Rcpp::export]] 
// NumericMatrix flagleadmCpp(NumericMatrix x, IntegerVector n = 1, double fill = NA_REAL, 
//                            int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) { 
//   int l = x.nrow(), col = x.ncol(), ns = n.size();
//   
//   NumericMatrix out = no_init_matrix(l, col*ns);
//   if(ng == 0) { // No groups 
//     if(Rf_isNull(t)) { // Ordered data
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         for(int p = ns; p--; ) {
//           NumericMatrix::Column outj = out( _ , j*ns+p);
//           int np = n[p];
//           if(np>0) {
//             int i = 0;
//             while(i != np) outj[i++] = fill;
//             for( ; i != l; ++i) outj[i] = column[i - np];
//           } else if(np<0) {
//             int i = l, st = l+np;
//             while(i != st) outj[--i] = fill;
//             for( ; i--; ) outj[i] = column[i - np];
//           } else outj = column; 
//         }
//       }
//     } else { // Unordered data: Timevar Provided
//       IntegerVector ord = t; 
//       if(ord.size() != l) stop("nrow(x) must match length(t)");
//       LogicalVector ocheck(l, true); // check needed ???
//       for(int i = 0; i != l; ++i) {
//         if(ocheck[ord[i]-1]) ocheck[ord[i]-1] = false;
//         else stop("Repeated values in timevar");
//       }
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         for(int p = ns; p--; ) {
//           NumericMatrix::Column outj = out( _ , j*ns+p);
//           int np = n[p];
//           if(np>0) {
//             for(int i = 0; i != l; ++i) {
//               if(ord[i] > np) {
//                 outj[i] = column[ord[i]-np-1]; 
//               } else {
//                 outj[i] = fill;
//               }
//             }
//           } else if(np<0) {
//             int st = l+np;
//             for(int i = 0; i != l; ++i) { // best loop ?? (i.e. fastest ??)
//               if(ord[i] <= st) {
//                 outj[i] = column[ord[i]-np-1]; 
//               } else {
//                 outj[i] = fill;
//               }
//             }
//           } else outj = column; 
//         }
//       }
//     }
//   } else { // With groups
//     if(g.size() != l) stop("nrow(x) must match length(g)");
//     if(Rf_isNull(t)) { // Ordered data
//       int seen[ng], memsize = sizeof(int)*ng;
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         for(int p = ns; p--; ) {
//           NumericMatrix::Column outj = out( _ , j*ns+p);
//           int np = n[p];
//           if(np>0) {
//             memset(seen, 0, memsize); // fastest !!
//             for(int i = 0; i != l; ++i) {  
//               if(seen[g[i]-1] == np) {
//                 outj[i] = column[i-np];
//               } else {
//                 outj[i] = fill;
//                 ++seen[g[i]-1];
//               }
//             }
//           } else if(np<0) {
//             memset(seen, 0, memsize);
//             for(int i = l; i--; ) { // good?? 
//               if(seen[g[i]-1] == np) {
//                 outj[i] = column[i-np];
//               } else {
//                 outj[i] = fill;
//                 --seen[g[i]-1];
//               }
//             }
//           } else outj = column; 
//         }
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t; 
//       if(ord.size() != l) stop("nrow(x) must match length(t)");
//       IntegerVector min(ng, INT_MAX);
//       IntegerVector gsv = no_init_vector(ng); // No real improvements here by using C++ arrays !!
//       IntegerVector ord2 = no_init_vector(l); 
//       IntegerVector g2 = no_init_vector(l);  
//       if(Rf_isNull(gs)) {
//         std::fill(gsv.begin(), gsv.end(), 0); 
//         for(int i = 0; i != l; ++i) {
//           g2[i] = g[i]-1; 
//           ++gsv[g2[i]];
//           if(ord[i] < min[g2[i]]) min[g2[i]] = ord[i]; 
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)"); 
//         for(int i = 0; i != l; ++i) {
//           g2[i] = g[i]-1; 
//           if(ord[i] < min[g2[i]]) min[g2[i]] = ord[i];
//         }
//       }
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) omap[i] = new int[gsv[i]]{}; // or () // better using vector of vectors (faster initialization) ??
//       for(int i = 0; i != l; ++i) {
//         ord2[i] = ord[i] - min[g2[i]]; 
//         if(ord2[i] >= gsv[g2[i]]) stop("Gaps in timevar within one or more groups");
//         if(omap[g2[i]][ord2[i]] == 0) omap[g2[i]][ord2[i]] = i; // fastest ??
//         else stop("Repeated values of timevar within one or more groups"); 
//       }
//       for(int j = col; j--; ) {
//         NumericMatrix::Column column = x( _ , j);
//         for(int p = ns; p--; ) {
//           NumericMatrix::Column outj = out( _ , j*ns+p);
//           int np = n[p];
//           if(np>0) {
//             for(int i = 0; i != l; ++i) {
//               if(ord2[i] >= np) {
//                 outj[i] = column[omap[g2[i]][ord2[i]-np]]; 
//               } else {
//                 outj[i] = fill;
//               }
//             }
//           } else if(np<0) {
//             for(int i = 0; i != l; ++i) { // best loop ??
//               if(ord2[i] < gsv[g2[i]]+np) { 
//                 outj[i] = column[omap[g2[i]][ord2[i]-np]]; 
//               } else {
//                 outj[i] = fill;
//               }
//             }
//           } else outj = column; 
//         }
//       }
//     }
//   }
//   return out;
// }

// PREVIOUS VERSION
// // Still do: What if fill = NULL -> Delete !!
// // Also use matrix column subsetting, or better one-pass through the matrix (treat as vector ????????)
// // [[Rcpp::export]]
// NumericMatrix flagmCpp(NumericMatrix x, int n = 1, double fill = NA_REAL, 
//                        int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, SEXP o = R_NilValue) { 
//   int l = x.nrow(), col = x.ncol();
//   NumericMatrix out = no_init_matrix(l, col);
//   
//   if(ng == 0) { // No groups !!
//     for(int j = 0; j != col; ++j) {
//       NumericMatrix::Column column = x( _ , j);
//       NumericMatrix::Column outj = out( _ , j);
//       int i = 0;
//       while(i != n) outj[i++] = fill;
//       for( ; i != l; ++i) outj[i] = column[i - n];
//     }
//   } else {
//     if(l != g.size()) stop("nrow(x) must match length(g)");
//     if (Rf_isNull(o)) { // Ordered data
//       // warning("Panel-lag computed without timevar: Assuming ordered data"); -> Do in main function !!
//       for(int j = 0; j != col; ++j) {
//         NumericMatrix::Column column = x( _ , j);
//         NumericMatrix::Column outj = out( _ , j);
//         IntegerVector seenj(ng); // Faster going row-wise here ????????????????????? or in one loop ??
//         for(int i = 0; i != l; ++i) { 
//           if(seenj[g[i]-1] == n) {
//             outj[i] = column[i-n];
//           } else {
//             outj[i] = fill;
//             ++seenj[g[i]-1];
//           }
//         }
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = o; 
//       if(l != ord.size()) stop("length(x) must match length(o)");
//       if(ng != gs.size()) stop("ng must match length(gs)");
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//       for(int i = 0; i != l; ++i) omap[g[i]-1][ord[i]-1] = i;
//       for(int j = 0; j != col; ++j) {
//         NumericMatrix::Column column = x( _ , j);
//         NumericMatrix::Column outj = out( _ , j);
//         for(int i = 0; i != l; ++i) {
//           if(ord[i] > n) {
//             outj[i] = column[omap[g[i]-1][ord[i]-n-1]]; // good ??? 
//           } else {
//             outj[i] = fill;
//           }
//         }
//       }
//     }
//   }
//   return out;
// }

// Test runs !!!!!!!!!!!!
// // [[Rcpp::export]] 
// NumericMatrix lagmCpp(NumericMatrix x, int n = 1) { 
//   int l = x.nrow(), col = x.ncol(), i = 0;
//   NumericMatrix out = no_init_matrix(l, col);
//   while(i != n) out.row(i++) = rep(NA_REAL, col);
//   for( ; i != l; ++i) out.row(i) = x.row(i - n);
//   return out;
// }
// 
// // [[Rcpp::export]] 
// NumericMatrix lagm2Cpp(NumericMatrix x, int n = 1) { 
//   int l = x.nrow(), col = x.ncol(), i = 0;
//   NumericMatrix out = no_init_matrix(l, col);
//   while(i != n) out(i++, _ ) = rep(NA_REAL, col);
//   for( ; i != l; ++i) out(i, _ ) = x(i - n, _ );
//   return out;
// }
// 
// // [[Rcpp::export]] 
// NumericMatrix lagm3Cpp(NumericMatrix x, int n = 1) { 
//   int l = x.nrow(), col = x.ncol(), i = 0;
//   NumericMatrix out = no_init_matrix(l, col);
//   while(i != n) out(i++, _ ) = rep(NA_REAL, col);
//   for( ; i != l; ++i) {
//     NumericMatrix::Row row = x.row(i - n);
//     NumericMatrix::Row outi = out.row(i);
//     outi = row;
//   }
//   return out;
// }
// 
// // [[Rcpp::export]] 
// NumericMatrix lagm4Cpp(NumericMatrix x, int n = 1) { // This is the best !!
//   int l = x.nrow(), col = x.ncol();
//   NumericMatrix out = no_init_matrix(l, col);
//   for(int j = col; j--; ) {
//     NumericMatrix::Column column = x( _ , j);
//     NumericMatrix::Column outj = out( _ , j);
//     int i = 0;
//     while(i != n) outj[i++] = NA_REAL;
//     for( ; i != l; ++i) outj[i] = column[i - n];
//   }
//   return out;
// }
// 
// // [[Rcpp::export]] 
// NumericMatrix lagm5Cpp(NumericMatrix x, int n = 1) { 
//   int l = x.nrow(), col = x.ncol();
//   NumericMatrix out = no_init_matrix(l, col);
//   for(int j = col; j--; ) {
//     NumericMatrix::Column column = x.column(j);
//     NumericMatrix::Column outj = out.column(j);
//     int i = 0;
//     while(i != n) outj[i++] = NA_REAL;
//     for( ; i != l; ++i) outj[i] = column[i - n];
//   }
//   return out;
// }
