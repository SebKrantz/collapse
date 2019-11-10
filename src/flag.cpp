// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Still do: What if fill = NULL -> Delete !!
// Make compatible with any kind of data ??

// 6th version: Type Dispatch and names argument !!
template <int RTYPE>
Vector<RTYPE> flagleadCppImpl(const Vector<RTYPE>& x, const IntegerVector& n, const SEXP& fill,
                              int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  
  // typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
  // storage_t fil;
  Vector<RTYPE> fil(1);
  if(Rf_isNull(fill)) { // fill != fill // Not necessary !!
    fil = Vector<RTYPE>::get_na();
  } else {
    fil = as<Vector<RTYPE> >(fill); //as<storage_t>(fill); -> doesn't work for Character vector fill !!
  }
  auto ff = fil[0];

  int l = x.size(), ns = n.size();
  IntegerVector absn = no_init_vector(ns);
  for(int i = 0; i != ns; ++i) {
    if(n[i]<0) absn[i] = -n[i];
    else absn[i] = n[i];
  }
  if(ns == 1) names = false;
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING;  // NumericVector(abs(n))
  CharacterVector colnam = names ? no_init_vector(ns) : no_init_vector(1); 
  Matrix<RTYPE> out = no_init_matrix(l, ns); 
  if(ng == 0) { // No groups
    if(Rf_isNull(t)) { // Ordered data
      for(int p = ns; p--; ) {
        int np = n[p];
        if(absn[p] > l) stop("lag-length exceeds length of vector");
        MatrixColumn<RTYPE> outp = out( _ , p);
        if(np>0) {
          if(names) colnam[p] = ".L" + nc[p]; 
          int i = 0;
          while(i != np) outp[i++] = ff;
          for( ; i != l; ++i) outp[i] = x[i - np];
        } else if(np<0) {
          if(names) colnam[p] = ".F" + nc[p]; 
          int i = l, st = l+np;
          while(i != st) outp[--i] = ff;
          for( ; i--; ) outp[i] = x[i - np];
        } else {
          if(names) colnam[p] = ".--"; 
          outp = x;
        }
      }
    } else { // Unordered data: Timevar provided 
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
      for(int p = ns; p--; ) {
        int np = n[p];
        if(absn[p] > l) stop("lag-length exceeds length of vector");
        MatrixColumn<RTYPE> outp = out( _ , p);
        if(np>0) {
          if(names) colnam[p] = ".L" + nc[p]; 
          int i = 0;
          while(i != np) outp[omap[i++]] = ff; 
          for( ; i != l; ++i) outp[omap[i]] = x[omap[i - np]]; 
        } else if(np<0) {
          if(names) colnam[p] = ".F" + nc[p]; 
          int st = l+np, i = l; 
          while(i != st) outp[omap[--i]] = ff;
          for( ; i--; ) outp[omap[i]] = x[omap[i - np]]; 
        } else {
          if(names) colnam[p] = ".--";
          outp = x;
        }
      }
    }
  } else { // With groups
    if(l != g.size()) stop("length(x) must match length(g)");
    int ags = l/ng, ngp = ng+1;
    if(Rf_isNull(t)) { // Ordered data
      int seen[ngp], memsize = sizeof(int)*ngp;
      for(int p = ns; p--; ) {
        int np = n[p];
        if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
        MatrixColumn<RTYPE> outp = out( _ , p);
        if(np>0) {
          if(names) colnam[p] = ".L" + nc[p]; 
          memset(seen, 0, memsize); 
          for(int i = 0; i != l; ++i) {
            if(seen[g[i]] == np) {
              outp[i] = x[i-np];
            } else {
              outp[i] = ff;
              ++seen[g[i]];
            }
          }
        } else if(np<0) {
          memset(seen, 0, memsize); 
          if(names) colnam[p] = ".F" + nc[p]; 
          for(int i = l; i--; ) { // good??
            if(seen[g[i]] == np) {
              outp[i] = x[i-np];
            } else {
              outp[i] = ff;
              --seen[g[i]];
            }
          }
        } else {
          if(names) colnam[p] = ".--"; 
          outp = x;
        }
      }
    } else { // Unordered data: Timevar provided
      IntegerVector ord = t;
      if(l != ord.size()) stop("length(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX); // INFINITY gives bug !!!!!!!!
      IntegerVector gsv = NULL; //no_init_vector(ng); // No real improvements here by using C++ arrays !!
      IntegerVector ord2 = no_init_vector(l); // use array ????????????????? 
      if(Rf_isNull(gs)) {
        // std::fill(gsv.begin(), gsv.end(), 0);
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
      int cgs[ngp]; 
      cgs[1] = 0;
      for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; // or get "starts from forderv"
      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]]; // still room for speed improvement ?? ???????????
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        if(omap[cgs[g[i]]+ord2[i]] == 0) omap[cgs[g[i]]+ord2[i]] = i; 
        else stop("Repeated values of timevar within one or more groups"); 
      }
      for(int p = ns; p--; ) {
        int np = n[p];
        if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
        MatrixColumn<RTYPE> outp = out( _ , p);
        if(np>0) {
          if(names) colnam[p] = ".L" + nc[p]; 
          for(int i = 0; i != l; ++i) {
            if(ord2[i] >= np) {
              outp[i] = x[omap[cgs[g[i]]+ord2[i]-np]];
            } else {
              outp[i] = ff;
            }
          }
        } else if(np<0) {
          if(names) colnam[p] = ".F" + nc[p]; 
          for(int i = 0; i != l; ++i) { 
            if(ord2[i] < gsv[g[i]-1]+np) {
              outp[i] = x[omap[cgs[g[i]]+ord2[i]-np]]; 
            } else {
              outp[i] = ff;
            }
          }
        } else {
          if(names) colnam[p] = ".--"; 
          outp = x;
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x); 
  if(ns != 1) { 
    if(x.hasAttribute("names")) out.attr("names") = R_NilValue; // fastest ??
    out.attr("dim") = Dimension(l, ns);
    if(names) out.attr("dimnames") = List::create(x.attr("names"), colnam);
    // out.attr("class") = CharacterVector::create(x.attr("class"),"matrix");
  }
  return out;
}

// template <>
// Vector<CPLXSXP> flagleadCppImpl(Vector<CPLXSXP> x, IntegerVector n, SEXP fill,
//                               int ng, IntegerVector g, SEXP gs, SEXP t) {
//   stop("Not supported SEXP type!");
// }

template <>
Vector<VECSXP> flagleadCppImpl(const Vector<VECSXP>& x, const IntegerVector& n, const SEXP& fill,
                               int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Vector<RAWSXP> flagleadCppImpl(const Vector<RAWSXP>& x, const IntegerVector& n, const SEXP& fill,
                               int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Vector<EXPRSXP> flagleadCppImpl(const Vector<EXPRSXP>& x, const IntegerVector& n, const SEXP& fill,
                                int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP flagleadCpp(SEXP x, IntegerVector n = 1, SEXP fill = R_NilValue, 
                 int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue, bool names = true){
  RCPP_RETURN_VECTOR(flagleadCppImpl, x, n, fill, ng, g, gs, t, names);
}


// 5th version = 3rd version + naming !! -> 1D omap + ng+1 trick -> fastest !!
// // [[Rcpp::export]]
// NumericVector flagleadCpp(NumericVector x, IntegerVector n = 1, double fill = NA_REAL,
//                           int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) {
//   
//   int l = x.size(), ns = n.size();
//   IntegerVector absn = no_init_vector(ns);
//   for(int i = 0; i != ns; ++i) {
//     if(n[i]<0) absn[i] = -n[i];
//     else absn[i] = n[i];
//   }
//   CharacterVector nc = Rf_coerceVector(absn, STRSXP);  // NumericVector(abs(n))
//   CharacterVector colnam = no_init_vector(ns); 
//   NumericMatrix out = no_init_matrix(l, ns); 
//   if(ng == 0) { // No groups
//     if(Rf_isNull(t)) { // Ordered data
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         if(absn[p] > l) stop("lag-length exceeds length of vector");
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           colnam[p] = ".L" + nc[p]; 
//           int i = 0;
//           while(i != np) outp[i++] = fill;
//           for( ; i != l; ++i) outp[i] = x[i - np];
//         } else if(np<0) {
//           colnam[p] = ".F" + nc[p]; 
//           int i = l, st = l+np;
//           while(i != st) outp[--i] = fill;
//           for( ; i--; ) outp[i] = x[i - np];
//         } else {
//           colnam[p] = ".--"; 
//           outp = x;
//         }
//       }
//     } else { // Unordered data: Timevar provided 
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       LogicalVector ocheck(l, true); 
//       IntegerVector omap = no_init_vector(l); 
//       for(int i = 0; i != l; ++i) { 
//         if(ord[i] > l) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
//         if(ocheck[ord[i]-1]) {
//           ocheck[ord[i]-1] = false;
//           omap[ord[i]-1] = i; // Note: omap is the same as order(ord) !!
//         } else {
//           stop("Repeated values in timevar");
//         }
//       }
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         if(absn[p] > l) stop("lag-length exceeds length of vector");
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           colnam[p] = ".L" + nc[p]; 
//           int i = 0;
//           while(i != np) outp[omap[i++]] = fill; 
//           for( ; i != l; ++i) outp[omap[i]] = x[omap[i - np]]; 
//         } else if(np<0) {
//           colnam[p] = ".F" + nc[p]; 
//           int st = l+np, i = l; 
//           while(i != st) outp[omap[--i]] = fill;
//           for( ; i--; ) outp[omap[i]] = x[omap[i - np]]; 
//         } else {
//           colnam[p] = ".--";
//           outp = x;
//         }
//       }
//     }
//   } else { // With groups
//     if(l != g.size()) stop("length(x) must match length(g)");
//     int ags = l/ng, ngp = ng+1;
//     if(Rf_isNull(t)) { // Ordered data
//       int seen[ngp], memsize = sizeof(int)*ngp;
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           colnam[p] = ".L" + nc[p]; 
//           memset(seen, 0, memsize); 
//           for(int i = 0; i != l; ++i) {
//             if(seen[g[i]] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               ++seen[g[i]];
//             }
//           }
//         } else if(np<0) {
//           memset(seen, 0, memsize); 
//           colnam[p] = ".F" + nc[p]; 
//           for(int i = l; i--; ) { // good??
//             if(seen[g[i]] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               --seen[g[i]];
//             }
//           }
//         } else {
//           colnam[p] = ".--"; 
//           outp = x;
//         }
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       IntegerVector min(ngp, INT_MAX);
//       IntegerVector gsv = NULL; //no_init_vector(ng); // No real improvements here by using C++ arrays !!
//       IntegerVector ord2 = no_init_vector(l);
//       if(Rf_isNull(gs)) {
//         // std::fill(gsv.begin(), gsv.end(), 0);
//         gsv = IntegerVector(ng);
//         for(int i = 0; i != l; ++i) {
//           ++gsv[g[i]-1];
//           if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)");
//         for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//       }
//       IntegerVector omap(l); 
//       int cgs[ngp]; 
//       cgs[1] = 0;
//       for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; // or get "starts from forderv"
//       for(int i = 0; i != l; ++i) {
//         ord2[i] = ord[i] - min[g[i]]; // still room for speed improvement ?? ???????????
//         if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//         if(omap[cgs[g[i]]+ord2[i]] == 0) omap[cgs[g[i]]+ord2[i]] = i; 
//         else stop("Repeated values of timevar within one or more groups"); 
//       }
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           colnam[p] = ".L" + nc[p]; 
//           for(int i = 0; i != l; ++i) {
//             if(ord2[i] >= np) {
//               outp[i] = x[omap[cgs[g[i]]+ord2[i]-np]];
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else if(np<0) {
//           colnam[p] = ".F" + nc[p]; 
//           for(int i = 0; i != l; ++i) { 
//             if(ord2[i] < gsv[g[i]-1]+np) {
//               outp[i] = x[omap[cgs[g[i]]+ord2[i]-np]]; 
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else {
//           colnam[p] = ".--"; 
//           outp = x;
//         }
//       }
//     }
//   }
//   if(ns == 1) SHALLOW_DUPLICATE_ATTRIB(out, x);
//   else colnames(out) = colnam;
//   return out;
// }


// 4th version: 1D omap and saving index trick. -> but less efficient than 3rd version on most data !!
// // [[Rcpp::export]]
// NumericVector flagleadCpp3(NumericVector x, IntegerVector n = 1, double fill = NA_REAL,
//                            int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) {
//   
//   int l = x.size(), ns = n.size();
//   NumericMatrix out = no_init_matrix(l, ns); // init better ??
//   if(ng == 0) { // No groups
//     if(Rf_isNull(t)) { // Ordered data
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > l) stop("lag-length exceeds length of vector");
//           int i = 0;
//           while(i != np) outp[i++] = fill;
//           for( ; i != l; ++i) outp[i] = x[i - np];
//         } else if(np<0) {
//           if(-np > l) stop("lag-length exceeds length of vector");
//           int i = l, st = l+np;
//           while(i != st) outp[--i] = fill;
//           for( ; i--; ) outp[i] = x[i - np];
//         } else outp = x;
//       }
//     } else { // Unordered data: Timevar provided -> works ???
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       LogicalVector ocheck(l, true); // Could do more efficiently if only one value -> but you don't use this anyway!!
//       IntegerVector omap = no_init_vector(l); // faster with C++ arrays ?? // int omap[l]; -> The omap array way the source of instability !!!!!!!!!
//       for(int i = 0; i != l; ++i) { // integrated below !! -> fastest ??
//         if(ord[i] > l) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
//         if(ocheck[ord[i]-1]) {
//           ocheck[ord[i]-1] = false;
//           omap[ord[i]-1] = i; // Note: omap is the same as order(ord) !!
//         } else {
//           stop("Repeated values in timevar");
//         }
//       }
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > l) stop("lag-length exceeds length of vector");
//           int i = 0;
//           while(i != np) outp[omap[i++]] = fill; //seems instable !!
//           for( ; i != l; ++i) outp[omap[i]] = x[omap[i - np]]; // faster than below ??
//         } else if(np<0) {
//           if(-np > l) stop("lag-length exceeds length of vector");
//           int st = l+np, i = l; // seems instable
//           while(i != st) outp[omap[--i]] = fill;
//           for( ; i--; ) outp[omap[i]] = x[omap[i - np]]; // faster than below ??
//         } else outp = x;
//       }
//     }
//   } else { // With groups
//     if(l != g.size()) stop("length(x) must match length(g)");
//     int ags = l/ng, ngp = ng+1;
//     if(Rf_isNull(t)) { // Ordered data
//       int seen[ngp], memsize = sizeof(int)*ngp;
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > ags) stop("lag-length exceeds average group-size");
//           memset(seen, 0, memsize); // fastest !!
//           for(int i = 0; i != l; ++i) {
//             if(seen[g[i]] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               ++seen[g[i]];
//             }
//           }
//         } else if(np<0) {
//           if(-np > ags) stop("lag-length exceeds average group-size");
//           memset(seen, 0, memsize); // fastest !!
//           for(int i = l; i--; ) { // good??
//             if(seen[g[i]] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               --seen[g[i]];
//             }
//           }
//         } else outp = x;
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       IntegerVector min(ngp, INT_MAX);
//       IntegerVector gsv = NULL; //no_init_vector(ng); // No real improvements here by using C++ arrays !!
//       IntegerVector ord2 = no_init_vector(l);
//       if(Rf_isNull(gs)) {
//         // std::fill(gsv.begin(), gsv.end(), 0);
//         gsv = IntegerVector(ng);
//         for(int i = 0; i != l; ++i) {
//           ++gsv[g[i]-1];
//           if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)");
//         for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//       }
//       IntegerVector omap(l), index = no_init_vector(l); 
//       int cgs[ngp]; // Stable ???? //, index = 0; // seen[ngp], memsize = sizeof(int)*(ngp) // i2 = 0; using i2 is not faster
//       cgs[1] = 0;
//       for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; // or get "starts from forderv"
//       for(int i = 0; i != l; ++i) {
//         ord2[i] = ord[i] - min[g[i]]; 
//         if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//         index[i] = cgs[g[i]]+ord2[i]; // fastest ?? -> Nope !!! 
//         if(omap[index[i]] == 0) omap[index[i]] = i; // fastest ?? -> What about saving  cgs[g[i]]+ord2[i] in another vector ??? 
//         else stop("Repeated values of timevar within one or more groups"); 
//       }
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > ags) stop("lag-length exceeds average group-size");
//           for(int i = 0; i != l; ++i) {
//             if(ord2[i] >= np) {
//               outp[i] = x[omap[index[i]-np]]; // old way of access !! -> much better !!
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else if(np<0) {
//           if(-np > ags) stop("lag-length exceeds average group-size");
//           for(int i = 0; i != l; ++i) { // best loop ??
//             if(ord2[i] < gsv[g[i]-1]+np) {
//               outp[i] = x[omap[index[i]-np]]; // old way of access !!
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else outp = x;
//       }
//     }
//   }
//   if(ns == 1) SHALLOW_DUPLICATE_ATTRIB(out, x);
//   return out;
// }

// // 3rd version -> 1D omap + ng+1 trick -> fastest !!
// // [[Rcpp::export]]
// NumericVector flagleadCpp2(NumericVector x, IntegerVector n = 1, double fill = NA_REAL,
//                            int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) {
//   
//   int l = x.size(), ns = n.size();
//   NumericMatrix out = no_init_matrix(l, ns); // init better ??
//   if(ng == 0) { // No groups
//     if(Rf_isNull(t)) { // Ordered data
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > l) stop("lag-length exceeds length of vector");
//           int i = 0;
//           while(i != np) outp[i++] = fill;
//           for( ; i != l; ++i) outp[i] = x[i - np];
//         } else if(np<0) {
//           if(-np > l) stop("lag-length exceeds length of vector");
//           int i = l, st = l+np;
//           while(i != st) outp[--i] = fill;
//           for( ; i--; ) outp[i] = x[i - np];
//         } else outp = x;
//       }
//     } else { // Unordered data: Timevar provided -> works ???
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       LogicalVector ocheck(l, true); // Could do more efficiently if only one value -> but you don't use this anyway!!
//       IntegerVector omap = no_init_vector(l); // faster with C++ arrays ?? // int omap[l]; -> The omap array way the source of instability !!!!!!!!!
//       for(int i = 0; i != l; ++i) { // integrated below !! -> fastest ??
//         if(ord[i] > l) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
//         if(ocheck[ord[i]-1]) {
//           ocheck[ord[i]-1] = false;
//           omap[ord[i]-1] = i; // Note: omap is the same as order(ord) !!
//         } else {
//           stop("Repeated values in timevar");
//         }
//       }
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > l) stop("lag-length exceeds length of vector");
//           
//           // check out what is faster ?? 
//           int i = 0;
//           while(i != np) outp[omap[i++]] = fill; //seems instable !!
//           for( ; i != l; ++i) outp[omap[i]] = x[omap[i - np]]; // faster than below ??
//           // for(int i = np; i != l; ++i) outp[omap[i]] = x[omap[i - np]]; 
//           // for(int i = 0; i != l; ++i) { // This is also instable with large data !!!
//           //   if(i >= np) {
//           //     outp[omap[i]] = x[omap[i - np]];
//           //   } else {
//           //     outp[omap[i]] = fill;
//           //   }
//           // }
//           // for(int i = 0; i != l; ++i) { // also works with small data but gives error on large !!
//           //   if(ord[i] > np) {
//           //     outp[i] = x[omap[ord[i]-np-1]];
//           //   } else {
//           //     outp[i] = fill;
//           //   }
//           // }
//         } else if(np<0) {
//           if(-np > l) stop("lag-length exceeds length of vector");
//           int st = l+np, i = l; // seems instable
//           while(i != st) outp[omap[--i]] = fill;
//           for( ; i--; ) outp[omap[i]] = x[omap[i - np]]; // faster than below ??
//           // for(int i = 0; i != l; ++i) { // best loop ?? (i.e. fastest ??)
//           //   if(i < st) { // right ?
//           //     outp[omap[i]] = x[omap[i - np]];
//           //   } else {
//           //     outp[omap[i]] = fill;
//           //   }
//           // }
//           // for(int i = 0; i != l; ++i) { // also works with small data but gives error on large !!
//           //   if(ord[i] <= st) {
//           //     outp[i] = x[omap[ord[i]-np-1]];
//           //   } else {
//           //     outp[i] = fill;
//           //   }
//           // }
//         } else outp = x;
//       }
//     }
//   } else { // With groups
//     if(l != g.size()) stop("length(x) must match length(g)");
//     int ags = l/ng, ngp = ng+1;
//     if(Rf_isNull(t)) { // Ordered data
//       int seen[ngp], memsize = sizeof(int)*ngp;
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > ags) stop("lag-length exceeds average group-size");
//           memset(seen, 0, memsize); // fastest !!
//           for(int i = 0; i != l; ++i) {
//             if(seen[g[i]] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               ++seen[g[i]];
//             }
//           }
//         } else if(np<0) {
//           if(-np > ags) stop("lag-length exceeds average group-size");
//           memset(seen, 0, memsize); // fastest !!
//           for(int i = l; i--; ) { // good??
//             if(seen[g[i]] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               --seen[g[i]];
//             }
//           }
//         } else outp = x;
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       IntegerVector min(ngp, INT_MAX);
//       IntegerVector gsv = NULL; //no_init_vector(ng); // No real improvements here by using C++ arrays !!
//       IntegerVector ord2 = no_init_vector(l);
//       if(Rf_isNull(gs)) {
//         // std::fill(gsv.begin(), gsv.end(), 0);
//         gsv = IntegerVector(ng);
//         for(int i = 0; i != l; ++i) {
//           ++gsv[g[i]-1];
//           if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)");
//         for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//       }
//       IntegerVector omap(l); 
//       int cgs[ngp]; // Stable ???? //, index = 0; // seen[ngp], memsize = sizeof(int)*(ngp) // i2 = 0; using i2 is not faster
//       cgs[1] = 0;
//       for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; // or get "starts from forderv"
//       for(int i = 0; i != l; ++i) {
//         ord2[i] = ord[i] - min[g[i]]; 
//         if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//         // index = cgs[g[i]]+ord2[i]; // fastest ?? -> Nope !!! 
//         if(omap[cgs[g[i]]+ord2[i]] == 0) omap[cgs[g[i]]+ord2[i]] = i; // fastest ?? -> What about saving  cgs[g[i]]+ord2[i] in another vector ??? 
//         else stop("Repeated values of timevar within one or more groups"); 
//       }
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > ags) stop("lag-length exceeds average group-size");
//           for(int i = 0; i != l; ++i) {
//             if(ord2[i] >= np) {
//               outp[i] = x[omap[cgs[g[i]]+ord2[i]-np]]; // old way of access !! -> much better !!
//             } else {
//               outp[i] = fill;
//             }
//           }
//           // memset(seen, 0, memsize); 
//           // for(int i = 0; i != l; ++i) {
//           //   if(seen[g[omap[i]]] == np) {
//           //     outp[omap[i]] = x[omap[i - np]];
//           //   } else {
//           //     outp[omap[i]] = fill;
//           //     ++seen[g[omap[i]]];
//           //   }
//           // }
//         } else if(np<0) {
//           if(-np > ags) stop("lag-length exceeds average group-size");
//           for(int i = 0; i != l; ++i) { // best loop ??
//             if(ord2[i] < gsv[g[i]-1]+np) {
//               outp[i] = x[omap[cgs[g[i]]+ord2[i]-np]]; // old way of access !!
//             } else {
//               outp[i] = fill;
//             }
//           }
//           // memset(seen, 0, memsize); 
//           // for(int i = l; i--; ) { // good??
//           //   if(seen[g[omap[i]]] == np) {
//           //     outp[omap[i]] = x[omap[i - np]];
//           //   } else {
//           //     outp[omap[i]] = fill;
//           //     --seen[g[omap[i]]];
//           //   }
//           // }
//         } else outp = x;
//       }
//     }
//   }
//   if(ns == 1) SHALLOW_DUPLICATE_ATTRIB(out, x);
//   return out;
// }


// Second versin unordered panel lag checking the data, 2d omap !!
// // [[Rcpp::export]]
// NumericVector flagleadCpp1(NumericVector x, IntegerVector n = 1, double fill = NA_REAL,
//                            int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) {
//   
//   int l = x.size(), ns = n.size();
//   NumericMatrix out = no_init_matrix(l, ns); // init better ??
//   if(ng == 0) { // No groups
//     if(Rf_isNull(t)) { // Ordered data
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > l) stop("lag-length exceeds length of vector");
//           int i = 0;
//           while(i != np) outp[i++] = fill;
//           for( ; i != l; ++i) outp[i] = x[i - np];
//         } else if(np<0) {
//           if(-np > l) stop("lag-length exceeds length of vector");
//           int i = l, st = l+np;
//           while(i != st) outp[--i] = fill;
//           for( ; i--; ) outp[i] = x[i - np];
//         } else outp = x;
//       }
//     } else { // Unordered data: Timevar provided -> works ???
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       LogicalVector ocheck(l, true); // Could do more efficiently if only one value -> but you don't use this anyway!!
//       for(int i = 0; i != l; ++i) {
//         if(ocheck[ord[i]-1]) ocheck[ord[i]-1] = false;
//         else stop("Repeated values in timevar");
//       }
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > l) stop("lag-length exceeds length of vector");
//           for(int i = 0; i != l; ++i) {
//             if(ord[i] > np) {
//               outp[i] = x[ord[i]-np-1];
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else if(np<0) {
//           if(-np > l) stop("lag-length exceeds length of vector");
//           int st = l+np;
//           for(int i = 0; i != l; ++i) { // best loop ?? (i.e. fastest ??)
//             if(ord[i] <= st) {
//               outp[i] = x[ord[i]-np-1];
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else outp = x;
//       }
//     }
//   } else { // With groups
//     if(l != g.size()) stop("length(x) must match length(g)");
//     int ags = l/ng, ngp = ng+1;
//     if(Rf_isNull(t)) { // Ordered data
//       int seen[ngp], memsize = sizeof(int)*ngp;
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > ags) stop("lag-length exceeds average group-size");
//           memset(seen, 0, memsize); // fastest !!
//           for(int i = 0; i != l; ++i) {
//             if(seen[g[i]] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               ++seen[g[i]];
//             }
//           }
//         } else if(np<0) {
//           if(-np > ags) stop("lag-length exceeds average group-size");
//           memset(seen, 0, memsize); // fastest !!
//           for(int i = l; i--; ) { // good??
//             if(seen[g[i]] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               --seen[g[i]];
//             }
//           }
//         } else outp = x;
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       IntegerVector min(ngp, INT_MAX);
//       IntegerVector gsv = NULL; // No real improvements here by using C++ arrays !!
//       IntegerVector ord2 = no_init_vector(l);
//       if(Rf_isNull(gs)) {
//         // std::fill(gsv.begin(), gsv.end(), 0);
//         gsv = IntegerVector(ng);
//         for(int i = 0; i != l; ++i) {
//           ++gsv[g[i]-1];
//           if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)");
//         for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//       }
//       int **omap = new int*[ngp];
//       for(int i = 0; i != ng; ++i) omap[i+1] = new int[gsv[i]]{}; // or () 
//       for(int i = 0; i != l; ++i) {
//         ord2[i] = ord[i] - min[g[i]];
//         if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//         if(omap[g[i]][ord2[i]] == 0) omap[g[i]][ord2[i]] = i; // fastest ??
//         else stop("Repeated values of timevar within one or more groups");
//       }
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > ags) stop("lag-length exceeds average group-size");
//           for(int i = 0; i != l; ++i) {
//             if(ord2[i] >= np) {
//               outp[i] = x[omap[g[i]][ord2[i]-np]];
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else if(np<0) {
//           if(-np > ags) stop("lag-length exceeds average group-size");
//           for(int i = 0; i != l; ++i) { // best loop ??
//             if(ord2[i] < gsv[g[i]-1]+np) {
//               outp[i] = x[omap[g[i]][ord2[i]-np]];
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else outp = x;
//       }
//     }
//   }
//   if(ns == 1) SHALLOW_DUPLICATE_ATTRIB(out, x);
//   return out;
// }


// First version with unordered panel lag checking the data , but 2d omap !!
// // [[Rcpp::export]]
// NumericVector flagleadCpp(NumericVector x, IntegerVector n = 1, double fill = NA_REAL,
//                           int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) {
// 
//   int l = x.size(), ns = n.size();
//   NumericMatrix out = no_init_matrix(l, ns); // init better ??
//   if(ng == 0) { // No groups
//     if(Rf_isNull(t)) { // Ordered data
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > l) stop("lag-length exceeds length of vector");
//           int i = 0;
//           while(i != np) outp[i++] = fill;
//           for( ; i != l; ++i) outp[i] = x[i - np];
//         } else if(np<0) {
//           if(-np > l) stop("lag-length exceeds length of vector");
//           int i = l, st = l+np;
//           while(i != st) outp[--i] = fill;
//           for( ; i--; ) outp[i] = x[i - np];
//         } else outp = x;
//       }
//     } else { // Unordered data: Timevar provided -> works ???
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       LogicalVector ocheck(l, true); // Could do more efficiently if only one value -> but you don't use this anyway!!
//       for(int i = 0; i != l; ++i) {
//         if(ocheck[ord[i]-1]) ocheck[ord[i]-1] = false;
//         else stop("Repeated values in timevar");
//       }
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > l) stop("lag-length exceeds length of vector");
//           for(int i = 0; i != l; ++i) {
//             if(ord[i] > np) {
//               outp[i] = x[ord[i]-np-1];
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else if(np<0) {
//           if(-np > l) stop("lag-length exceeds length of vector");
//           int st = l+np;
//           for(int i = 0; i != l; ++i) { // best loop ?? (i.e. fastest ??)
//             if(ord[i] <= st) {
//               outp[i] = x[ord[i]-np-1];
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else outp = x;
//       }
//     }
//   } else { // With groups
//     if(l != g.size()) stop("length(x) must match length(g)");
//     int ags = l/ng;
//     if(Rf_isNull(t)) { // Ordered data
//       int seen[ng], memsize = sizeof(int)*ng;
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > ags) stop("lag-length exceeds average group-size");
//           memset(seen, 0, memsize); // fastest !!
//           for(int i = 0; i != l; ++i) {
//             if(seen[g[i]-1] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               ++seen[g[i]-1];
//             }
//           }
//         } else if(np<0) {
//           if(-np > ags) stop("lag-length exceeds average group-size");
//           memset(seen, 0, memsize); // fastest !!
//           for(int i = l; i--; ) { // good??
//             if(seen[g[i]-1] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               --seen[g[i]-1];
//             }
//           }
//         } else outp = x;
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       IntegerVector min(ng, INT_MAX);
//       IntegerVector gsv = no_init_vector(ng); // No real improvements here by using C++ arrays !!
//       IntegerVector ord2 = no_init_vector(l);
//       if(Rf_isNull(gs)) {
//         std::fill(gsv.begin(), gsv.end(), 0);
//         for(int i = 0; i != l; ++i) {
//           ++gsv[g[i]-1];
//           if(ord[i] < min[g[i]-1]) min[g[i]-1] = ord[i];
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)");
//         for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]-1]) min[g[i]-1] = ord[i];
//       }
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) omap[i] = new int[gsv[i]]{}; // or () // better using vector of vectors (faster initialization) ??
//       for(int i = 0; i != l; ++i) {
//         ord2[i] = ord[i] - min[g[i]-1];
//         if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//         if(omap[g[i]-1][ord2[i]] == 0) omap[g[i]-1][ord2[i]] = i; // fastest ??
//         else stop("Repeated values of timevar within one or more groups");
//       }
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           if(np > ags) stop("lag-length exceeds average group-size");
//           for(int i = 0; i != l; ++i) {
//             if(ord2[i] >= np) {
//               outp[i] = x[omap[g[i]-1][ord2[i]-np]];
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else if(np<0) {
//           if(-np > ags) stop("lag-length exceeds average group-size");
//           for(int i = 0; i != l; ++i) { // best loop ??
//             if(ord2[i] < gsv[g[i]-1]+np) {
//               outp[i] = x[omap[g[i]-1][ord2[i]-np]];
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else outp = x;
//       }
//     }
//   }
//   if(ns == 1) SHALLOW_DUPLICATE_ATTRIB(out, x);
//   return out;
// }



// Previous Version: Old Unordered panel-lag!! - without min !!
// // [[Rcpp::export]]
// NumericVector flagleadCpp(NumericVector x, IntegerVector n = 1, double fill = NA_REAL, 
//                           int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, SEXP t = R_NilValue) { 
//   
//   int l = x.size(), ns = n.size();
//   NumericMatrix out = no_init_matrix(l, ns); // init better ??
//   if(ng == 0) { // No groups 
//     if (Rf_isNull(t)) { // Ordered data
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           int i = 0; 
//           while(i != np) outp[i++] = fill;
//           for( ; i != l; ++i) outp[i] = x[i - np];
//         } else if(np<0) {
//           int i = l, st = l+np; 
//           while(i != st) outp[--i] = fill;
//           for( ; i--; ) outp[i] = x[i - np];
//         } else outp = x;
//       }
//     } else { // Unordered data: Timevar provided -> works ???
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           for(int i = 0; i != l; ++i) {
//             if(ord[i] > np) {
//               outp[i] = x[ord[i]-np-1]; 
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else if(np<0) {
//           int st = l+np;
//           for(int i = 0; i != l; ++i) { // best loop ?? (i.e. fastest ??)
//             if(ord[i] <= st) {
//               outp[i] = x[ord[i]-np-1]; 
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else outp = x;
//       }
//     }
//   } else { // With groups
//     if(l != g.size()) stop("length(x) must match length(g)");
//     if(Rf_isNull(t)) { // Ordered data
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         IntegerVector seenp(ng); // No init possible ??
//         if(np>0) {
//           for(int i = 0; i != l; ++i) { 
//             if(seenp[g[i]-1] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               ++seenp[g[i]-1];
//             }
//           }
//         } else if(np<0) {
//           for(int i = l; i--; ) { // good?? 
//             if(seenp[g[i]-1] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               --seenp[g[i]-1];
//             }
//           }
//         } else outp = x;
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t; 
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       if(ng != gs.size()) stop("ng must match length(gs)");
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//       for(int i = 0; i != l; ++i) omap[g[i]-1][ord[i]-1] = i;
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           for(int i = 0; i != l; ++i) {
//             if(ord[i] > np) {
//               outp[i] = x[omap[g[i]-1][ord[i]-np-1]]; 
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else if(np<0) {
//           // int tgt = gs+np+1; // good ??? Best way??, could do without this vector !!
//           for(int i = 0; i != l; ++i) { // best loop ??
//             if(ord[i] <= gs[g[i]-1]+np) {
//               outp[i] = x[omap[g[i]-1][ord[i]-np-1]]; 
//             } else {
//               outp[i] = fill;
//             }
//           }
//         } else outp = x;
//       }
//     }
//   }
//   return out;
// }


// Previous Version: distinguishing between one and more lags !! -> Does not give spped improvement !!
// // [[Rcpp::export]]
// NumericVector flagleadCpp(NumericVector x, IntegerVector n = 1, double fill = NA_REAL, 
//                            int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, SEXP t = R_NilValue) { 
//   int l = x.size(), ns = n.size();
//   
//   if(ns == 1) { // Only one lag !! (Need extra method ??)
//     int nn = n[0]; // good ??
//     NumericVector out = no_init_vector(l); // init better ??
//     if(ng == 0) { // No groups !!
//       if(nn>0) {
//         int i = 0; 
//         while(i != nn) out[i++] = fill;
//         for( ; i != l; ++i) out[i] = x[i - nn];
//       } else if(nn<0) {
//         int i = l, st = l+nn; 
//         while(i != st) out[--i] = fill;
//         for( ; i--; ) out[i] = x[i - nn];
//       } else return x;
//     } else {
//       if(l != g.size()) stop("length(x) must match length(g)");
//       if (Rf_isNull(t)) { // Ordered data
//         IntegerVector seen(ng);
//         if(nn>0) {
//           for(int i = 0; i != l; ++i) {  
//             if(seen[g[i]-1] == nn) {
//               out[i] = x[i-nn];
//             } else {
//               out[i] = fill;
//               ++seen[g[i]-1];
//             }
//           }
//         } else if(nn<0) {
//           for(int i = l; i--; ) { // good?? 
//             if(seen[g[i]-1] == nn) {
//               out[i] = x[i-nn];
//             } else {
//               out[i] = fill;
//               --seen[g[i]-1];
//             }
//           }
//         } else return x;
//       } else { // Unordered data: Timevar provided
//         IntegerVector ord = t; 
//         if(l != ord.size()) stop("length(x) must match length(t)");
//         if(ng != gs.size()) stop("ng must match length(gs)");
//         int **omap = new int*[ng];
//         for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//         for(int i = 0; i != l; ++i) omap[g[i]-1][ord[i]-1] = i;
//         if(nn>0) {
//           for(int i = 0; i != l; ++i) {
//             if(ord[i] > nn) {
//               out[i] = x[omap[g[i]-1][ord[i]-nn-1]]; // good ??? -> Yes, great and efficient !!!
//             } else {
//               out[i] = fill;
//             }
//           }
//         } else if(nn<0) {
//           // IntegerVector tgt = gs+nn+1; // good ??? Best way
//           for(int i = 0; i != l; ++i) { // best loop ??
//             if(ord[i] <= gs[g[i]-1]+nn) { // easier way, rather a counting vector by group ???
//               out[i] = x[omap[g[i]-1][ord[i]-nn-1]]; // good ??? -> Yes, great and efficient !!!
//             } else {
//               out[i] = fill;
//             }
//           }
//         } else return x;
//       }
//     }
//     return out;
//   } else { // Multiple lags !!!
//     NumericMatrix out = no_init_matrix(l, ns); // init better ??
//     if(ng == 0) { // No groups !!
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           int i = 0; 
//           while(i != np) outp[i++] = fill;
//           for( ; i != l; ++i) outp[i] = x[i - np];
//         } else if(np<0) {
//           int i = l, st = l+np; 
//           while(i != st) outp[--i] = fill;
//           for( ; i--; ) outp[i] = x[i - np];
//         } else outp = x;
//       }
//     } else {
//       if(l != g.size()) stop("length(x) must match length(g)");
//       if (Rf_isNull(t)) { // Ordered data
//         for(int p = ns; p--; ) {
//           int np = n[p];
//           NumericMatrix::Column outp = out( _ , p);
//           IntegerVector seenp(ng); // No init possible ??
//           if(np>0) {
//             for(int i = 0; i != l; ++i) { 
//               if(seenp[g[i]-1] == np) {
//                 outp[i] = x[i-np];
//               } else {
//                 outp[i] = fill;
//                 ++seenp[g[i]-1];
//               }
//             }
//           } else if(np<0) {
//             for(int i = l; i--; ) { // good?? 
//               if(seenp[g[i]-1] == np) {
//                 outp[i] = x[i-np];
//               } else {
//                 outp[i] = fill;
//                 --seenp[g[i]-1];
//               }
//             }
//           } else outp = x;
//         }
//       } else { // Unordered data: Timevar provided
//         IntegerVector ord = t; 
//         if(l != ord.size()) stop("length(x) must match length(t)");
//         if(ng != gs.size()) stop("ng must match length(gs)");
//         int **omap = new int*[ng];
//         for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//         for(int i = 0; i != l; ++i) omap[g[i]-1][ord[i]-1] = i;
//         for(int p = ns; p--; ) {
//           int np = n[p];
//           NumericMatrix::Column outp = out( _ , p);
//           if(np>0) {
//             for(int i = 0; i != l; ++i) {
//               if(ord[i] > np) {
//                 outp[i] = x[omap[g[i]-1][ord[i]-np-1]]; 
//               } else {
//                 outp[i] = fill;
//               }
//             }
//           } else if(np<0) {
//             // int tgt = gs+np+1; // good ??? Best way??, could do without this vector !!
//             for(int i = 0; i != l; ++i) { // best loop ??
//               if(ord[i] <= gs[g[i]-1]+np) {
//                 outp[i] = x[omap[g[i]-1][ord[i]-np-1]]; 
//               } else {
//                 outp[i] = fill;
//               }
//             }
//           } else outp = x;
//         }
//       }
//     }
//     return out;
//   }
// }


// Previous Version: Using a vector for gs if t supplied !!
// // [[Rcpp::export]]
// NumericVector flagleadCpp(NumericVector x, IntegerVector n = 1, double fill = NA_REAL, 
//                       int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, SEXP t = R_NilValue) { 
//   int l = x.size(), ns = n.size();
//   
//   if(ns == 1) { // Only one lag !! (Need extra method ??)
//     int nn = n[0]; // good ??
//     NumericVector out = no_init_vector(l); // init better ??
//     if(ng == 0) { // No groups !!
//       if(nn>0) {
//         int i = 0; 
//         while(i != nn) out[i++] = fill;
//         for( ; i != l; ++i) out[i] = x[i - nn];
//       } else if(nn<0) {
//         int i = l, st = l+nn; 
//         while(i != st) out[--i] = fill;
//         for( ; i--; ) out[i] = x[i - nn];
//       } else return x;
//     } else {
//       if(l != g.size()) stop("length(x) must match length(g)");
//       if (Rf_isNull(t)) { // Ordered data
//         IntegerVector seen(ng);
//         if(nn>0) {
//           for(int i = 0; i != l; ++i) {  
//             if(seen[g[i]-1] == nn) {
//               out[i] = x[i-nn];
//             } else {
//               out[i] = fill;
//               ++seen[g[i]-1];
//             }
//           }
//         } else if(nn<0) {
//           for(int i = l; i--; ) { // good?? 
//             if(seen[g[i]-1] == nn) {
//               out[i] = x[i-nn];
//             } else {
//               out[i] = fill;
//               --seen[g[i]-1];
//             }
//           }
//         } else return x;
//       } else { // Unordered data: Timevar provided
//         IntegerVector ord = t; 
//         if(l != ord.size()) stop("length(x) must match length(t)");
//         if(ng != gs.size()) stop("ng must match length(gs)");
//         int **omap = new int*[ng];
//         for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//         for(int i = 0; i != l; ++i) omap[g[i]-1][ord[i]-1] = i;
//         if(nn>0) {
//           for(int i = 0; i != l; ++i) {
//             if(ord[i] > nn) {
//               out[i] = x[omap[g[i]-1][ord[i]-nn-1]]; // good ??? -> Yes, great and efficient !!!
//             } else {
//               out[i] = fill;
//             }
//           }
//         } else if(nn<0) {
//           IntegerVector tgt = gs+nn+1; // good ??? Best way ??
//           for(int i = 0; i != l; ++i) { // best loop ??
//             if(ord[i] < tgt[g[i]-1]) { // easier way, rather a counting vector by group ???
//               out[i] = x[omap[g[i]-1][ord[i]-nn-1]]; // good ??? -> Yes, great and efficient !!!
//             } else {
//               out[i] = fill;
//             }
//           }
//         } else return x;
//       }
//     }
//     return out;
//   } else { // Multiple lags !!!
//     NumericMatrix out = no_init_matrix(l, ns); // init better ??
//     if(ng == 0) { // No groups !!
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         if(np>0) {
//           int i = 0; 
//           while(i != np) outp[i++] = fill;
//           for( ; i != l; ++i) outp[i] = x[i - np];
//         } else if(np<0) {
//           int i = l, st = l+np; 
//           while(i != st) outp[--i] = fill;
//           for( ; i--; ) outp[i] = x[i - np];
//         } else outp = x;
//       }
//     } else {
//       if(l != g.size()) stop("length(x) must match length(g)");
//       if (Rf_isNull(t)) { // Ordered data
//         for(int p = ns; p--; ) {
//           int np = n[p];
//           NumericMatrix::Column outp = out( _ , p);
//           IntegerVector seenp(ng); // No init possible ??
//           if(np>0) {
//             for(int i = 0; i != l; ++i) { 
//               if(seenp[g[i]-1] == np) {
//                 outp[i] = x[i-np];
//               } else {
//                 outp[i] = fill;
//                 ++seenp[g[i]-1];
//               }
//             }
//           } else if(np<0) {
//             for(int i = l; i--; ) { // good?? 
//               if(seenp[g[i]-1] == np) {
//                 outp[i] = x[i-np];
//               } else {
//                 outp[i] = fill;
//                 --seenp[g[i]-1];
//               }
//             }
//           } else outp = x;
//         }
//       } else { // Unordered data: Timevar provided
//         IntegerVector ord = t; 
//         if(l != ord.size()) stop("length(x) must match length(t)");
//         if(ng != gs.size()) stop("ng must match length(gs)");
//         int **omap = new int*[ng];
//         for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//         for(int i = 0; i != l; ++i) omap[g[i]-1][ord[i]-1] = i;
//         for(int p = ns; p--; ) {
//           int np = n[p];
//           NumericMatrix::Column outp = out( _ , p);
//           if(np>0) {
//             for(int i = 0; i != l; ++i) {
//               if(ord[i] > np) {
//                 outp[i] = x[omap[g[i]-1][ord[i]-np-1]]; 
//               } else {
//                 outp[i] = fill;
//               }
//             }
//           } else if(np<0) {
//             IntegerVector tgt = gs+np+1; // good ??? Best way??, could do without this vector !!
//             for(int i = 0; i != l; ++i) { // best loop ??
//               if(ord[i] < tgt[g[i]-1]) {
//                 outp[i] = x[omap[g[i]-1][ord[i]-np-1]]; 
//               } else {
//                 outp[i] = fill;
//               }
//             }
//           } else outp = x;
//         }
//       }
//     }
//     return out;
//   }
// }


// Previous version: Multiple but only lags !!
// // [[Rcpp::export]]
// NumericVector flagCpp(NumericVector x, IntegerVector n = 1, double fill = NA_REAL, 
//                       int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, SEXP t = R_NilValue) { 
//   int l = x.size(), ns = n.size();
//   
//   if(ns == 1) { // Only one lag !! (Need extra method ??)
//     int nn = n[0]; // good ??
//     NumericVector out = no_init_vector(l); // init better ??
//     if(ng == 0) { // No groups !!
//         int i = 0; // best in front ??
//         while(i != nn) out[i++] = fill;
//         for( ; i != l; ++i) out[i] = x[i - nn];
//     } else {
//       if(l != g.size()) stop("length(x) must match length(g)");
//       if (Rf_isNull(t)) { // Ordered data
//         // warning("Panel-lag computed without timevar: Assuming ordered data"); // This takes time ??????? -> Do in main function !!
//         IntegerVector seen(ng);
//         for(int i = 0; i != l; ++i) { // This only works if groups are in order (sorted) !!!! 
//           if(seen[g[i]-1] == nn) {
//             out[i] = x[i-nn];
//           } else {
//             out[i] = fill;
//             ++seen[g[i]-1];
//           }
//         }
//       } else { // Unordered data: Timevar provided
//         IntegerVector ord = t; 
//         if(l != ord.size()) stop("length(x) must match length(t)");
//         if(ng != gs.size()) stop("ng must match length(gs)");
//         int **omap = new int*[ng];
//         for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//         for(int i = 0; i != l; ++i) omap[g[i]-1][ord[i]-1] = i;
//         for(int i = 0; i != l; ++i) {
//           if(ord[i] > nn) {
//             out[i] = x[omap[g[i]-1][ord[i]-nn-1]]; // good ??? -> Yes, great and efficient !!!
//           } else {
//             out[i] = fill;
//           }
//         }
//       }
//     }
//     return out;
//   } else { // Multiple lags !!!
//     NumericMatrix out = no_init_matrix(l, ns); // init better ??
//     if(ng == 0) { // No groups !!
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         int i = 0;
//         while(i != np) outp[i++] = fill;
//         for( ; i != l; ++i) outp[i] = x[i - np];
//       }
//     } else {
//       if(l != g.size()) stop("length(x) must match length(g)");
//       if (Rf_isNull(t)) { // Ordered data
//         for(int p = ns; p--; ) {
//           int np = n[p];
//           NumericMatrix::Column outp = out( _ , p);
//           IntegerVector seenp(ng); // No init possible ??
//           for(int i = 0; i != l; ++i) { 
//             if(seenp[g[i]-1] == np) {
//               outp[i] = x[i-np];
//             } else {
//               outp[i] = fill;
//               ++seenp[g[i]-1];
//             }
//           }
//         }
//       } else { // Unordered data: Timevar provided
//         IntegerVector ord = t; 
//         if(l != ord.size()) stop("length(x) must match length(t)");
//         if(ng != gs.size()) stop("ng must match length(gs)");
//         int **omap = new int*[ng];
//         for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//         for(int i = 0; i != l; ++i) omap[g[i]-1][ord[i]-1] = i;
//         for(int p = ns; p--; ) {
//           int np = n[p];
//           NumericMatrix::Column outp = out( _ , p);
//           for(int i = 0; i != l; ++i) {
//             if(ord[i] > np) {
//               out[i] = x[omap[g[i]-1][ord[i]-np-1]]; // good ??? -> Yes, great and efficient !!!
//             } else {
//               out[i] = fill;
//             }
//           }
//         }
//       }
//     }
//     return out;
//   }
// }

// // [[Rcpp::export]] // Better -> same speed as above and less code
// unc <- function(x) {res <- flag2Cpp(x); attributes(res) = NULL; res}
// NumericVector flag2Cpp(NumericVector x, IntegerVector n = 1, double fill = NA_REAL, 
//                        int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, SEXP t = R_NilValue) { 
//   int l = x.size(), ns = n.size();
//   NumericMatrix out = no_init_matrix(l, ns); // init better ??
//   if(ng == 0) { // No groups !!
//     for(int p = ns; p--; ) {
//       int np = n[p];
//       NumericMatrix::Column outp = out( _ , p);
//       int i = 0;
//       while(i != np) outp[i++] = fill;
//       for( ; i != l; ++i) outp[i] = x[i - np];
//     }
//   } else {
//     if(l != g.size()) stop("length(x) must match length(g)");
//     if (Rf_isNull(t)) { // Ordered data
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         IntegerVector seenp(ng); // No init possible ??
//         for(int i = 0; i != l; ++i) { 
//           if(seenp[g[i]-1] == np) {
//             outp[i] = x[i-np];
//           } else {
//             outp[i] = fill;
//             ++seenp[g[i]-1];
//           }
//         }
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t; 
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       if(ng != gs.size()) stop("ng must match length(gs)");
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//       for(int i = 0; i != l; ++i) omap[g[i]-1][ord[i]-1] = i;
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         NumericMatrix::Column outp = out( _ , p);
//         for(int i = 0; i != l; ++i) {
//           if(ord[i] > np) {
//             out[i] = x[omap[g[i]-1][ord[i]-np-1]]; // good ??? -> Yes, great and efficient !!!
//           } else {
//             out[i] = fill;
//           }
//         }
//       }
//     }
//   }
//   return out;
// }


// Previous Version: Only one lag
// // [[Rcpp::export]]
// NumericVector flagCpp(NumericVector x, int n = 1, double fill = NA_REAL, 
//                       int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, SEXP t = R_NilValue) { 
//   int l = x.size(), i = 0;
//   NumericVector out = no_init_vector(l); // init better ??
//   
//   if(ng == 0) { // No groups !!
//     // int i = 0; // best in front ??
//     while(i != n) out[i++] = fill;
//     for( ; i != l; ++i) out[i] = x[i - n];
//   } else {
//     if(l != g.size()) stop("length(x) must match length(g)");
//     if (Rf_isNull(t)) { // Ordered data
//     // warning("Panel-lag computed without timevar: Assuming ordered data"); // This takes time ??????? -> Do in main function !!
//     IntegerVector seen(ng);
//      for( ; i != l; ++i) { // This only works if groups are in order (sorted) !!!! 
//        if(seen[g[i]-1] == n) {
//            out[i] = x[i-n];
//        } else {
//            out[i] = fill;
//            ++seen[g[i]-1];
//        }
//      }
//     } else { // Unordered data: Timevar provided
//     IntegerVector ord = t; 
//     if(l != ord.size()) stop("length(x) must match length(t)");
//     if(ng != gs.size()) stop("ng must match length(gs)");
//     int **omap = new int*[ng];
//     for(     ; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//     for(i = 0; i != l; ++i) omap[g[i]-1][ord[i]-1] = i;
//     for(i = 0; i != l; ++i) {
//         if(ord[i] > n) {
//             out[i] = x[omap[g[i]-1][ord[i]-n-1]]; // good ??? -> Yes, great and efficient !!!
//         } else {
//             out[i] = fill;
//         }
//     }
//     }
//   }
//   return out;
//  }
