// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Implemented smarter copy names !!

template <int RTYPE>
inline bool isnaNUM(typename Rcpp::traits::storage_type<RTYPE>::type x) { 
  return x != x;
}

template <int RTYPE>
inline bool isnaOTH(typename Rcpp::traits::storage_type<RTYPE>::type x) { 
  return x == Vector<RTYPE>::get_na();
}

template <int RTYPE>
Vector<RTYPE> flastCppImpl(const Vector<RTYPE>& x, int ng, const IntegerVector& g, bool narm) {
  int l = x.size();
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;
  if (ng == 0) {
    if(narm) {
      int j = l-1;
      auto last = x[j]; 
      while(isnanT(last) && j!=0) last = x[--j]; 
      Vector<RTYPE> out(1, last); // faster using create ??
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        out.attr("names") = CharacterVector::create(names[j]);
      }
      return out;
    } else {
      Vector<RTYPE> out(1, x[l-1]);
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        out.attr("names") = CharacterVector::create(names[l-1]);
      }
      return out;
    }
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    int ngs = 0;
    Vector<RTYPE> last = no_init_vector(ng); 
    DUPLICATE_ATTRIB(last, x);
    if(narm) {
      std::fill(last.begin(), last.end(), Vector<RTYPE>::get_na());
      if(Rf_getAttrib(x, R_NamesSymbol) == R_NilValue) {
        for(int i = l; i--; ) {
          if(!isnanT(x[i])) { 
            if(isnanT(last[g[i]-1])) { 
              last[g[i]-1] = x[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
      } else {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        if(names.size() != l) stop("x has a names attribute of length != length(x)");
        CharacterVector newnames = no_init_vector(ng);
        for(int i = l; i--; ) {
          if(!isnanT(x[i])) { 
            if(isnanT(last[g[i]-1])) { 
              last[g[i]-1] = x[i];
              newnames[g[i]-1] = names[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
        last.attr("names") = newnames;
      }
    } else {
      LogicalVector gl(ng, true); // std::vector<bool> glj(ng, true);? -> Nope, not faster (see matrix method)
      if(Rf_getAttrib(x, R_NamesSymbol) == R_NilValue) {
        for(int i = l; i--; ) {
          if(gl[g[i]-1]) {
            gl[g[i]-1] = false;
            last[g[i]-1] = x[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
      } else {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        if(names.size() != l) stop("x has a names attribute of length != length(x)");
        CharacterVector newnames = no_init_vector(ng);
        for(int i = l; i--; ) {
          if(gl[g[i]-1]) {
            gl[g[i]-1] = false;
            last[g[i]-1] = x[i];
            newnames[g[i]-1] = names[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
        last.attr("names") = newnames;
      }
    }
    return last;
  }
}

template <>
Vector<CPLXSXP> flastCppImpl(const Vector<CPLXSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<VECSXP> flastCppImpl(const Vector<VECSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<RAWSXP> flastCppImpl(const Vector<RAWSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<EXPRSXP> flastCppImpl(const Vector<EXPRSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP flastCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true){
  RCPP_RETURN_VECTOR(flastCppImpl, x, ng, g, narm);
}



// Previous Version: before isnan template:
// 
// #define isnan2(x) ((x) != (x)) // http://www.ebyte.it/library/codesnippets/WritingCppMacros.html
// 
// template <int RTYPE>
// Vector<RTYPE> flastCppImpl(Vector<RTYPE> x, int ng, IntegerVector g, bool narm) {
//   
//   int l = x.size();
//   if (ng == 0) {
//     if(narm) {
//       int j = l-1;
//       auto last = x[j]; 
//       if(TYPEOF(x) == REALSXP) {
//         while(last != last && j!=0) last = x[--j]; 
//       } else {
//         while(last == Vector<RTYPE>::get_na() && j!=0) last = x[--j]; 
//       }
//       Vector<RTYPE> out(1, last); // faster using create ??
//       DUPLICATE_ATTRIB(out, x);
//       return out;
//     } else {
//       Vector<RTYPE> out(1, x[l-1]);
//       DUPLICATE_ATTRIB(out, x);
//       return out;
//     }
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     int ngs = 0;
//     Vector<RTYPE> last = no_init_vector(ng); 
//     if(narm) {
//       std::fill(last.begin(), last.end(), Vector<RTYPE>::get_na());
//       if(TYPEOF(x) == REALSXP) { 
//         for(int i = l; i--; ) {
//           if(!isnan2(x[i])) { // use macro ??? -> Yes, not bad !! // 
//             if(isnan2(last[g[i]-1])) { // isnan
//               last[g[i]-1] = x[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }
//         }
//       } else {
//         for(int i = l; i--; ) {
//           if(x[i] != Vector<RTYPE>::get_na()) { 
//             if(last[g[i]-1] == Vector<RTYPE>::get_na()) { 
//               last[g[i]-1] = x[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }
//         }
//       }
//     } else {
//       LogicalVector gl(ng, true); 
//       for(int i = l; i--; ) {
//         if(gl[g[i]-1]) {
//           gl[g[i]-1] = false;
//           last[g[i]-1] = x[i];
//           ++ngs;
//           if(ngs == ng) break;
//         }
//       }
//     }
//     DUPLICATE_ATTRIB(last, x);
//     return last;
//   }
// }
// 
// template <>
// Vector<CPLXSXP> flastCppImpl(Vector<CPLXSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// Vector<VECSXP> flastCppImpl(Vector<VECSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// Vector<RAWSXP> flastCppImpl(Vector<RAWSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// Vector<EXPRSXP> flastCppImpl(Vector<EXPRSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// // [[Rcpp::export]]
// SEXP flastCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true){
//   RCPP_RETURN_VECTOR(flastCppImpl, x, ng, g, narm);
// }



// Previous Version: Before efficient attribute copy !!
// template <int RTYPE>
// Vector<RTYPE> flastCppImpl(Vector<RTYPE> x, int ng, IntegerVector g, bool narm) {
//   
//   int l = x.size();
//   if (ng == 0) {
//     if(narm) {
//       int j = l-1;
//       auto last = x[j]; 
//       if(TYPEOF(x) == REALSXP) {
//         while(last != last && j!=0) last = x[--j]; 
//       } else {
//         while(last == Vector<RTYPE>::get_na() && j!=0) last = x[--j]; 
//       }
//       if(x.hasAttribute("class")) {
//         Vector<RTYPE> out(1, last);
//         CharacterVector an = wrap(x.attributeNames());
//         for(int i = 0; i != an.size(); ++i) {
//           String s(an[i]);
//           out.attr(s) = x.attr(s);
//         }
//         return out;
//       } else {
//         return Vector<RTYPE>::create(last);
//       }
//     } else {
//       if(x.hasAttribute("class")) {
//         Vector<RTYPE> out(1, x[l-1]);
//         CharacterVector an = wrap(x.attributeNames());
//         for(int i = 0; i != an.size(); ++i) {
//           String s(an[i]);
//           out.attr(s) = x.attr(s);
//         }
//         return out;
//       } else {
//         return Vector<RTYPE>::create(x[l-1]);
//       }
//     }
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     int ngs = 0;
//     Vector<RTYPE> last = no_init_vector(ng); 
//     if(narm) {
//       std::fill(last.begin(), last.end(), Vector<RTYPE>::get_na());
//       if(TYPEOF(x) == REALSXP) { 
//         for(int i = l; i--; ) {
//           if(x[i] == x[i]) { // !isnan
//             if(last[g[i]-1] != last[g[i]-1]) { // isnan
//               last[g[i]-1] = x[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }
//         }
//       } else {
//         for(int i = l; i--; ) {
//           if(x[i] != Vector<RTYPE>::get_na()) { 
//             if(last[g[i]-1] == Vector<RTYPE>::get_na()) { 
//               last[g[i]-1] = x[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }
//         }
//       }
//     } else {
//       LogicalVector gl(ng, true); 
//       for(int i = l; i--; ) {
//         if(gl[g[i]-1]) {
//           gl[g[i]-1] = false;
//           last[g[i]-1] = x[i];
//           ++ngs;
//           if(ngs == ng) break;
//         }
//       }
//     }
//     if(x.hasAttribute("class")) {
//       CharacterVector an = wrap(x.attributeNames());
//       for(int i = 0; i != an.size(); ++i) {
//         String s(an[i]);
//         last.attr(s) = x.attr(s);
//       }
//     } 
//     return last;
//   }
// }
// 
// template <>
// Vector<CPLXSXP> flastCppImpl(Vector<CPLXSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// Vector<VECSXP> flastCppImpl(Vector<VECSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// Vector<RAWSXP> flastCppImpl(Vector<RAWSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// Vector<EXPRSXP> flastCppImpl(Vector<EXPRSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// // [[Rcpp::export]]
// SEXP flastCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true){
//   RCPP_RETURN_VECTOR(flastCppImpl, x, ng, g, narm);
// }

// Only Numeric Method
// // [[Rcpp::export]]
// NumericVector flastCpp(NumericVector x, int ng = 0, IntegerVector g = 0, 
//                         bool narm = true) {
//   int l = x.size();
//   if (ng == 0) {
//     if(narm) {
//       int j = l-1;
//       double last = x[j]; 
//       while(std::isnan(last) && j!=0) last = x[--j]; 
//       return NumericVector::create(last);
//     } else {
//       return NumericVector::create(x[l-1]);
//     }
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     int ngs = 0;
//     if(narm) {
//       NumericVector last(ng, NA_REAL); // Other way ??
//       for(int i = l; i--; ) {
//         if(!std::isnan(x[i])) {
//           if(std::isnan(last[g[i]-1])) {
//             last[g[i]-1] = x[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         } 
//       }
//       return last;
//     } else {
//       NumericVector last = no_init_vector(ng); // Other way ?? Stable ?? -> Yes !!!
//       LogicalVector gl(ng, true); // Other way around ?? -> Nope, this is faster !!
//       for(int i = l; i--; ) {
//         if(gl[g[i]-1]) {
//           gl[g[i]-1] = false;
//           last[g[i]-1] = x[i];
//           ++ngs;
//           if(ngs == ng) break;
//         }
//       }
//       return last;
//     }
//   }
// }
