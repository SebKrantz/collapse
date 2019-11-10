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
Vector<RTYPE> ffirstCppImpl(const Vector<RTYPE>& x, int ng, const IntegerVector& g, bool narm) {
  int l = x.size(), end = l-1;
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;
  if (ng == 0) {
    if(narm) {
      int j = 0;
      auto first = x[j]; 
      while(isnanT(first) && j!=end) first = x[++j]; 
      Vector<RTYPE> out(1, first); // faster using create ??
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        out.attr("names") = CharacterVector::create(names[j]);
      }
      return out;
    } else {
      Vector<RTYPE> out(1, x[0]);
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        out.attr("names") = CharacterVector::create(names[0]);
      }
      return out;
    }
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    int ngs = 0;
    Vector<RTYPE> first = no_init_vector(ng); 
    DUPLICATE_ATTRIB(first, x);
    if(narm) {
      std::fill(first.begin(), first.end(), Vector<RTYPE>::get_na());
      if(Rf_getAttrib(x, R_NamesSymbol) == R_NilValue) {
        for(int i = 0; i != l; ++i) {
          if(!isnanT(x[i])) { 
            if(isnanT(first[g[i]-1])) { 
              first[g[i]-1] = x[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
      } else {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        if(names.size() != l) stop("x has a names attribute of length != length(x)");
        CharacterVector newnames = no_init_vector(ng);
        for(int i = 0; i != l; ++i) {
          if(!isnanT(x[i])) { 
            if(isnanT(first[g[i]-1])) { 
              first[g[i]-1] = x[i];
              newnames[g[i]-1] = names[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
        first.attr("names") = newnames;
      }
    } else {
      LogicalVector gl(ng, true); // std::vector<bool> glj(ng, true);? -> Nope, not faster (see matrix method)
      if(Rf_getAttrib(x, R_NamesSymbol) == R_NilValue) {
        for(int i = 0; i != l; ++i) {
          if(gl[g[i]-1]) {
            gl[g[i]-1] = false;
            first[g[i]-1] = x[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
      } else {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        if(names.size() != l) stop("x has a names attribute of length != length(x)");
        CharacterVector newnames = no_init_vector(ng);
        for(int i = 0; i != l; ++i) {
          if(gl[g[i]-1]) {
            gl[g[i]-1] = false;
            first[g[i]-1] = x[i];
            newnames[g[i]-1] = names[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
        first.attr("names") = newnames;
      }
    }
    return first;
  }
}

template <>
Vector<CPLXSXP> ffirstCppImpl(const Vector<CPLXSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<VECSXP> ffirstCppImpl(const Vector<VECSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<RAWSXP> ffirstCppImpl(const Vector<RAWSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<EXPRSXP> ffirstCppImpl(const Vector<EXPRSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}


// [[Rcpp::export]]
SEXP ffirstCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true){
  RCPP_RETURN_VECTOR(ffirstCppImpl, x, ng, g, narm);
}


// Previous Version: before isnan template (see also a bit of trial and error below):

// #define isnan2(x) ((x) != (x)) // http://www.ebyte.it/library/codesnippets/WritingCppMacros.html

// template <int RTYPE>
// Vector<RTYPE> ffirstCppImpl(Vector<RTYPE> x, int ng, IntegerVector g, bool narm) {
//   
//   int l = x.size();
//   if (ng == 0) {
//     if(narm) {
//       int j = 0;
//       auto first = x[j]; 
//       if(TYPEOF(x) == REALSXP) {
//         while(first != first && j!=l-1) first = x[++j]; 
//       } else {
//         while(first == Vector<RTYPE>::get_na() && j!=l-1) first = x[++j]; 
//       }
//       Vector<RTYPE> out(1, first); // faster using create ??
//       DUPLICATE_ATTRIB(out, x);
//       return out;
//     } else {
//       Vector<RTYPE> out(1, x[0]);
//       DUPLICATE_ATTRIB(out, x);
//       return out;
//     }
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     int ngs = 0;
//     Vector<RTYPE> first = no_init_vector(ng); 
//     if(narm) {
//       std::fill(first.begin(), first.end(), Vector<RTYPE>::get_na());
//       if(TYPEOF(x) == REALSXP) { 
//         for(int i = 0; i != l; ++i) {
//           if(!isnan2(x[i])) { // use macro ??? -> Yes, not bad !!
//             if(isnan2(first[g[i]-1])) { 
//               first[g[i]-1] = x[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }
//         }
//       } else {
//         for(int i = 0; i != l; ++i) {
//           if(x[i] != Vector<RTYPE>::get_na()) { 
//             if(first[g[i]-1] == Vector<RTYPE>::get_na()) { 
//               first[g[i]-1] = x[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }
//         }
//       }
//     } else {
//       LogicalVector gl(ng, true); // std::vector<bool> glj(ng, true);?
//       for(int i = 0; i != l; ++i) {
//         if(gl[g[i]-1]) {
//           gl[g[i]-1] = false;
//           first[g[i]-1] = x[i];
//           ++ngs;
//           if(ngs == ng) break;
//         }
//       }
//     }
//     DUPLICATE_ATTRIB(first, x);
//     return first;
//   }
// }




// Playing around: trying to create an isnana macro. 
// template <int RTYPE>
// bool isnan2Impl(typename traits::storage_type<RTYPE>::type x) { // https://en.cppreference.com/w/cpp/numeric/math/isnan
//   return x != x;
// }
// // Not working properly !!!
// // [[Rcpp::export]]
// bool isnan2(SEXP x){
//   RCPP_RETURN_VECTOR(isnan2Impl, x);
// }
// inline bool isnan2(SEXP x) { // https://en.cppreference.com/w/cpp/numeric/math/isnan
//   switch(TYPEOF(x)) {
//    case REALSXP: return as<double>(x) != as<double>(x);
//    case INTSXP: return as<int>(x) == NA_INTEGER;
//    case STRSXP: return as<String>(x) == NA_INTEGER;
//    case LGLSXP: return as<bool>(x) == NA_INTEGER;
//    default: stop("Not supported SEXP type");
//   }
// }


// Previous Version: Before efficient copy of attributes
// template <int RTYPE>
// Vector<RTYPE> ffirstCppImpl(Vector<RTYPE> x, int ng, IntegerVector g, bool narm) {
// 
//   int l = x.size();
//   if (ng == 0) {
//     if(narm) {
//       int j = 0;
//       // typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t; // or auto in cpp11 // https://stackoverflow.com/questions/42363850/template-for-class-returning-elements-of-rcppvectorrtype/42373077#42373077
//       auto first = x[j]; // auto is same speed as above
//       if(TYPEOF(x) == REALSXP) {
//         while(first != first && j!=l-1) first = x[++j]; // or Vector<RTYPE>::is_na( ); ??? -> very slow !!! first != first is even faster than std::isnan !!!
//       } else {
//         while(first == Vector<RTYPE>::get_na() && j!=l-1) first = x[++j]; // fast ?? 2step ?? 
//       }
//       if(x.hasAttribute("class")) {
//         Vector<RTYPE> out(1, first);
//         CharacterVector an = wrap(x.attributeNames());
//         for(int i = 0; i != an.size(); ++i) {
//           String s(an[i]);
//           out.attr(s) = x.attr(s);
//         }
//         return out;
//       } else {
//         return Vector<RTYPE>::create(first);
//       }
//     } else {
//       if(x.hasAttribute("class")) {
//         Vector<RTYPE> out(1, x[0]);
//         CharacterVector an = wrap(x.attributeNames());
//         for(int i = 0; i != an.size(); ++i) {
//           String s(an[i]);
//           out.attr(s) = x.attr(s);
//         }
//         return out;
//       } else {
//         return Vector<RTYPE>::create(x[0]);
//       }
//     }
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     int ngs = 0;
//     Vector<RTYPE> first = no_init_vector(ng); // Other way ?? Stable ?? -> Yes !!!
//     if(narm) {
//       std::fill(first.begin(), first.end(), Vector<RTYPE>::get_na());
//       // Vector<RTYPE> first(ng, Vector<RTYPE>::get_na()); // Other way ??  Na_Proxy // https://stackoverflow.com/questions/23744114/return-na-via-rcpp
//       if(TYPEOF(x) == REALSXP) { // Here iterator method is not faster !!
//         for(int i = 0; i != l; ++i) {
//           if(x[i] == x[i]) { // !isnan
//             if(first[g[i]-1] != first[g[i]-1]) { // isnan
//               first[g[i]-1] = x[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }
//         }
//       } else {
//         for(int i = 0; i != l; ++i) {
//           if(x[i] != Vector<RTYPE>::get_na()) { // 2step faster ??
//             if(first[g[i]-1] == Vector<RTYPE>::get_na()) { 
//               first[g[i]-1] = x[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }
//         }
//       }
//     } else {
//       LogicalVector gl(ng, true); // Other way around ?? -> Nope, this is faster !!
//       for(int i = 0; i != l; ++i) {
//         if(gl[g[i]-1]) {
//           gl[g[i]-1] = false;
//           first[g[i]-1] = x[i];
//           ++ngs;
//           if(ngs == ng) break;
//         }
//       }
//     }
//     if(x.hasAttribute("class")) {
//       CharacterVector an = wrap(x.attributeNames());
//       for(int i = 0; i != an.size(); ++i) {
//         String s(an[i]);
//         first.attr(s) = x.attr(s);
//       }
//     } 
//     return first;
//   }
// }
// 
// // These are the types not supported by std::isnan:
// 
// template <>
// Vector<CPLXSXP> ffirstCppImpl(Vector<CPLXSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// // template <>
// // Vector<STRSXP> ffirstCppImpl(Vector<STRSXP> x, int ng, IntegerVector, bool narm) {
// //   stop("Not supported SEXP type!");
// // }
// 
// template <>
// Vector<VECSXP> ffirstCppImpl(Vector<VECSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// Vector<RAWSXP> ffirstCppImpl(Vector<RAWSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// Vector<EXPRSXP> ffirstCppImpl(Vector<EXPRSXP> x, int ng, IntegerVector, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// 
// // [[Rcpp::export]]
// SEXP ffirstCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true){
//   RCPP_RETURN_VECTOR(ffirstCppImpl, x, ng, g, narm);
// }



// Only Numeric Vectors Version 
// // [[Rcpp::export]]
// NumericVector ffirstCpp(NumericVector x, int ng = 0, IntegerVector g = 0, 
//                        bool narm = true) {
//   int l = x.size();
//   if (ng == 0) {
//     if(narm) {
//       int j = 0;
//       double first = x[j]; 
//       while(std::isnan(first) && j!=l-1) first = x[++j]; 
//       return NumericVector::create(first);
//     } else {
//       return NumericVector::create(x[0]);
//     }
//   } else { // with groups
//     if(g.size() != l) stop("length(g) must match nrow(X)");
//     int ngs = 0;
//     if(narm) {
//       NumericVector first(ng, NA_REAL); // Other way ??
//       for(int i = 0; i != l; ++i) {
//         if(!std::isnan(x[i])) {
//           if(std::isnan(first[g[i]-1])) {
//             first[g[i]-1] = x[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         } 
//       }
//       return first;
//     } else {
//       NumericVector first = no_init_vector(ng); // Other way ?? Stable ?? -> Yes !!!
//       LogicalVector gl(ng, true); // Other way around ?? -> Nope, this is faster !!
//       for(int i = 0; i != l; ++i) {
//         if(gl[g[i]-1]) {
//           gl[g[i]-1] = false;
//           first[g[i]-1] = x[i];
//           ++ngs;
//           if(ngs == ng) break;
//         }
//       }
//       return first;
//     }
//   }
// }

// Test: Logical vectors are initialized to false !!!
// // [[Rcpp::export]]
// LogicalVector tf(int k = 3) {
//   LogicalVector g(k);
//   return g;
// }