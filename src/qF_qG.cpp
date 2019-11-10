// // [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Old, non Cpp11 solution -> problem:: logical vectors not working !!!
template <int RTYPE>
IntegerVector qFCppImpl( const Vector<RTYPE>& x , bool ordered = false) {
    Vector<RTYPE> levs = (ordered) ? sort_unique(x) : unique(x); 
    IntegerVector out = match(x, levs);
    if(TYPEOF(x) == STRSXP) { // slightly better 
      out.attr("levels") = levs; 
    } else {
      out.attr("levels") = Rf_coerceVector(levs, STRSXP); 
    }
    // out.attr("levels") = Rf_coerceVector(levs, STRSXP); 
    out.attr("class") = (ordered) ? CharacterVector::create("ordered","factor") : "factor";
    return out;
}

template <int RTYPE>
IntegerVector qGCppImpl( const Vector<RTYPE>& x , bool ordered = false) {
    Vector<RTYPE> levs = (ordered) ? sort_unique(x) : unique(x); 
    IntegerVector out = match(x, levs); // faster than just match ??
    out.attr("N.groups") = levs.size();
    out.attr("class") = "qG";
    return out;
}

// [[Rcpp::export]]   // do Cpp 11 solution using return macro ?? 
SEXP qFCpp( SEXP x , bool ordered = false) {
  switch( TYPEOF(x) ) {
  case INTSXP: return qFCppImpl<INTSXP>(x, ordered);
  case REALSXP: return qFCppImpl<REALSXP>(x, ordered);
  case STRSXP: return qFCppImpl<STRSXP>(x, ordered);
  // case LGLSXP: return qFCppImpl<LGLSXP>(x, ordered);
  default: stop("Not Supported SEXP Type");
  }
  return R_NilValue;
}

// [[Rcpp::export]]   // do Cpp 11 solution using return macro ?? 
SEXP qGCpp( SEXP x , bool ordered = false) {
  switch( TYPEOF(x) ) {
  case INTSXP: return qGCppImpl<INTSXP>(x, ordered);
  case REALSXP: return qGCppImpl<REALSXP>(x, ordered);
  case STRSXP: return qGCppImpl<STRSXP>(x, ordered);
  // case LGLSXP: return qGCppImpl<LGLSXP>(x, ordered);
  default: stop("Not Supported SEXP Type");
  }
  return R_NilValue;
}

// Cpp 11 Solution: Not working for some reason !!
// template <int RTYPE>
// IntegerVector qFCppImpl( Vector<RTYPE> x , bool ordered) {
//   if(ordered) {
//     Vector<RTYPE> levs = sort_unique(x); // (ordered) ? sort_unique(x) : unique(x);
//     IntegerVector out = match(x, levs);
//     SET_ATTRIB(out, List::create(Rf_coerceVector(levs, STRSXP), CharacterVector::create("ordered","factor")));
//     // out.attr("levels") = Rf_coerceVector(levs, STRSXP); // as<CharacterVector>(levs); // -> same speed on large data
//     // out.attr("class") = (ordered) ? CharacterVector::create("ordered","factor") : "factor";
//     return out;
//   } else {
//     Vector<RTYPE> levs = unique(x); // (ordered) ? sort_unique(x) : unique(x);
//     IntegerVector out = match(x, levs);
//     SET_ATTRIB(out, List::create(Rf_coerceVector(levs, STRSXP), CharacterVector::create("factor")));
//     return out;
//   }
// }
// 
// 
// template <>
// IntegerVector qFCppImpl( Vector<CPLXSXP> x , bool ordered) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// IntegerVector qFCppImpl( Vector<VECSXP> x , bool ordered) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// IntegerVector qFCppImpl( Vector<RAWSXP> x , bool ordered) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// IntegerVector qFCppImpl( Vector<EXPRSXP> x , bool ordered) {
//   stop("Not supported SEXP type!");
// }
// 
// // [[Rcpp::export]]
// IntegerVector qFCpp(SEXP x , bool ordered = false) { // need const ???
//   RCPP_RETURN_VECTOR(qFCppImpl, x, ordered);
// }
// 
