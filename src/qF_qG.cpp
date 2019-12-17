// // [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>
IntegerVector qFCppImpl(const Vector<RTYPE>& x, bool ordered = true, bool na_exclude = true) {
    Vector<RTYPE> levs = (ordered) ? sort_unique(x) : unique(x);
    IntegerVector out = no_init_vector(levs.size());
    if(na_exclude) {
      out = match(x, levs);
    } else {
      out = Rf_match(levs, x, NA_INTEGER); // Rf_matchE(levs, x, ls, wrap(Vector<RTYPE>::get_na()));
    }
    SHALLOW_DUPLICATE_ATTRIB(out, x); // works for all atomic objects ??
    if(TYPEOF(x) == STRSXP) { // slightly better
      out.attr("levels") = levs;
    } else {
      out.attr("levels") = Rf_coerceVector(levs, STRSXP);
    }
    out.attr("class") = (ordered && !na_exclude) ? CharacterVector::create("ordered","factor","na.included") :
                         ordered ? CharacterVector::create("ordered","factor") :
                         (!na_exclude) ? CharacterVector::create("factor","na.included") : "factor";
    return out;
}

template <int RTYPE>
IntegerVector qGCppImpl(const Vector<RTYPE>& x, bool ordered = true, bool na_exclude = true) {
    Vector<RTYPE> levs = (ordered) ? sort_unique(x) : unique(x);
    IntegerVector out = no_init_vector(levs.size());
    if(na_exclude) {
      out = match(x, levs);
    } else {
      out = Rf_match(levs, x, NA_INTEGER); // Rf_matchE(levs, x, ls, wrap(Vector<RTYPE>::get_na()));
    }
    // SHALLOW_DUPLICATE_ATTRIB(out, x); // needed ??
    out.attr("N.groups") = levs.size();
    out.attr("class") = (ordered && !na_exclude) ? CharacterVector::create("ordered","qG","na.included") :
                           ordered ? CharacterVector::create("ordered","qG") :
                        (!na_exclude) ? CharacterVector::create("qG","na.included") : "qG";
    return out;
}

// [[Rcpp::export]]   // do Cpp 11 solution using return macro ??
SEXP qFCpp(SEXP x, bool ordered = true, bool na_exclude = true) {
  switch(TYPEOF(x)) {
  case INTSXP: return qFCppImpl<INTSXP>(x, ordered, na_exclude);
  case REALSXP: return qFCppImpl<REALSXP>(x, ordered, na_exclude);
  case STRSXP: return qFCppImpl<STRSXP>(x, ordered, na_exclude);
  case LGLSXP: {
    LogicalVector xl = x;
    int l = xl.size();
    LogicalVector nd(3);
    IntegerVector out = no_init_vector(l);
    if(na_exclude) {
      for(int i = 0; i != l; ++i) {
        if(xl[i] == NA_LOGICAL) {
          out[i] = NA_INTEGER;
          nd[0] = true;
        } else if(xl[i]) {
          out[i] = 2;
          nd[2] = true;
        } else {
          out[i] = 1;
          nd[1] = true;
        }
      }
    } else {
      for(int i = 0; i != l; ++i) {
        if(xl[i] == NA_LOGICAL) {
          out[i] = 3;
          nd[0] = true;
        } else if(xl[i]) {
          out[i] = 2;
          nd[2] = true;
        } else {
          out[i] = 1;
          nd[1] = true;
        }
      }
    }
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    out.attr("levels") = CharacterVector::create("NA", "FALSE", "TRUE")[nd];
    out.attr("class") = (ordered && !na_exclude) ? CharacterVector::create("ordered","factor","na.included") :
                        ordered ? CharacterVector::create("ordered","factor") :
                        (!na_exclude) ? CharacterVector::create("factor","na.included") : "factor";
    return out;
  }
  default: stop("Not Supported SEXP Type");
  }
  return R_NilValue;
}

// [[Rcpp::export]]   // do Cpp 11 solution using return macro ??
SEXP qGCpp(SEXP x, bool ordered = true, bool na_exclude = true) {
  switch(TYPEOF(x)) {
  case INTSXP: return qGCppImpl<INTSXP>(x, ordered, na_exclude);
  case REALSXP: return qGCppImpl<REALSXP>(x, ordered, na_exclude);
  case STRSXP: return qGCppImpl<STRSXP>(x, ordered, na_exclude);
  case LGLSXP: {
    LogicalVector xl = x;
    int l = xl.size();
    LogicalVector nd(3);
    IntegerVector out = no_init_vector(l);
    if(na_exclude) {
      for(int i = 0; i != l; ++i) {
        if(xl[i] == NA_LOGICAL) {
          out[i] = NA_INTEGER;
          nd[0] = true;
        } else if(xl[i]) {
          out[i] = 2;
          nd[2] = true;
        } else {
          out[i] = 1;
          nd[1] = true;
        }
      }
    } else {
      for(int i = 0; i != l; ++i) {
        if(xl[i] == NA_LOGICAL) {
          out[i] = 3;
          nd[0] = true;
        } else if(xl[i]) {
          out[i] = 2;
          nd[2] = true;
        } else {
          out[i] = 1;
          nd[1] = true;
        }
      }
    }
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    out.attr("N.groups") = int(nd[0]+nd[1]+nd[2]);
    out.attr("class") = (ordered && !na_exclude) ? CharacterVector::create("ordered","qG","na.included") :
                         ordered ? CharacterVector::create("ordered","qG") :
                        (!na_exclude) ? CharacterVector::create("qG","na.included") : "qG";
    return out;
  }
  default: stop("Not Supported SEXP Type");
  }
  return R_NilValue;
}


template <int RTYPE>
Vector<RTYPE> funiqueImpl(const Vector<RTYPE>& x, bool ordered = true) {
  if(ordered) {
    Vector<RTYPE> out = sort_unique(x);
    DUPLICATE_ATTRIB(out, x);
    return out;
  } else {
    Vector<RTYPE> out = unique(x);
    DUPLICATE_ATTRIB(out, x);
    return out;
  }
}

// [[Rcpp::export]]
SEXP funique(SEXP x, bool ordered = true) {
  switch(TYPEOF(x)) {
  case INTSXP: return funiqueImpl<INTSXP>(x, ordered);
  case REALSXP: return funiqueImpl<REALSXP>(x, ordered);
  case STRSXP: return funiqueImpl<STRSXP>(x, ordered);
  case LGLSXP: {
    LogicalVector xl = x;
    LogicalVector nd(3);
    int ndc = 0;
    for(int i = xl.size(); i--; ) {
      if(nd[0] == false && xl[i] == NA_LOGICAL) {
        nd[0] = true;
        ++ndc;
      } else if(nd[2] == false && xl[i]) {
        nd[2] = true;
        ++ndc;
      } else if(nd[1] == false) {
        nd[1] = true;
        ++ndc;
      }
      if(ndc == 3) break;
    }
    LogicalVector out = LogicalVector::create(NA_LOGICAL, false, true)[nd];
    DUPLICATE_ATTRIB(out, x);
    return out;
  }
  default: stop("Not Supported SEXP Type");
  }
  return R_NilValue;
}

// template <int RTYPE>
// Vector<RTYPE> funiqueImpl(const Vector<RTYPE>& x , bool ordered = false) {
//   Vector<RTYPE> out = (ordered) ? sort_unique(x) : unique(x);
//   return out;
// }
//
// template <int RTYPE>
// IntegerVector qGCppImpl( const Vector<RTYPE>& x , bool ordered = false) {
//   Vector<RTYPE> levs = (ordered) ? sort_unique(x) : unique(x);
//   IntegerVector out = match(x, levs); // faster than just match ??
//   out.attr("N.groups") = levs.size();
//   out.attr("class") = "qG";
//   return out;
// }


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
