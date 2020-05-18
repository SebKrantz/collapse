// // [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>
IntegerVector qFCppImpl(const Vector<RTYPE>& x, bool sort = true, bool ordered = true, bool na_exclude = true) {
    Vector<RTYPE> levs = (sort && na_exclude) ? na_omit(sort_unique(x)) :
                          sort ? sort_unique(x) :
                          na_exclude ? na_omit(unique(x)) : unique(x);
    IntegerVector out = (na_exclude || RTYPE != REALSXP) ? match(x, levs) : as<IntegerVector>(Rf_match(levs, x, NA_INTEGER));
    SHALLOW_DUPLICATE_ATTRIB(out, x); // works for all atomic objects ?
    if(RTYPE == STRSXP) {
      out.attr("levels") = levs;
    } else {
      out.attr("levels") = Rf_coerceVector(levs, STRSXP); // What about date objects...
    }
    out.attr("class") = (ordered && !na_exclude) ? CharacterVector::create("ordered","factor","na.included") :
                         ordered ? CharacterVector::create("ordered","factor") :
                        (!na_exclude) ? CharacterVector::create("factor","na.included") : "factor";
    return out;
}

template <int RTYPE>
IntegerVector qGCppImpl(const Vector<RTYPE>& x, bool sort = true, bool ordered = true, bool na_exclude = true, bool retgrp = false) {
  Vector<RTYPE> levs = (sort && na_exclude) ? na_omit(sort_unique(x)) :
                        sort ? sort_unique(x) :
                        na_exclude ? na_omit(unique(x)) : unique(x);
  IntegerVector out = (na_exclude || RTYPE != REALSXP) ? match(x, levs) : as<IntegerVector>(Rf_match(levs, x, NA_INTEGER));
  // SHALLOW_DUPLICATE_ATTRIB(out, x); // needed ? -> Nah, it's a programmers function..
  out.attr("N.groups") = levs.size();
  if(retgrp) {
    DUPLICATE_ATTRIB(levs, x);
    out.attr("groups") = levs;
  }
  out.attr("class") = (ordered && !na_exclude) ? CharacterVector::create("ordered","qG","na.included") :
                       ordered ? CharacterVector::create("ordered","qG") :
                      (!na_exclude) ? CharacterVector::create("qG","na.included") : "qG";
  return out;
}

// [[Rcpp::export]]   // do Cpp 11 solution using return macro ?
SEXP qFCpp(SEXP x, bool sort = true, bool ordered = true, bool na_exclude = true) {
  switch(TYPEOF(x)) {
  case INTSXP: return qFCppImpl<INTSXP>(x, sort, ordered, na_exclude);
  case REALSXP: return qFCppImpl<REALSXP>(x, sort, ordered, na_exclude);
  case STRSXP: return qFCppImpl<STRSXP>(x, sort, ordered, na_exclude);
  case LGLSXP: {
    LogicalVector xl = x;
    int l = xl.size();
    LogicalVector nd(3);
    IntegerVector out = no_init_vector(l);
    if(na_exclude) {
      for(int i = 0; i != l; ++i) {
        if(xl[i] == NA_LOGICAL) {
          out[i] = NA_INTEGER;
        } else if(xl[i]) {
          out[i] = 2;
          nd[1] = true;
        } else {
          out[i] = 1;
          nd[0] = true;
        }
      }
    } else {
      for(int i = 0; i != l; ++i) {
        if(xl[i] == NA_LOGICAL) {
          out[i] = 3;
          nd[2] = true;
        } else if(xl[i]) {
          out[i] = 2;
          nd[1] = true;
        } else {
          out[i] = 1;
          nd[0] = true;
        }
      }
    }
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    out.attr("levels") = CharacterVector::create("FALSE", "TRUE", NA_STRING)[nd];
    out.attr("class") = (ordered && !na_exclude) ? CharacterVector::create("ordered","factor","na.included") :
                        ordered ? CharacterVector::create("ordered","factor") :
                        (!na_exclude) ? CharacterVector::create("factor","na.included") : "factor";
    return out;
  }
  default: stop("Not Supported SEXP Type");
  }
  return R_NilValue;
}

// [[Rcpp::export]]   // do Cpp 11 solution using return macro ?
SEXP qGCpp(SEXP x, bool sort = true, bool ordered = true, bool na_exclude = true, bool retgrp = false) {
  switch(TYPEOF(x)) {
  case INTSXP: return qGCppImpl<INTSXP>(x, sort, ordered, na_exclude, retgrp);
  case REALSXP: return qGCppImpl<REALSXP>(x, sort, ordered, na_exclude, retgrp);
  case STRSXP: return qGCppImpl<STRSXP>(x, sort, ordered, na_exclude, retgrp);
  case LGLSXP: {
    LogicalVector xl = x;
    int l = xl.size();
    LogicalVector nd(3);
    IntegerVector out = no_init_vector(l);
    if(na_exclude) {
      for(int i = 0; i != l; ++i) {
        if(xl[i] == NA_LOGICAL) {
          out[i] = NA_INTEGER;
        } else if(xl[i]) {
          out[i] = 2;
          nd[1] = true;
        } else {
          out[i] = 1;
          nd[0] = true;
        }
      }
    } else {
      for(int i = 0; i != l; ++i) {
        if(xl[i] == NA_LOGICAL) {
          out[i] = 3;
          nd[2] = true;
        } else if(xl[i]) {
          out[i] = 2;
          nd[1] = true;
        } else {
          out[i] = 1;
          nd[0] = true;
        }
      }
    }
    // SHALLOW_DUPLICATE_ATTRIB(out, x);
    out.attr("N.groups") = int(nd[0]+nd[1]+nd[2]);
    if(retgrp) out.attr("groups") = CharacterVector::create("FALSE", "TRUE", NA_STRING)[nd];
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
      if(!nd[2] && xl[i] == NA_LOGICAL) {
        nd[2] = true;
        ++ndc;
      } else if(!nd[1] && xl[i]) {
        nd[1] = true;
        ++ndc;
      } else if(!nd[0]) {
        nd[0] = true;
        ++ndc;
      }
      if(ndc == 3) break;
    }
    LogicalVector out = LogicalVector::create(false, true, NA_LOGICAL)[nd];
    DUPLICATE_ATTRIB(out, x);
    return out;
  }
  default: stop("Not Supported SEXP Type");
  }
  return R_NilValue;
}

