/*
 This code is adapted from the data.table package: http://r-datatable.com
 and licensed under a Mozilla Public License 2.0 (MPL-2.0) license.
*/

#include "data.table.h"


int need2utf8(SEXP x) {
  const int xlen = length(x);
  const SEXP *xd = STRING_PTR_RO(x);
  // for (int i=0; i<xlen; i++) {
  //   if (NEED2UTF8(xd[i]))
  //     return(true);
  // }
  // return(false);
  if (xlen <= 1) return xlen == 1 ? NEED2UTF8(xd[0]) : 0;
  return NEED2UTF8(xd[0]) || NEED2UTF8(xd[xlen/2]) || NEED2UTF8(xd[xlen-1]);
}

SEXP coerceUtf8IfNeeded(SEXP x) {
  if (!need2utf8(x))
    return(x);
  const int xlen = length(x);
  SEXP ans = PROTECT(allocVector(STRSXP, xlen));
  const SEXP *xd = STRING_PTR_RO(x);
  for (int i=0; i<xlen; i++) {
    SET_STRING_ELT(ans, i, ENC2UTF8(xd[i]));
  }
  UNPROTECT(1);
  return(ans);
}


SEXP setnames(SEXP x, SEXP nam) {
  if(TYPEOF(nam) != STRSXP) error("names need to be character typed");
  if(INHERITS(x, char_datatable)) {
    int n = TRUELENGTH(x), l = LENGTH(nam);
    if(n < l) { // error("Invalid data.table (underallocated), use qDT(data) to make valid.");
      setAttrib(x, R_NamesSymbol, nam);
      // setselfref(x);
      return x;
    }
    SEXP newnam = PROTECT(allocVector(STRSXP, n)),
      *pnn = SEXPPTR(newnam), *pn = SEXPPTR(nam);
    for(int i = 0; i < l; ++i) pnn[i] = pn[i];
    SETLENGTH(newnam, l);
    SET_TRUELENGTH(newnam, n);
    setAttrib(x, R_NamesSymbol, newnam);
    setselfref(x);
    UNPROTECT(1);
  } else setAttrib(x, R_NamesSymbol, nam);
  return x;
}

bool allNA(SEXP x, bool errorForBadType) {
  // less space and time than all(is.na(x)) at R level because that creates full size is.na(x) first before all()
  // whereas this allNA can often return early on testing the first value without reading the rest
  const int n = length(x);
  if (n==0) // empty vectors (including raw(), NULL, and list()) same as R's all(is.na()) true result; tests 2116.*
    return true;
  switch (TYPEOF(x)) {
  case RAWSXP: // raw doesn't support NA so always false (other than length 0 case above)
    return false;
  case LGLSXP:
  case INTSXP: {
    const int *xd = INTEGER(x);
    for (int i=0; i != n; ++i)    if (xd[i]!=NA_INTEGER) {
      return false;
    }
    return true;
  }
  case REALSXP:
    if (INHERITS(x,char_integer64)) {
      const int64_t *xd = (int64_t *)REAL(x);
      for (int i=0; i != n; ++i)  if (xd[i]!=NA_INTEGER64) {
        return false;
      }
    } else {
      const double *xd = REAL(x);
      for (int i=0; i != n; ++i)  if (!ISNAN(xd[i])) {
        return false;
      }
    }
    return true;
  case STRSXP: {
    const SEXP *xd = SEXPPTR(x);
    for (int i=0; i != n; ++i)    if (xd[i]!=NA_STRING) {
      return false;
    }
    return true;
  }}
  if (!errorForBadType) return false;
  error("Unsupported type '%s' passed to allNA()", type2char(TYPEOF(x)));  // e.g. VECSXP; tests 2116.16-18
  // turned off allNA list support for now to avoid accidentally using it internally where we did not intend; allNA not yet exported
  //   https://github.com/Rdatatable/data.table/pull/3909#discussion_r329065950
}

SEXP allNAv(SEXP x, SEXP errorForBadType) {
  return ScalarLogical(allNA(x, asLogical(errorForBadType)));
}

inline bool INHERITS(SEXP x, SEXP char_) {
  // Thread safe inherits() by pre-calling install() in init.c and then
  // passing those char_* in here for simple and fast non-API pointer compare.
  // The thread-safety aspect here is only currently actually needed for list columns in
  // fwrite() where the class of the cell's vector is tested; the class of the column
  // itself is pre-stored by fwrite (for example in isInteger64[] and isITime[]).
  // Thread safe in the limited sense of correct and intended usage :
  // i) no API call such as install() or mkChar() must be passed in.
  // ii) no attrib writes must be possible in other threads.
  SEXP klass;
  if (isString(klass = getAttrib(x, R_ClassSymbol))) {
    for (int i=0; i<LENGTH(klass); ++i) {
      if (STRING_ELT(klass, i) == char_) return true;
    }
    if (char_==char_integer64) {
      // package:nanotime is S4 and inherits from integer64 via S3 extends; i.e. integer64 does not appear in its R_ClassSymbol
      // R's C API inherits() does not cover S4 and returns FALSE for nanotime
      // R's R-level inherits() calls objects.c:inherits2 which calls attrib.c:R_data_class2 and
      // then attrib.c:S4_extends which itself calls R level methods:::.extendsForS3 which then calls R level methods::extends.
      // Since that chain of calls is so complicated and involves evaluating R level (not thread-safe) we
      // special case nanotime here. We used to have Rinherits() as well which did call R level but couldn't be called from
      // parallel regions. That became too hard to reason about two functions, #4752.
      // If any other classes come to light that, like nanotime, S4 inherit from integer64, we can i) encourage them to change
      // to regular S3, or ii) state we simply don't support that; i.e. nanotime was an exception, or iii) add a function that
      // gets called on C entry points which loops through columns and if any are S4 calls the old Rinherits() to see if they S4
      // inherit from integer64, and if so add that class to a vector that gets looped through here. That way we isolate the
      // non-TS call into argument massage header code, and we can continue to use INHERITS() throughout the code base.
      for (int i=0; i<LENGTH(klass); ++i) {
        if (STRING_ELT(klass, i) == char_nanotime) return true;
      }
    }
  }
  return false;
}

// Enhanced version of the original
SEXP dt_na(SEXP x, SEXP cols, SEXP Rprop, SEXP Rcount) {
  int n = 0, elem, ncol = LENGTH(cols), count = asLogical(Rcount);
  double prop = asReal(Rprop);
  if(ISNAN(prop) || prop < 0.0 || prop > 1.0) error("prop needs to be a proportion [0, 1]");

  if(!isNewList(x)) error("Internal error. Argument 'x' to missing_cases is type '%s' not 'list'", type2char(TYPEOF(x))); // # nocov
  if(!isInteger(cols)) error("Internal error. Argument 'cols' to missing_cases is type '%s' not 'integer'", type2char(TYPEOF(cols))); // # nocov
  for (int i = 0; i < ncol; ++i) {
    elem = INTEGER(cols)[i];
    if(elem < 1 || elem > LENGTH(x))
      error("Item %d of 'cols' is %d which is outside 1-based range [1,ncol(x)=%d]", i+1, elem, LENGTH(x));
    if(!n) n = length(VECTOR_ELT(x, elem-1));
  }

  SEXP ans = PROTECT(allocVector(LGLSXP, n));
  int *ians = LOGICAL(ans);
  memset(ians, 0, sizeof(int) * n);  // for (int i=0; i != n; ++i) ians[i]=0;

  if(count || prop > 0.0) { // More than 1 missing row, or counting mising values
    // if(prop == 1) { // Not sensible: better skip lists...
    //   // Preliminary check for early return
    //   for (int i = 0, tv; i < ncol; ++i) {
    //     tv = TYPEOF(VECTOR_ELT(x, INTEGER(cols)[i]-1));
    //     if(tv != LGLSXP && tv != INTSXP && tv != REALSXP && tv != STRSXP && tv != CPLXSXP && tv != NILSXP) {
    //       UNPROTECT(1);
    //       return(ans);
    //     }
    //   }
    // }

    // Counting the missing values
    int len = ncol;
    for (int i = 0; i < ncol; ++i) {
      SEXP v = VECTOR_ELT(x, INTEGER(cols)[i]-1);
      if (!length(v) || isNewList(v) || isList(v) || TYPEOF(v) == RAWSXP) {
        --len; continue;
      }
      if (n != length(v))
        error("Column %d of input list x is length %d, inconsistent with first column of that item which is length %d.", i+1,length(v),n);
      switch (TYPEOF(v)) {
      case LGLSXP: {
        const int *iv = LOGICAL(v);
        for (int j=0; j != n; ++j) ians[j] += (iv[j] == NA_LOGICAL);
      } break;
      case INTSXP: {
        const int *iv = INTEGER(v);
        for (int j=0; j != n; ++j) ians[j] += (iv[j] == NA_INTEGER);
      } break;
      case STRSXP: {
        const SEXP *sv = SEXPPTR(v);
        for (int j=0; j != n; ++j) ians[j] += (sv[j] == NA_STRING);
      } break;
      case REALSXP: {
        const double *dv = REAL(v);
        if (INHERITS(v, char_integer64)) {
          for (int j=0; j != n; ++j) ians[j] += (dv[j] == NA_INT64_D);
        } else {
          for (int j=0; j != n; ++j) ians[j] += ISNAN(dv[j]);
        }
      } break;
      case CPLXSXP: {
        const Rcomplex *dv = COMPLEX(v);
        for (int j=0; j != n; ++j) ians[j] += (ISNAN(dv[j].r) || ISNAN(dv[j].i));
      } break;
      default:
        error("Unsupported column type '%s'", type2char(TYPEOF(v)));
      }
    }
    if(count) {
      SETTOF(ans, INTSXP);
    } else {
      // This computes the result
      if(prop < 1.0) {
        len = (int)((double)len * prop);
        if(len < 1) len = 1;
      }
      for (int j = 0; j != n; ++j) ians[j] = ians[j] >= len;
    }
  } else { // Any missing (default)
    for (int i = 0; i < ncol; ++i) {
      SEXP v = VECTOR_ELT(x, INTEGER(cols)[i]-1);
      if (!length(v) || isNewList(v) || isList(v)) continue; // like stats:::na.omit.data.frame, skip list/pairlist columns
      if (n != length(v))
        error("Column %d of input list x is length %d, inconsistent with first column of that item which is length %d.", i+1,length(v),n);
      switch (TYPEOF(v)) {
      case LGLSXP: {
        const int *iv = LOGICAL(v);
        for (int j=0; j != n; ++j) ians[j] |= (iv[j] == NA_LOGICAL);
      } break;
      case INTSXP: {
        const int *iv = INTEGER(v);
        for (int j=0; j != n; ++j) ians[j] |= (iv[j] == NA_INTEGER);
      } break;
      case STRSXP: {
        const SEXP *sv = SEXPPTR(v);
        for (int j=0; j != n; ++j) ians[j] |= (sv[j] == NA_STRING);
      } break;
      case REALSXP: {
        const double *dv = REAL(v);
        if (INHERITS(v, char_integer64)) {
          for (int j=0; j != n; ++j) ians[j] |= (dv[j] == NA_INT64_D);
        } else {
          for (int j=0; j != n; ++j) ians[j] |= ISNAN(dv[j]);
        }
      } break;
      case RAWSXP: {
        // no such thing as a raw NA
        // vector already initialised to all 0's
      } break;
      case CPLXSXP: {
        // taken from https://github.com/wch/r-source/blob/d75f39d532819ccc8251f93b8ab10d5b83aac89a/src/main/coerce.c
        const Rcomplex *dv = COMPLEX(v);
        for (int j=0; j != n; ++j) ians[j] |= (ISNAN(dv[j].r) || ISNAN(dv[j].i));
      } break;
      default:
        error("Unsupported column type '%s'", type2char(TYPEOF(v)));
      }
    }
  }

  UNPROTECT(1);
  return(ans);
}

// from data.table_frank.c -> simplified frank, only dense method !!

SEXP frankds(SEXP xorderArg, SEXP xstartArg, SEXP xlenArg, SEXP dns) {
  int i=0, j=0, k=0, end=0, n, ng;
  int *xstart = INTEGER(xstartArg), *xlen = INTEGER(xlenArg), *xorder = INTEGER(xorderArg);
  n = length(xorderArg);
  ng = length(xstartArg);
  SEXP ans = PROTECT(allocVector(INTSXP, n));
  int *ians = INTEGER(ans);
  if(n > 0) {
    switch(asInteger(dns)) {
    case 0: // Not Sorted
      k=1;
      for (i = 0; i != ng; i++) {
        for (j = xstart[i]-1, end = xstart[i]+xlen[i]-1; j < end; j++)
          ians[xorder[j]-1] = k;
        k++;
      }
      break;
    case 1: // Sorted
      k=1;
      for (i = 0; i != ng; i++) {
        for (j = xstart[i]-1, end = xstart[i]+xlen[i]-1; j < end; j++) ians[j] = k;
        k++;
      }
      break;
    case 2: // This is basically run-length type group-id
      for (i = 0; i != ng; i++) {
        k=1;
        for (j = xstart[i]-1, end = xstart[i]+xlen[i]-1; j < end; j++)
          ians[xorder[j]-1] = k++;
      }
      break;
    default: error("dns must be 0, 1 or 2");
    }
  }
  UNPROTECT(1);
  return(ans);
}

// from data.table_assign.c:
SEXP setcolorder(SEXP x, SEXP o) {
  SEXP names = getAttrib(x, R_NamesSymbol);
  const int *od = INTEGER(o), ncol=LENGTH(x);
  if (isNull(names)) error("list passed to setcolorder has no names");
  if (ncol != LENGTH(names))
    error("Internal error: dt passed to setcolorder has %d columns but %d names", ncol, LENGTH(names));  // # nocov

  // Double-check here at C level that o[] is a strict permutation of 1:ncol. Reordering columns by reference makes no
  // difference to generations/refcnt so we can write behind barrier in this very special case of strict permutation.
  bool *seen = R_Calloc(ncol, bool);
  for (int i=0; i != ncol; ++i) {
    if (od[i]==NA_INTEGER || od[i]<1 || od[i]>ncol)
      error("Internal error: o passed to Csetcolorder contains an NA or out-of-bounds");  // # nocov
    if (seen[od[i]-1])
      error("Internal error: o passed to Csetcolorder contains a duplicate");             // # nocov
    seen[od[i]-1] = true;
  }
  R_Free(seen);

  SEXP *tmp = R_Calloc(ncol, SEXP), *namesd = SEXPPTR(names);
  const SEXP *xd = SEXPPTR_RO(x);
  for (int i=0; i != ncol; ++i) tmp[i] = xd[od[i]-1];
  for (int i=0; i != ncol; ++i) SET_VECTOR_ELT(x, i, tmp[i]);
  // SEXP *xd = SEXPPTR(x);
  // for (int i=0; i != ncol; ++i) tmp[i] = xd[od[i]-1];
  // memcpy(xd, tmp, ncol*sizeof(SEXP)); // sizeof is type size_t so no overflow here
  for (int i=0; i != ncol; ++i) tmp[i] = namesd[od[i]-1];
  memcpy(namesd, tmp, ncol*sizeof(SEXP));
  // No need to change key (if any); sorted attribute is column names not positions
  R_Free(tmp);
  return(R_NilValue);
}




