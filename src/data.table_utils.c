#include "data.table.h"

bool isRealReallyInt(SEXP x) {
  if (!isReal(x)) return(false);
  R_xlen_t n=xlength(x), i=0;
  double *dx = REAL(x);
  while (i<n &&
         ( ISNA(dx[i]) ||
         ( R_FINITE(dx[i]) && dx[i] == (int)(dx[i])))) {
    i++;
  }
  return i==n;
}

SEXP isReallyReal(SEXP x) {
  SEXP ans = PROTECT(allocVector(INTSXP, 1));
  INTEGER(ans)[0] = 0;
  // return 0 (FALSE) when not type double, or is type double but contains integers
  // used to error if not passed type double but this needed extra is.double() calls in calling R code
  // which needed a repeat of the argument. Hence simpler and more robust to return 0 when not type double.
  if (isReal(x)) {
    int n=length(x), i=0;
    double *dx = REAL(x);
    while (i<n &&
           ( ISNA(dx[i]) ||
           ( R_FINITE(dx[i]) && dx[i] == (int)(dx[i])))) {
      i++;
    }
    if (i<n) INTEGER(ans)[0] = i+1;  // return the location of first element which is really real; i.e. not an integer
  }
  UNPROTECT(1);
  return(ans);
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
    for (int i=0; i<n; ++i)    if (xd[i]!=NA_INTEGER) {
      return false;
    }
    return true;
  }
  case REALSXP:
    if (Rinherits(x,char_integer64)) {
      const int64_t *xd = (int64_t *)REAL(x);
      for (int i=0; i<n; ++i)  if (xd[i]!=NA_INTEGER64) {
        return false;
      }
    } else {
      const double *xd = REAL(x);
      for (int i=0; i<n; ++i)  if (!ISNAN(xd[i])) {
        return false;
      }
    }
    return true;
  case STRSXP: {
    const SEXP *xd = STRING_PTR(x);
    for (int i=0; i<n; ++i)    if (xd[i]!=NA_STRING) {
      return false;
    }
    return true;
  }}
  if (!errorForBadType) return false;
  error("Unsupported type '%s' passed to allNA()", type2char(TYPEOF(x)));  // e.g. VECSXP; tests 2116.16-18
  // turned off allNA list support for now to avoid accidentally using it internally where we did not intend; allNA not yet exported
  //   https://github.com/Rdatatable/data.table/pull/3909#discussion_r329065950
}

/* colnamesInt
 * for provided data.table (or a list-like) and a subset of its columns, it returns integer positions of those columns in DT
 * handle columns input as: integer, double, character and NULL (handled as seq_along(x))
 * adds validation for:
 *   correct range [1,ncol], and if type real checks whole integer
 *   existing columns for character
 *   optionally check for no duplicates
 */

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
    for (int i=0; i<LENGTH(klass); i++) {
      if (STRING_ELT(klass, i) == char_) return true;
    }
  }
  return false;
}

bool Rinherits(SEXP x, SEXP char_) {
  // motivation was nanotime which is S4 and inherits from integer64 via S3 extends
  // R's C API inherits() does not cover S4 and returns FALSE for nanotime, as does our own INHERITS above.
  // R's R-level inherits() calls objects.c:inherits2 which calls attrib.c:R_data_class2 and
  // then attrib.c:S4_extends which itself calls R level methods:::.extendsForS3 which then calls R level methods::extends.
  // Since that chain of calls is so complicated and involves evaluating R level anyway, let's just reuse it.
  // Rinherits prefix with 'R' to signify i) it may call R level and is therefore not thread safe, and ii) includes R level inherits which covers S4.
  bool ans = INHERITS(x, char_);        // try standard S3 class character vector first
  if (!ans && char_==char_integer64)    // save the eval() for known S4 classes that inherit from integer64
    ans = INHERITS(x, char_nanotime);   // comment this out to test the eval() works for nanotime
  if (!ans && IS_S4_OBJECT(x)) {        // if it's not S4 we can save the overhead of R eval()
    SEXP vec = PROTECT(ScalarString(char_));           // TODO: cover this branch by making two new test S4 classes: one that
    SEXP call = PROTECT(lang3(sym_inherits, x, vec));  //       does inherit from integer64 and one that doesn't
    ans = LOGICAL(eval(call, R_GlobalEnv))[0]==1;
    UNPROTECT(2);
  }
  return ans;
}

SEXP copyAsPlain(SEXP x) {
  // v1.12.2 and before used standard R duplicate() to do this. But that's not guaranteed to not return an ALTREP.
  // e.g. ALTREP 'wrapper' on factor column (with materialized INTSXP) in package VIM under example(hotdeck)
  //      .Internal(inspect(x[[5]]))
  //      @558adf4d9508 13 INTSXP g0c0 [OBJ,NAM(7),ATT]  wrapper [srt=-2147483648,no_na=0]
  // 'AsPlain' is intended to convey unALTREP-ing; i.e. materializing and removing any ALTREP attributes too
  // For non-ALTREP this should do the same as R's duplicate(); but doesn't quite currently, so has to divert to duplicated() for now
  // Intended for use on columns; to either un-ALTREP them or duplicate shared memory columns; see copySharedColumns() below
  // Not intended to be called on a DT VECSXP where a concept of 'deep' might refer to whether the columns are copied

  if (!ALTREP(x)) return duplicate(x);
  // would prefer not to have this line, but without it test 1639.064 fails :
  //   Running test id 1639.064      Error in `[.data.table`(r, -ii) :
  //   Item 2 of i is -1 and item 1 is NA. Cannot mix negatives and NA.
  //   Calls: test.data.table ... FUN -> make.levels -> rbindlist -> [ -> [.data.table
  // Perhaps related to row names and the copyMostAttrib() below is not quite sufficient

  size_t n = XLENGTH(x);
  SEXP ans = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
  switch (TYPEOF(ans)) {
  case RAWSXP:
    memcpy(RAW(ans),     RAW(x),     n*sizeof(Rbyte));           // # nocov; add coverage when ALTREP is turned on for all types
    break;                                                       // # nocov
  case LGLSXP:
    memcpy(LOGICAL(ans), LOGICAL(x), n*sizeof(Rboolean));        // # nocov
    break;                                                       // # nocov
  case INTSXP:
    memcpy(INTEGER(ans), INTEGER(x), n*sizeof(int));             // covered by 10:1 after test 178
    break;
  case REALSXP:
    memcpy(REAL(ans),    REAL(x),    n*sizeof(double));          // covered by as.Date("2013-01-01")+seq(1,1000,by=10) after test 1075
    break;
  case CPLXSXP:
    memcpy(COMPLEX(ans), COMPLEX(x), n*sizeof(Rcomplex));        // # nocov
    break;                                                       // # nocov
  case STRSXP: {
    const SEXP *xp=STRING_PTR(x);                                // covered by as.character(as.hexmode(1:500)) after test 642
    for (R_xlen_t i=0; i<n; ++i) SET_STRING_ELT(ans, i, xp[i]);
  } break;
  case VECSXP: {
    const SEXP *xp=VECTOR_PTR(x);                                // # nocov
    for (R_xlen_t i=0; i<n; ++i) SET_VECTOR_ELT(ans, i, xp[i]);  // # nocov
  } break;                                                       // # nocov
  default:
    error("Internal error: unsupported type '%s' passed to copyAsPlain()", type2char(TYPEOF(x))); // # nocov
  }
  copyMostAttrib(x, ans); // e.g. factor levels, class etc, but not names, dim or dimnames
  if (ALTREP(ans))
    error("Internal error: type '%s' passed to copyAsPlain() but it seems copyMostAttrib() retains ALTREP attributes", type2char(TYPEOF(x))); // # nocov
  UNPROTECT(1);
  return ans;
}

