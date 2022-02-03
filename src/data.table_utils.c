/*
 This code is adapted from the data.table package: http://r-datatable.com
 and licensed under a Mozilla Public License 2.0 (MPL-2.0) license.
*/

#include "data.table.h"
#include <Rdefines.h>

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
    const SEXP *xd = STRING_PTR(x);
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
    for (R_xlen_t i=0; i != n; ++i) SET_STRING_ELT(ans, i, xp[i]);
  } break;
  case VECSXP: {
    const SEXP *xp=SEXPPTR(x);                                // # nocov
    for (R_xlen_t i=0; i != n; ++i) SET_VECTOR_ELT(ans, i, xp[i]);  // # nocov
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

/* No longer needed in GRP.default ! -> group sizes are now directly calculated by radixsort!
SEXP uniqlengths(SEXP x, SEXP n) {
  // seems very similar to rbindlist.c:uniq_lengths. TODO: centralize into common function
  if (TYPEOF(x) != INTSXP) error("Input argument 'x' to 'uniqlengths' must be an integer vector");
  if (TYPEOF(n) != INTSXP || length(n) != 1) error("Input argument 'n' to 'uniqlengths' must be an integer vector of length 1");
  R_len_t len = length(x);
  SEXP ans = PROTECT(allocVector(INTSXP, len));
  for (R_len_t i=1; i<len; i++) {
    INTEGER(ans)[i-1] = INTEGER(x)[i] - INTEGER(x)[i-1];
  }
  if (len>0) INTEGER(ans)[len-1] = INTEGER(n)[0] - INTEGER(x)[len-1] + 1;
  UNPROTECT(1);
  return(ans);
}
 */

// from data.table_frank.c -> simplified frank, only dense method !!

SEXP dt_na(SEXP x, SEXP cols) {
  int n=0, elem;

  if (!isNewList(x)) error("Internal error. Argument 'x' to Cdt_na is type '%s' not 'list'", type2char(TYPEOF(x))); // # nocov
  if (!isInteger(cols)) error("Internal error. Argument 'cols' to Cdt_na is type '%s' not 'integer'", type2char(TYPEOF(cols))); // # nocov
  for (int i=0; i<LENGTH(cols); ++i) {
    elem = INTEGER(cols)[i];
    if (elem<1 || elem>LENGTH(x))
      error("Item %d of 'cols' is %d which is outside 1-based range [1,ncol(x)=%d]", i+1, elem, LENGTH(x));
    if (!n) n = length(VECTOR_ELT(x, elem-1));
  }
  SEXP ans = PROTECT(allocVector(LGLSXP, n));
  int *ians = LOGICAL(ans);
  for (int i=0; i != n; ++i) ians[i]=0;
  for (int i=0; i<LENGTH(cols); ++i) {
    SEXP v = VECTOR_ELT(x, INTEGER(cols)[i]-1);
    if (!length(v) || isNewList(v) || isList(v)) continue; // like stats:::na.omit.data.frame, skip list/pairlist columns
    if (n != length(v))
      error("Column %d of input list x is length %d, inconsistent with first column of that item which is length %d.", i+1,length(v),n);
    switch (TYPEOF(v)) {
    case LGLSXP: {
      const int *iv = LOGICAL(v);
      for (int j=0; j != n; ++j) ians[j] |= (iv[j] == NA_LOGICAL);
    }
      break;
    case INTSXP: {
      const int *iv = INTEGER(v);
      for (int j=0; j != n; ++j) ians[j] |= (iv[j] == NA_INTEGER);
    }
      break;
    case STRSXP: {
      const SEXP *sv = STRING_PTR(v);
      for (int j=0; j != n; ++j) ians[j] |= (sv[j] == NA_STRING);
    }
      break;
    case REALSXP: {
      const double *dv = REAL(v);
      if (INHERITS(v, char_integer64)) {
        for (int j=0; j != n; ++j) {
          ians[j] |= (DtoLL(dv[j]) == NA_INT64_LL);   // TODO: can be == NA_INT64_D directly
        }
      } else {
        for (int j=0; j != n; ++j) ians[j] |= ISNAN(dv[j]);
      }
    }
      break;
    case RAWSXP: {
      // no such thing as a raw NA
      // vector already initialised to all 0's
    }
      break;
    case CPLXSXP: {
      // taken from https://github.com/wch/r-source/blob/d75f39d532819ccc8251f93b8ab10d5b83aac89a/src/main/coerce.c
      for (int j=0; j != n; ++j) ians[j] |= (ISNAN(COMPLEX(v)[j].r) || ISNAN(COMPLEX(v)[j].i));
    }
      break;
    default:
      error("Unsupported column type '%s'", type2char(TYPEOF(v)));
    }
  }
  UNPROTECT(1);
  return(ans);
}

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

// extern SEXP char_integer64;

// internal version of anyNA for data.tables
// SEXP anyNA(SEXP x, SEXP cols) {
//   int i, j, n=0, elem;
//
//   if (!isNewList(x)) error("Internal error. Argument 'x' to CanyNA is type '%s' not 'list'", type2char(TYPEOF(x))); // #nocov
//   if (!isInteger(cols)) error("Internal error. Argument 'cols' to CanyNA is type '%s' not 'integer'", type2char(TYPEOF(cols))); // # nocov
//   for (i=0; i<LENGTH(cols); i++) {
//     elem = INTEGER(cols)[i];
//     if (elem<1 || elem>LENGTH(x))
//       error("Item %d of 'cols' is %d which is outside 1-based range [1,ncol(x)=%d]", i+1, elem, LENGTH(x));
//     if (!n) n = length(VECTOR_ELT(x, elem-1));
//   }
//   SEXP ans = PROTECT(allocVector(LGLSXP, 1));
//   LOGICAL(ans)[0]=0;
//   for (i=0; i<LENGTH(cols); i++) {
//     SEXP v = VECTOR_ELT(x, INTEGER(cols)[i]-1);
//     if (!length(v) || isNewList(v) || isList(v)) continue; // like stats:::na.omit.data.frame, skip list/pairlist columns
//     if (n != length(v))
//       error("Column %d of input list x is length %d, inconsistent with first column of that item which is length %d.", i+1,length(v),n);
//     j=0;
//     switch (TYPEOF(v)) {
//     case LGLSXP: {
//       const int *iv = LOGICAL(v);
//       while(j < n && iv[j] != NA_LOGICAL) j++;
//       if (j < n) LOGICAL(ans)[0] = 1;
//     }
//       break;
//     case INTSXP: {
//       const int *iv = INTEGER(v);
//       while(j < n && iv[j] != NA_INTEGER) j++;
//       if (j < n) LOGICAL(ans)[0] = 1;
//     }
//       break;
//     case STRSXP: {
//       while (j < n && STRING_ELT(v, j) != NA_STRING) j++;
//       if (j < n) LOGICAL(ans)[0] = 1;
//     }
//       break;
//     case REALSXP: {
//       const double *dv = REAL(v);
//       if (INHERITS(v, char_integer64)) {
//         for (j=0; j != n; j++) {
//           if (DtoLL(dv[j]) == NA_INT64_LL) {
//             LOGICAL(ans)[0] = 1;
//             break;
//           }
//         }
//       } else {
//         while(j < n && !ISNAN(dv[j])) j++;
//         if (j < n) LOGICAL(ans)[0] = 1;
//       }
//     }
//       break;
//     case RAWSXP: {
//       // no such thing as a raw NA
//       // vector already initialised to all 0's
//     }
//       break;
//     case CPLXSXP: {
//       // taken from https://github.com/wch/r-source/blob/d75f39d532819ccc8251f93b8ab10d5b83aac89a/src/main/coerce.c
//       while (j < n && !ISNAN(COMPLEX(v)[j].r) && !ISNAN(COMPLEX(v)[j].i)) j++;
//       if (j < n) LOGICAL(ans)[0] = 1;
//     }
//       break;
//     default:
//       error("Unsupported column type '%s'", type2char(TYPEOF(v)));
//     }
//     if (LOGICAL(ans)[0]) break;
//   }
//   UNPROTECT(1);
//   return(ans);
// }

// from data.table_assign.c:
SEXP setcolorder(SEXP x, SEXP o)
{
  SEXP names = getAttrib(x, R_NamesSymbol);
  const int *od = INTEGER(o), ncol=LENGTH(x);
  if (isNull(names)) error("list passed to setcolorder has no names");
  if (ncol != LENGTH(names))
    error("Internal error: dt passed to setcolorder has %d columns but %d names", ncol, LENGTH(names));  // # nocov

  // Double-check here at C level that o[] is a strict permutation of 1:ncol. Reordering columns by reference makes no
  // difference to generations/refcnt so we can write behind barrier in this very special case of strict permutation.
  bool *seen = Calloc(ncol, bool);
  for (int i=0; i != ncol; ++i) {
    if (od[i]==NA_INTEGER || od[i]<1 || od[i]>ncol)
      error("Internal error: o passed to Csetcolorder contains an NA or out-of-bounds");  // # nocov
    if (seen[od[i]-1])
      error("Internal error: o passed to Csetcolorder contains a duplicate");             // # nocov
    seen[od[i]-1] = true;
  }
  Free(seen);

  SEXP *tmp = Calloc(ncol, SEXP);
  SEXP *xd = SEXPPTR(x), *namesd = STRING_PTR(names);
  for (int i=0; i != ncol; ++i) tmp[i] = xd[od[i]-1];
  memcpy(xd, tmp, ncol*sizeof(SEXP)); // sizeof is type size_t so no overflow here
  for (int i=0; i != ncol; ++i) tmp[i] = namesd[od[i]-1];
  memcpy(namesd, tmp, ncol*sizeof(SEXP));
  // No need to change key (if any); sorted attribute is column names not positions
  Free(tmp);
  return(R_NilValue);
}

SEXP keepattr(SEXP to, SEXP from)
{
  // Same as R_copyDFattr in src/main/attrib.c, but that seems not exposed in R's api
  // Only difference is that we reverse from and to in the prototype, for easier calling above
  SET_ATTRIB(to, ATTRIB(from));
  IS_S4_OBJECT(from) ?  SET_S4_OBJECT(to) : UNSET_S4_OBJECT(to);
  SET_OBJECT(to, OBJECT(from));
  return to;
}

SEXP growVector(SEXP x, R_len_t newlen)
{
  // Similar to EnlargeVector in src/main/subassign.c, with the following changes :
  // * replaced switch and loops with one memcpy for INTEGER and REAL, but need to age CHAR and VEC.
  // * no need to cater for names
  // * much shorter and faster
  SEXP newx;
  R_len_t i, len = length(x);
  if (isNull(x)) error("growVector passed NULL");
  PROTECT(newx = allocVector(TYPEOF(x), newlen));   // TO DO: R_realloc(?) here?
  if (newlen < len) len=newlen;   // i.e. shrink
  switch (TYPEOF(x)) {
  case STRSXP :
    for (i=0; i<len; i++)
      SET_STRING_ELT(newx, i, STRING_ELT(x, i));
    // TO DO. Using SET_ to ensure objects are aged, rather than memcpy. Perhaps theres a bulk/fast way to age CHECK_OLD_TO_NEW
    break;
  case VECSXP :
    for (i=0; i<len; i++)
      SET_VECTOR_ELT(newx, i, VECTOR_ELT(x, i));
    // TO DO: Again, is there bulk op to avoid this loop, which still respects older generations
    break;
  default :
    memcpy((char *)DATAPTR(newx), (char *)DATAPTR(x), len*SIZEOF(x));   // SIZEOF() returns size_t (just as sizeof()) so * shouldn't overflow // TODO remove DATAPTR
  }
  // if (verbose) Rprintf("Growing vector from %d to %d items of type '%s'\n", len, newlen, type2char(TYPEOF(x)));
  // Would print for every column if here. Now just up in dogroups (one msg for each table grow).
  keepattr(newx,x);
  UNPROTECT(1);
  return newx;
}


