/*
 This code is adapted from the data.table package: http://r-datatable.com
 and licensed under a Mozilla Public License 2.0 (MPL-2.0) license.
*/


#include "data.table.h"
#include <Rdefines.h>
#include <ctype.h>   // for isdigit
#include <Rmath.h> // from assign.c. needed ? - yes for functions like R_FINITE
// #include <stdint.h>

// Above code pasted from data.table_assign.c -> needed for rbindlist...

static bool anyNamed(SEXP x) {
  if (MAYBE_REFERENCED(x)) return true;
  if (isNewList(x)) { // fixed gcc10 issue through better indentation
    for (int i=0; i<LENGTH(x); i++) {
      if (anyNamed(VECTOR_ELT(x,i))) return true;
    }
  }
  return false;
}

#define MSGSIZE 1000
static char memrecycle_message[MSGSIZE+1]; // returned to rbindlist so it can prefix with which one of the list of data.table-like objects

const char *memrecycle(SEXP target, SEXP where, int start, int len, SEXP source, int colnum, const char *colname)
  // like memcpy but recycles single-item source
  // 'where' a 1-based INTEGER vector subset of target to assign to, or NULL or integer()
  // assigns to target[start:start+len-1] or target[where[start:start+len-1]] where start is 0-based
{
  if (len<1) return NULL;
  int slen = length(source);
  if (slen==0) return NULL;
  if (slen>1 && slen!=len && (!isNewList(target) || isNewList(source)))
    error("Internal error: recycle length error not caught earlier. slen=%d len=%d", slen, len); // # nocov
  // Internal error because the column has already been added to the DT, so length mismatch should have been caught before adding the column.
  // for 5647 this used to limit slen to len, but no longer
  if (colname==NULL)
    error("Internal error: memrecycle has received NULL colname"); // # nocov
  *memrecycle_message = '\0';
  int protecti=0;
  if (isNewList(source)) {
    // A list() column; i.e. target is a column of pointers to SEXPs rather than the much more common case
    // where memrecycle copies the DATAPTR data to the atomic target from the atomic source.
    // If any item within the list is NAMED then take a fresh copy. So far this has occurred from dogroups.c when
    // j returns .BY or similar specials as-is within a list(). Those specials are static inside
    // dogroups so if we don't copy now the last value written to them by dogroups becomes repeated in the result;
    // i.e. the wrong result.
    // If source is itself recycled later (many list() column items pointing to the same object) we are ok with that
    // since we now have a fresh copy and := will not assign with a list() column's cell value; := only changes the
    // SEXP pointed to.
    // If source is already not named (because j already created a fresh unnamed vector within a list()) we don't want to
    // duplicate unnecessarily, hence checking for named rather than duplicating always.
    // See #481, #1270 and tests 1341.* fail without this copy.
    if (anyNamed(source)) {
      source = PROTECT(copyAsPlain(source)); protecti++;
    }
  }
  const bool sourceIsFactor=isFactor(source), targetIsFactor=isFactor(target);
  const bool sourceIsI64=isReal(source) && Rinherits(source, char_integer64);
  const bool targetIsI64=isReal(target) && Rinherits(target, char_integer64);
  if (sourceIsFactor || targetIsFactor) {
    if (!targetIsFactor) {
      if (!isString(target) && !isNewList(target))
        error("Cannot assign 'factor' to '%s'. Factors can only be assigned to factor, character or list columns.", type2char(TYPEOF(target)));
      // else assigning factor to character is left to later below, avoiding wasteful asCharacterFactor
    } else if (!sourceIsFactor && !isString(source)) {
      // target is factor
      if (allNA(source, false)) {  // return false for list and other types that allNA does not support
        source = ScalarLogical(NA_LOGICAL); // a global constant in R and won't allocate; fall through to regular zero-copy coerce
      } else if (isInteger(source) || isReal(source)) {
        // allow assigning level numbers to factor columns; test 425, 426, 429 and 1945
        const int nlevel = length(getAttrib(target, R_LevelsSymbol));
        if (isInteger(source)) {
          const int *sd = INTEGER(source);
          for (int i=0; i<slen; ++i) {
            const int val = sd[i];
            if ((val<1 && val!=NA_INTEGER) || val>nlevel) {
              error("Assigning factor numbers to column %d named '%s'. But %d is outside the level range [1,%d]", colnum, colname, val, nlevel);
            }
          }
        } else {
          const double *sd = REAL(source);
          for (int i=0; i<slen; ++i) {
            const double val = sd[i];
            if (!ISNAN(val) && (!R_FINITE(val) || val!=(int)val || (int)val<1 || (int)val>nlevel)) {
              error("Assigning factor numbers to column %d named '%s'. But %f is outside the level range [1,%d], or is not a whole number.", colnum, colname, val, nlevel);
            }
          }
        }
        // Now just let the valid level numbers fall through to regular assign by BODY below
      } else {
        error("Cannot assign '%s' to 'factor'. Factor columns can be assigned factor, character, NA in any type, or level numbers.", type2char(TYPEOF(source)));
      }
    } else {
      // either factor or character being assigned to factor column
      SEXP targetLevels = PROTECT(getAttrib(target, R_LevelsSymbol)); protecti++;
      SEXP sourceLevels = source;  // character source
      if (sourceIsFactor) { sourceLevels=PROTECT(getAttrib(source, R_LevelsSymbol)); protecti++; }
      if (!sourceIsFactor || !R_compute_identical(sourceLevels, targetLevels, 0)) {  // !sourceIsFactor for test 2115.6
        const int nTargetLevels=length(targetLevels), nSourceLevels=length(sourceLevels);
        const SEXP *targetLevelsD=STRING_PTR(targetLevels), *sourceLevelsD=STRING_PTR(sourceLevels);
        SEXP newSource = PROTECT(allocVector(INTSXP, length(source))); protecti++;
        savetl_init();
        for (int k=0; k<nTargetLevels; ++k) {
          const SEXP s = targetLevelsD[k];
          const int tl = TRUELENGTH(s);
          if (tl>0) {
            savetl(s);
          } else if (tl<0) {
            // # nocov start
            for (int j=0; j<k; ++j) SET_TRUELENGTH(s, 0);  // wipe our negative usage and restore 0
            savetl_end();                                  // then restore R's own usage (if any)
            error("Internal error: levels of target are either not unique or have truelength<0");
            // # nocov end
          }
          SET_TRUELENGTH(s, -k-1);
        }
        int nAdd = 0;
        for (int k=0; k<nSourceLevels; ++k) {
          const SEXP s = sourceLevelsD[k];
          const int tl = TRUELENGTH(s);
          if (tl>=0) {
            if (!sourceIsFactor && s==NA_STRING) continue; // don't create NA factor level when assigning character to factor; test 2117
            if (tl>0) savetl(s);
            SET_TRUELENGTH(s, -nTargetLevels-(++nAdd));
          } // else, when sourceIsString, it's normal for there to be duplicates here
        }
        const int nSource = length(source);
        int *newSourceD = INTEGER(newSource);
        if (sourceIsFactor) {
          const int *sourceD = INTEGER(source);
          for (int i=0; i<nSource; ++i) {  // convert source integers to refer to target levels
            const int val = sourceD[i];
            newSourceD[i] = val==NA_INTEGER ? NA_INTEGER : -TRUELENGTH(sourceLevelsD[val-1]); // retains NA factor levels here via TL(NA_STRING); e.g. ordered factor
          }
        } else {
          const SEXP *sourceD = STRING_PTR(source);
          for (int i=0; i<nSource; ++i) {  // convert source integers to refer to target levels
            const SEXP val = sourceD[i];
            newSourceD[i] = val==NA_STRING ? NA_INTEGER : -TRUELENGTH(val);
          }
        }
        source = newSource;
        for (int k=0; k<nTargetLevels; ++k) SET_TRUELENGTH(targetLevelsD[k], 0);  // don't need those anymore
        if (nAdd) {
          // cannot grow the levels yet as that would be R call which could fail to alloc and we have no hook to clear up
          SEXP *temp = (SEXP *)malloc(nAdd * sizeof(SEXP *));
          if (!temp) {
            // # nocov start
            for (int k=0; k<nSourceLevels; ++k) SET_TRUELENGTH(sourceLevelsD[k], 0);
            savetl_end();
            error("Unable to allocate working memory of %d bytes to combine factor levels", nAdd*sizeof(SEXP *));
            // # nocov end
          }
          for (int k=0, thisAdd=0; thisAdd<nAdd; ++k) {   // thisAdd<nAdd to stop early when the added ones are all reached
            SEXP s = sourceLevelsD[k];
            int tl = TRUELENGTH(s);
            if (tl) {  // tl negative here
              if (tl != -nTargetLevels-thisAdd-1) error("Internal error: extra level check sum failed"); // # nocov
              temp[thisAdd++] = s;
              SET_TRUELENGTH(s,0);
            }
          }
          savetl_end();
          setAttrib(target, R_LevelsSymbol, targetLevels=growVector(targetLevels, nTargetLevels + nAdd));
          for (int k=0; k<nAdd; ++k) {
            SET_STRING_ELT(targetLevels, nTargetLevels+k, temp[k]);
          }
          free(temp);
        } else {
          // all source levels were already in target levels, but not with the same integers; we're done
          savetl_end();
        }
        // now continue, but with the mapped integers in the (new) source
      }
    }
  } else if (isString(source) && !isString(target) && !isNewList(target)) {
    warning("Coercing 'character' RHS to '%s' to match the type of the target column (column %d named '%s').",
            type2char(TYPEOF(target)), colnum, colname);
    // this "Coercing ..." warning first to give context in case coerceVector warns 'NAs introduced by coercion'
    source = PROTECT(coerceVector(source, TYPEOF(target))); protecti++;
  } else if (isNewList(source) && !isNewList(target)) {
    if (targetIsI64) {
      error("Cannot coerce 'list' RHS to 'integer64' to match the type of the target column (column %d named '%s').", colnum, colname);
      // because R's coerceVector doesn't know about integer64
    }
    // as in base R; e.g. let as.double(list(1,2,3)) work but not as.double(list(1,c(2,4),3))
    // relied on by NNS, simstudy and table.express; tests 1294.*
    warning("Coercing 'list' RHS to '%s' to match the type of the target column (column %d named '%s').",
            type2char(TYPEOF(target)), colnum, colname);
    source = PROTECT(coerceVector(source, TYPEOF(target))); protecti++;
  }
/*  Removing these checks because of too many -Wformat errors which I don't know how to solve... !!!
  else if ((TYPEOF(target)!=TYPEOF(source) || targetIsI64!=sourceIsI64) && !isNewList(target)) {
    // The following checks are up front here, otherwise we'd need them twice in the two branches
    //   inside BODY that cater for 'where' or not. Maybe there's a way to merge the two macros in future.
    // The idea is to do these range checks without calling coerceVector() (which allocates)
    // before "%"FMT" (typ .....

#define CHECK_RANGE(STYPE, RFUN, COND, FMT, TO) {{                                                                                \
    const STYPE *sd = (const STYPE *)RFUN(source);                                                                                \
    for (int i=0; i<slen; ++i) {                                                                                                  \
      const STYPE val = sd[i];                                                                                                    \
      if (COND) {                                                                                                                 \
        const char *sType = sourceIsI64 ? "integer64" : type2char(TYPEOF(source));                                                \
        const char *tType = targetIsI64 ? "integer64" : type2char(TYPEOF(target));                                                \
        int n = snprintf(memrecycle_message, MSGSIZE,                                                                             \
                         "%" FMT "' (type '%s') at RHS position %d " TO " when assigning to type '%s'", val, sType, i+1, tType);  \
        if (colnum>0 && n>0 && n<MSGSIZE)                                                                                         \
          snprintf(memrecycle_message+n, MSGSIZE-n, " (column %d named '%s')", colnum, colname);                                  \
        // string returned so that rbindlist/dogroups can prefix it with which item of its list this refers to                    \
        break;                                                                                                                    \
      }                                                                                                                           \
    }                                                                                                                             \
    } break; }

  switch(TYPEOF(target)) {
  case LGLSXP:
    switch (TYPEOF(source)) {
    case RAWSXP:  CHECK_RANGE(Rbyte, RAW,    val!=0 && val!=1,                                        "d",    "taken as TRUE")
    case INTSXP:  CHECK_RANGE(int, INTEGER,  val!=0 && val!=1 && val!=NA_INTEGER,                     "d",    "taken as TRUE")
    case REALSXP: if (sourceIsI64)
            CHECK_RANGE(int64_t, REAL, val!=0 && val!=1 && val!=NA_INTEGER64,                   PRId64, "taken as TRUE")
      else  CHECK_RANGE(double, REAL,  !ISNAN(val) && val!=0.0 && val!=1.0,                     "f",   "taken as TRUE")
    } break;
  case RAWSXP:
    switch (TYPEOF(source)) {
    case INTSXP:  CHECK_RANGE(int, INTEGER,  val<0 || val>255,                                        "d",    "taken as 0")
    case REALSXP: if (sourceIsI64)
            CHECK_RANGE(int64_t, REAL, val<0 || val>255,                                        PRId64, "taken as 0")
      else  CHECK_RANGE(double, REAL,  !R_FINITE(val) || val<0.0 || val>256.0 || (int)val!=val, "f",    "either truncated (precision lost) or taken as 0")
    } break;
  case INTSXP:
    if (TYPEOF(source)==REALSXP) {
      if (sourceIsI64)
              CHECK_RANGE(int64_t, REAL, val!=NA_INTEGER64 && (val<=NA_INTEGER || val>INT_MAX),   PRId64,  "out-of-range (NA)")
        else  CHECK_RANGE(double, REAL,  !ISNAN(val) && (!R_FINITE(val) || (int)val!=val),        "f",     "truncated (precision lost)")
    } break;
  case REALSXP:
    if (targetIsI64 && isReal(source) && !sourceIsI64) {
      CHECK_RANGE(double, REAL,  !ISNAN(val) && (!R_FINITE(val) || (int)val!=val),        "f",     "truncated (precision lost)")
    }
  }
  }
*/

#undef BODY
#define BODY(STYPE, RFUN, CTYPE, CAST, ASSIGN) {{    \
  const STYPE *sd = (const STYPE *)RFUN(source);     \
  if (length(where)) {                               \
    if (slen==1) {                                   \
      const STYPE val = sd[0];                       \
      const CTYPE cval = CAST;                       \
      for (int wi=0; wi<len; ++wi) {                 \
        const int w = wd[wi];                        \
        if (w<1) continue; /*0 or NA*/               \
        const int i = w-1;                           \
        ASSIGN;                                      \
      }                                              \
    } else {                                         \
      for (int wi=0; wi<len; ++wi) {                 \
        const int w = wd[wi];                        \
        if (w<1) continue;                           \
        const STYPE val = sd[wi];                    \
        const CTYPE cval = CAST;                     \
        const int i = w-1;                           \
        ASSIGN;                                      \
      }                                              \
    }                                                \
  } else {                                           \
    if (slen==1) {                                   \
      const STYPE val = sd[0];                       \
      const CTYPE cval = CAST;                       \
      for (int i=0; i<len; ++i) {                    \
        ASSIGN;                                      \
      }                                              \
    } else {                                         \
      for (int i=0; i<len; i++) {                    \
        const STYPE val = sd[i];                     \
        const CTYPE cval = CAST;                     \
        ASSIGN;                                      \
      }                                              \
    }                                                \
  }                                                  \
  } break; }

#define COERCE_ERROR(targetType) error("type '%s' cannot be coerced to '%s'", type2char(TYPEOF(source)), targetType); // 'targetType' for integer64 vs double

const int off = length(where) ? 0 : start;  // off = target offset; e.g. called from rbindlist with where=R_NilValue and start!=0
const bool mc = length(where)==0 && slen==len;  // mc=memcpy; only used if types match too
const int *wd = length(where) ? INTEGER(where)+start : NULL;
switch (TYPEOF(target)) {
case RAWSXP: {
  Rbyte *td = RAW(target) + off;
  switch (TYPEOF(source)) {
  case RAWSXP:
    if (mc) {
      memcpy(td, RAW(source), slen*sizeof(Rbyte)); break;
    } else        BODY(Rbyte, RAW,    Rbyte, val,                                     td[i]=cval)
  case LGLSXP:    BODY(int, LOGICAL,  Rbyte, val==1,                                  td[i]=cval)
  case INTSXP:    BODY(int, INTEGER,  Rbyte, (val>255 || val<0) ? 0 : val,            td[i]=cval)
  case REALSXP:
    if (sourceIsI64)
      BODY(int64_t, REAL, Rbyte, (val>255 || val<0) ? 0 : val,            td[i]=cval)
      else          BODY(double, REAL,  Rbyte, (ISNAN(val)||val>255||val<0) ? 0 : val,  td[i]=cval)
  default:        COERCE_ERROR("raw");
  }
} break;
case LGLSXP: {
  int *td = LOGICAL(target) + off;
  switch (TYPEOF(source)) {
  case RAWSXP:    BODY(Rbyte, RAW,    int, val!=0,                                    td[i]=cval)
  case LGLSXP:
    if (mc) {
      memcpy(td, LOGICAL(source), slen*sizeof(Rboolean)); break;
    } else        BODY(int, LOGICAL,  int, val,                                       td[i]=cval)
  case INTSXP:    BODY(int, INTEGER,  int, val==NA_INTEGER ? NA_LOGICAL : val!=0,     td[i]=cval)
  case REALSXP:
    if (sourceIsI64)
      BODY(int64_t, REAL, int, val==NA_INTEGER64 ? NA_LOGICAL : val!=0,   td[i]=cval)
      else          BODY(double,  REAL, int, ISNAN(val) ? NA_LOGICAL : val!=0.0,        td[i]=cval)
  default:        COERCE_ERROR("logical");
  }
} break;
case INTSXP : {
  int *td = INTEGER(target) + off;
  switch (TYPEOF(source)) {
  case  RAWSXP:   BODY(Rbyte, RAW,    int, (int)val,                                  td[i]=cval)
  case  LGLSXP:   // same as INTSXP ...
  case  INTSXP:
    if (mc) {
      memcpy(td, INTEGER(source), slen*sizeof(int)); break;
    } else        BODY(int, INTEGER,  int, val,                                       td[i]=cval)
  case REALSXP:
    if (sourceIsI64)
      BODY(int64_t, REAL, int, (val==NA_INTEGER64||val>INT_MAX||val<=NA_INTEGER) ? NA_INTEGER : (int)val,  td[i]=cval)
      else          BODY(double, REAL,  int, ISNAN(val) ? NA_INTEGER : (int)val,        td[i]=cval)
  default:        COERCE_ERROR("integer"); // test 2005.4
  }
} break;
case REALSXP : {
  if (targetIsI64) {
  int64_t *td = (int64_t *)REAL(target) + off;
  switch (TYPEOF(source)) {
  case RAWSXP:  BODY(Rbyte, RAW,    int64_t, (int64_t)val,                          td[i]=cval)
  case LGLSXP:  // same as INTSXP
  case INTSXP:  BODY(int, INTEGER,  int64_t, val==NA_INTEGER ? NA_INTEGER64 : val,  td[i]=cval)
  case REALSXP:
    if (sourceIsI64) {
      if (mc) {
        memcpy(td, (int64_t *)REAL(source), slen*sizeof(int64_t)); break;
      } else    BODY(int64_t, REAL, int64_t, val,                                   td[i]=cval)
    } else      BODY(double, REAL,  int64_t, R_FINITE(val) ? val : NA_INTEGER64,    td[i]=cval)
  default:      COERCE_ERROR("integer64");
  }
} else {
  double *td = REAL(target) + off;
  switch (TYPEOF(source)) {
  case  RAWSXP: BODY(Rbyte, RAW,    double, (double)val,                            td[i]=cval)
  case  LGLSXP: // same as INTSXP
  case  INTSXP: BODY(int, INTEGER,  double, val==NA_INTEGER ? NA_REAL : val,        td[i]=cval)
  case REALSXP:
    if (!sourceIsI64) {
      if (mc) {
        memcpy(td, (double *)REAL(source), slen*sizeof(double)); break;
      } else    BODY(double, REAL,  double, val,                                    td[i]=cval)
    } else      BODY(int64_t, REAL, double, val==NA_INTEGER64 ? NA_REAL : val,      td[i]=cval)
  default:      COERCE_ERROR("double");
  }
}
} break;
case CPLXSXP: {
  Rcomplex *td = COMPLEX(target) + off;
  double im = 0.0;
  switch (TYPEOF(source)) {
  case  RAWSXP:   BODY(Rbyte, RAW,    double, (im=0.0,val),                                         td[i].r=cval;td[i].i=im)
  case  LGLSXP:   // same as INTSXP
  case  INTSXP:   BODY(int, INTEGER,  double, val==NA_INTEGER?(im=NA_REAL,NA_REAL):(im=0.0,val),    td[i].r=cval;td[i].i=im)
  case REALSXP:
    if (sourceIsI64)
      BODY(int64_t, REAL, double, val==NA_INTEGER64?(im=NA_REAL,NA_REAL):(im=0.0,val),  td[i].r=cval;td[i].i=im)
      else          BODY(double,  REAL, double, ISNAN(val)?(im=NA_REAL,NA_REAL):(im=0.0,val),         td[i].r=cval;td[i].i=im)
  case CPLXSXP:
    if (mc) {
      memcpy(td, COMPLEX(source), slen*sizeof(Rcomplex)); break;
    } else        BODY(Rcomplex, COMPLEX, Rcomplex, val,                                            td[i]=cval)
  default:        COERCE_ERROR("complex");
  }
} break;
case STRSXP :
  if (sourceIsFactor) {
    const SEXP *ld = STRING_PTR(PROTECT(getAttrib(source, R_LevelsSymbol))); protecti++;
    BODY(int, INTEGER, SEXP, val==NA_INTEGER ? NA_STRING : ld[val-1],  SET_STRING_ELT(target, off+i, cval))
  } else {
    if (!isString(source)) {
      if (allNA(source, true)) {  // saves common coercion of NA (logical) to NA_character_
        //              ^^ =errorForBadType; if type list, that was already an error earlier so we
        //                 want to be strict now otherwise list would get to coerceVector below
        if (length(where)) {
          for (int i=0; i<len; ++i) if (wd[i]>0) SET_STRING_ELT(target, wd[i]-1, NA_STRING);
        } else {
          for (int i=0; i<len; ++i) SET_STRING_ELT(target, start+i, NA_STRING);
        }
        break;
      }
      if (sourceIsI64)
        error("To assign integer64 to a character column, please use as.character() for clarity.");
      source = PROTECT(coerceVector(source, STRSXP)); protecti++;
    }
    BODY(SEXP, STRING_PTR, SEXP, val,  SET_STRING_ELT(target, off+i, cval))
  }
case VECSXP :
  if (TYPEOF(source)!=VECSXP)
    BODY(SEXP, &, SEXP, val,           SET_VECTOR_ELT(target, off+i, cval))
    else
      BODY(SEXP, SEXPPTR, SEXP, val,  SET_VECTOR_ELT(target, off+i, cval))
default :
      error("Unsupported column type in assign.c:memrecycle '%s'", type2char(TYPEOF(target)));  // # nocov
}
UNPROTECT(protecti);
return memrecycle_message[0] ? memrecycle_message : NULL;
}

void writeNA(SEXP v, const int from, const int n)
  // this is for use after allocVector() which does not initialize its result. It does write NA as you'd
  // think, other than for VECSXP which allocVector() already initializes with NULL.
{
  const int to = from-1+n;  // together with <=to below with writing NA to position 2147483647 in mind
  switch(TYPEOF(v)) {
  case RAWSXP:
    memset(RAW(v)+from, 0, n*SIZEOF(v));
    break;
  case LGLSXP : {
    Rboolean *vd = (Rboolean *)LOGICAL(v);
    for (int i=from; i<=to; ++i) vd[i] = NA_LOGICAL;
  } break;
  case INTSXP : {
    // same whether factor or not
    int *vd = INTEGER(v);
    for (int i=from; i<=to; ++i) vd[i] = NA_INTEGER;
  } break;
  case REALSXP : {
    if (INHERITS(v, char_integer64)) {
    int64_t *vd = (int64_t *)REAL(v);
    for (int i=from; i<=to; ++i) vd[i] = INT64_MIN;
  } else {
    double *vd = REAL(v);
    for (int i=from; i<=to; ++i) vd[i] = NA_REAL;
  }
  } break;
  case CPLXSXP: {
    Rcomplex *vd = COMPLEX(v);
    for (int i=from; i<=to; ++i) vd[i] = NA_CPLX;
  } break;
  case STRSXP :
    // character columns are initialized with blank string (""). So replace the all-"" with all-NA_character_
    // Since "" and NA_character_ are global constants in R, it should be ok to not use SET_STRING_ELT here. But use it anyway for safety (revisit if proved slow)
    // If there's ever a way added to R API to pass NA_STRING to allocVector() to tell it to initialize with NA not "", would be great
    for (int i=from; i<=to; ++i) SET_STRING_ELT(v, i, NA_STRING);
    break;
  case VECSXP : case EXPRSXP :
    // list & expression columns already have each item initialized to NULL
    break;
  default :
    error("Internal error: writeNA passed a vector of type '%s'", type2char(TYPEOF(v)));  // # nocov
  }
}

static SEXP *saveds=NULL;
static R_len_t *savedtl=NULL, nalloc=0, nsaved=0;

void savetl_init() {
  if (nsaved || nalloc || saveds || savedtl) {
    error("Internal error: savetl_init checks failed (%d %d %p %p). please report to data.table issue tracker.", nsaved, nalloc, saveds, savedtl); // # nocov
  }
  nsaved = 0;
  nalloc = 100;
  saveds = (SEXP *)malloc(nalloc * sizeof(SEXP));
  savedtl = (R_len_t *)malloc(nalloc * sizeof(R_len_t));
  if (saveds==NULL || savedtl==NULL) {
    savetl_end();                                                        // # nocov
    error("Failed to allocate initial %d items in savetl_init", nalloc); // # nocov
  }
}

void savetl(SEXP s)
{
  if (nsaved==nalloc) {
    if (nalloc==INT_MAX) {
      savetl_end();                                                                                                     // # nocov
      error("Internal error: reached maximum %d items for savetl. Please report to data.table issue tracker.", nalloc); // # nocov
    }
    nalloc = nalloc>(INT_MAX/2) ? INT_MAX : nalloc*2;
    char *tmp = (char *)realloc(saveds, nalloc*sizeof(SEXP));
    if (tmp==NULL) {
      // C spec states that if realloc() fails the original block is left untouched; it is not freed or moved. We rely on that here.
      savetl_end();                                                      // # nocov  free(saveds) happens inside savetl_end
      error("Failed to realloc saveds to %d items in savetl", nalloc);   // # nocov
    }
    saveds = (SEXP *)tmp;
    tmp = (char *)realloc(savedtl, nalloc*sizeof(R_len_t));
    if (tmp==NULL) {
      savetl_end();                                                      // # nocov
      error("Failed to realloc savedtl to %d items in savetl", nalloc);  // # nocov
    }
    savedtl = (R_len_t *)tmp;
  }
  saveds[nsaved] = s;
  savedtl[nsaved] = TRUELENGTH(s);
  nsaved++;
}

void savetl_end() {
  // Can get called if nothing has been saved yet (nsaved==0), or even if _init() hasn't been called yet (pointers NULL). Such
  // as to clear up before error. Also, it might be that nothing needed to be saved anyway.
  for (int i=0; i<nsaved; i++) SET_TRUELENGTH(saveds[i],savedtl[i]);
  free(saveds);  // possible free(NULL) which is safe no-op
  saveds = NULL;
  free(savedtl);
  savedtl = NULL;
  nsaved = nalloc = 0;
}

// now this was in rbindlist.c !!

SEXP rbindlist(SEXP l, SEXP usenamesArg, SEXP fillArg, SEXP idcolArg)
{
  if (!isLogical(fillArg) || LENGTH(fillArg) != 1 || LOGICAL(fillArg)[0] == NA_LOGICAL)
    error("fill= should be TRUE or FALSE");
  if (!isLogical(usenamesArg) || LENGTH(usenamesArg)!=1)
    error("use.names= should be TRUE, FALSE, or not used (\"check\" by default)");  // R levels converts "check" to NA
  if (!length(l)) return(l);
  if (TYPEOF(l) != VECSXP) error("Input to rbindlist must be a list. This list can contain data.tables, data.frames or plain lists.");
  Rboolean usenames = LOGICAL(usenamesArg)[0];
  const bool fill = LOGICAL(fillArg)[0];
  if (fill && usenames!=TRUE) {
    if (usenames==FALSE) warning("use.names= cannot be FALSE when fill is TRUE. Setting use.names=TRUE."); // else no warning if usenames==NA (default)
    usenames=TRUE;
  }
  const bool idcol = !isNull(idcolArg);
  if (idcol && (!isString(idcolArg) || LENGTH(idcolArg)!=1)) error("Internal error: rbindlist.c idcol is not a single string");  // # nocov
  int ncol=0, first=0;
  int64_t nrow=0, upperBoundUniqueNames=1;
  bool anyNames=false;
  int numZero=0, firstZeroCol=0, firstZeroItem=0;
  int *eachMax = (int *)R_alloc(LENGTH(l), sizeof(int));
  // pre-check for any errors here to save having to get cleanup right below when usenames
  for (int i=0; i<LENGTH(l); i++) {  // length(l)>0 checked above
    eachMax[i] = 0;
    SEXP li = VECTOR_ELT(l, i);
    if (isNull(li)) continue;
    if (TYPEOF(li) != VECSXP) error("Item %d of input is not a data.frame, data.table or list", i+1);
    const int thisncol = length(li);
    if (!thisncol) continue;
    // delete as now more flexible ... if (fill && isNull(getAttrib(li, R_NamesSymbol))) error("When fill=TRUE every item of the input must have column names. Item %d does not.", i+1);
    if (fill) {
      if (thisncol>ncol) ncol=thisncol;  // this section initializes ncol with max ncol. ncol may be increased when usenames is accounted for further down
    } else {
      if (ncol==0) { ncol=thisncol; first=i; }
      else if (thisncol!=ncol) error("Item %d has %d columns, inconsistent with item %d which has %d columns. To fill missing columns use fill=TRUE.", i+1, thisncol, first+1, ncol);
    }
    int nNames = length(getAttrib(li, R_NamesSymbol));
    if (nNames>0 && nNames!=thisncol) error("Item %d has %d columns but %d column names. Invalid object.", i+1, thisncol, nNames);
    if (nNames>0) anyNames=true;
    upperBoundUniqueNames += nNames;
    int maxLen=0, whichMax=0;
    for (int j=0; j<thisncol; ++j) { int tt=length(VECTOR_ELT(li,j)); if (tt>maxLen) { maxLen=tt; whichMax=j; } }
    for (int j=0; j<thisncol; ++j) {
      int tt = length(VECTOR_ELT(li, j));
      if (tt>1 && tt!=maxLen) error("Column %d of item %d is length %d inconsistent with column %d which is length %d. Only length-1 columns are recycled.", j+1, i+1, tt, whichMax+1, maxLen);
      if (tt==0 && maxLen>0 && numZero++==0) { firstZeroCol = j; firstZeroItem=i; }
    }
    eachMax[i] = maxLen;
    nrow += maxLen;
  }
  if (numZero) {  // #1871
    SEXP names = getAttrib(VECTOR_ELT(l, firstZeroItem), R_NamesSymbol);
    const char *ch = names==R_NilValue ? "" : CHAR(STRING_ELT(names, firstZeroCol));
    warning("Column %d ['%s'] of item %d is length 0. This (and %d other%s like it) has been filled with NA (NULL for list columns) to make each item uniform.",
            firstZeroCol+1, ch, firstZeroItem+1, numZero-1, numZero==2?"":"s");
  }
  if (nrow==0 && ncol==0) return(R_NilValue);
  if (nrow>INT32_MAX) error("Total rows in the list is %lld which is larger than the maximum number of rows, currently %d", nrow, INT32_MAX);
  if (usenames==TRUE && !anyNames) error("use.names=TRUE but no item of input list has any names");

  int *colMap=NULL; // maps each column in final result to the column of each list item
  if (usenames==TRUE || usenames==NA_LOGICAL) {
    // here we proceed as if fill=true for brevity (accounting for dups is tricky) and then catch any missings after this branch
    // when use.names==NA we also proceed here as if use.names was TRUE to save new code and then check afterwards the map is 1:ncol for every item
    // first find number of unique column names present; i.e. length(unique(unlist(lapply(l,names))))
    SEXP *uniq = (SEXP *)malloc(upperBoundUniqueNames * sizeof(SEXP));  // upperBoundUniqueNames was initialized with 1 to ensure this is defined (otherwise 0 when no item has names)
    if (!uniq) error("Failed to allocate upper bound of %lld unique column names [sum(lapply(l,ncol))]", upperBoundUniqueNames);
    savetl_init();
    int nuniq=0;
    for (int i=0; i<LENGTH(l); i++) {
      SEXP li = VECTOR_ELT(l, i);
      int thisncol=LENGTH(li);
      if (isNull(li) || !LENGTH(li)) continue;
      const SEXP cn = getAttrib(li, R_NamesSymbol);
      if (!length(cn)) continue;
      const SEXP *cnp = STRING_PTR(cn);
      for (int j=0; j<thisncol; j++) {
        SEXP s = cnp[j];
        if (TRUELENGTH(s)<0) continue;  // seen this name before
        if (TRUELENGTH(s)>0) savetl(s);
        uniq[nuniq++] = s;
        SET_TRUELENGTH(s,-nuniq);
      }
    }
    if (nuniq>0) {
      SEXP *tt = realloc(uniq, nuniq*sizeof(SEXP));  // shrink to only what we need to release the spare
      if (!tt) free(uniq);  // shrink never fails; just keep codacy happy
      uniq = tt;
    }
    // now count the dups (if any) and how they're distributed across the items
    int *counts = (int *)calloc(nuniq, sizeof(int)); // counts of names for each colnames
    int *maxdup = (int *)calloc(nuniq, sizeof(int)); // the most number of dups for any name within one colname vector
    if (!counts || !maxdup) {
      // # nocov start
      for (int i=0; i<nuniq; ++i) SET_TRUELENGTH(uniq[i], 0);
      free(uniq); free(counts); free(maxdup);
      savetl_end();
      error("Failed to allocate nuniq=%d items working memory in rbindlist.c", nuniq);
      // # nocov end
    }
    for (int i=0; i<LENGTH(l); i++) {
      SEXP li = VECTOR_ELT(l, i);
      int thisncol=length(li);
      if (thisncol==0) continue;
      const SEXP cn = getAttrib(li, R_NamesSymbol);
      if (!length(cn)) continue;
      const SEXP *cnp = STRING_PTR(cn);
      memset(counts, 0, nuniq*sizeof(int));
      for (int j=0; j<thisncol; j++) {
        SEXP s = cnp[j];
        counts[ -TRUELENGTH(s)-1 ]++;
      }
      for (int u=0; u<nuniq; u++) {
        if (counts[u] > maxdup[u]) maxdup[u] = counts[u];
      }
    }
    int ttncol = 0;
    for (int u=0; u<nuniq; ++u) ttncol+=maxdup[u];
    if (ttncol>ncol) ncol=ttncol;
    free(maxdup); maxdup=NULL;  // not needed again
    // ncol is now the final number of columns accounting for unique and dups across all colnames
    // allocate a matrix:  nrows==length(list)  each entry contains which column to fetch for that final column

    int *colMapRaw = (int *)malloc(LENGTH(l)*ncol * sizeof(int));  // the result of this scope used later
    int *uniqMap = (int *)malloc(ncol * sizeof(int)); // maps the ith unique string to the first time it occurs in the final result
    int *dupLink = (int *)malloc(ncol * sizeof(int)); // if a colname has occurred before (a dup) links from the 1st to the 2nd time in the final result, 2nd to 3rd, etc
    if (!colMapRaw || !uniqMap || !dupLink) {
      // # nocov start
      for (int i=0; i<nuniq; ++i) SET_TRUELENGTH(uniq[i], 0);
      free(uniq); free(counts); free(colMapRaw); free(uniqMap); free(dupLink);
      savetl_end();
      error("Failed to allocate ncol=%d items working memory in rbindlist.c", ncol);
      // # nocov end
    }
    for (int i=0; i<LENGTH(l)*ncol; ++i) colMapRaw[i]=-1;   // 0-based so use -1
    for (int i=0; i<ncol; ++i) {uniqMap[i] = dupLink[i] = -1;}
    int nextCol=0, lastDup=ncol-1;

    for (int i=0; i<LENGTH(l); ++i) {
      SEXP li = VECTOR_ELT(l, i);
      int thisncol=length(li);
      if (thisncol==0) continue;
      const SEXP cn = getAttrib(li, R_NamesSymbol);
      if (!length(cn)) {
        for (int j=0; j<thisncol; j++) colMapRaw[i*ncol + j] = j;
      } else {
        const SEXP *cnp = STRING_PTR(cn);
        memset(counts, 0, nuniq*sizeof(int));
        for (int j=0; j<thisncol; j++) {
          SEXP s = cnp[j];
          int w = -TRUELENGTH(s)-1;
          int wi = counts[w]++; // how many dups have we seen before of this name within this item
          if (uniqMap[w]==-1) {
            // first time seen this name across all items
            uniqMap[w] = nextCol++;
          } else {
            while (wi && dupLink[w]>0) { w=dupLink[w]; --wi; }  // hop through the dups
            if (wi && dupLink[w]==-1) {
              // first time we've seen this number of dups of this name
              w = dupLink[w] = lastDup--;
              uniqMap[w] = nextCol++;
            }
          }
          colMapRaw[i*ncol + uniqMap[w]] = j;
        }
      }
    }
    for (int i=0; i<nuniq; ++i) SET_TRUELENGTH(uniq[i], 0);  // zero out our usage of tl
    free(uniq); free(counts); free(uniqMap); free(dupLink);  // all local scope so no need to set to NULL
    savetl_end();  // restore R's usage

    // colMapRaw is still allocated. It was allocated with malloc because we needed to catch if the alloc failed.
    // move it to R's heap so it gets automatically free'd on exit, and on any error between now and the end of rbindlist.
    colMap = (int *)R_alloc(LENGTH(l)*ncol, sizeof(int));
    // This R_alloc could fail with out-of-memory but given it is very small it's very unlikely. If it does fail, colMapRaw will leak.
    //   But colMapRaw leaking now in this very rare situation is better than colMapRaw leaking in the more likely but still rare conditions later.
    //   And it's better than having to trap all exit point from here to the end of rbindlist, which may not be possible; e.g. writeNA() could error inside it with unsupported type.
    //   This very unlikely leak could be fixed by using an on.exit() at R level rbindlist(); R-exts$6.1.2 refers to pwilcox for example. However, that would not
    //   solve the (mere) leak if we ever call rbindlist internally from other C functions.
    memcpy(colMap, colMapRaw, LENGTH(l)*ncol*sizeof(int));
    free(colMapRaw);  // local scope in this branch to ensure can't be used below

    // to view map when debugging ...
    // for (int i=0; i<LENGTH(l); ++i) { for (int j=0; j<ncol; ++j) Rprintf("%2d ",colMap[i*ncol + j]);  Rprintf("\n"); }
  }

  if (fill && usenames==NA_LOGICAL) error("Internal error: usenames==NA but fill=TRUE. usenames should have been set to TRUE earlier with warning.");
  if (!fill && (usenames==TRUE || usenames==NA_LOGICAL)) {
    // Ensure no missings in both cases, and (when usenames==NA) all columns in same order too
    // We proceeded earlier as if fill was true, so varying ncol items will have missings here
    char buff[1001] = "";
    const char *extra = usenames==TRUE?"":" use.names='check' (default from v1.12.2) emits this message and proceeds as if use.names=FALSE for "\
                                          " backwards compatibility. See news item 5 in v1.12.2 for options to control this message.";
    for (int i=0; i<LENGTH(l); ++i) {
      SEXP li = VECTOR_ELT(l, i);
      if (!length(li) || !length(getAttrib(li, R_NamesSymbol))) continue;
      for (int j=0; j<ncol; ++j) {
        const int w = colMap[i*ncol + j];
        if (w==-1) {
          int missi = i;
          while (colMap[i*ncol + j]==-1 && i<LENGTH(l)) i++;
          if (i==LENGTH(l)) error("Internal error: could not find the first column name not present in earlier item");
          SEXP s = getAttrib(VECTOR_ELT(l, i), R_NamesSymbol);
          int w2 = colMap[i*ncol + j];
          const char *str = isString(s) ? CHAR(STRING_ELT(s,w2)) : "";
          snprintf(buff, 1000, "Column %d ['%s'] of item %d is missing in item %d. Use fill=TRUE to fill with NA (NULL for list columns), or use.names=FALSE to ignore column names.%s",
                        w2+1, str, i+1, missi+1, extra );
          if (usenames==TRUE) error(buff);
          i = LENGTH(l); // break from outer i loop
          break;         // break from inner j loop
        }
        if (w!=j && usenames==NA_LOGICAL) {
          SEXP s = getAttrib(VECTOR_ELT(l, i), R_NamesSymbol);
          if (!isString(s) || i==0) error("Internal error: usenames==NA but an out-of-order name has been found in an item with no names or the first item. [%d]", i);
          snprintf(buff, 1000, "Column %d ['%s'] of item %d appears in position %d in item %d. Set use.names=TRUE to match by column name, or use.names=FALSE to ignore column names.%s",
                               w+1, CHAR(STRING_ELT(s,w)), i+1, j+1, i, extra);
          i = LENGTH(l);
          break;
        }
      }
//      if (buff[0]) {
//        SEXP opt = GetOption(install("datatable.rbindlist.check"), R_NilValue);
//        if (!isNull(opt) && !(isString(opt) && length(opt)==1)) {
//          warning("options()$datatable.rbindlist.check is set but is not a single string. See news item 5 in v1.12.2.");
//          opt = R_NilValue;
//        }
//        const char *o = isNull(opt) ? "message" : CHAR(STRING_ELT(opt,0));
//        if      (strcmp(o,"message")==0) { eval(PROTECT(lang2(install("message"),PROTECT(ScalarString(mkChar(buff))))), R_GlobalEnv); UNPROTECT(2); }
//        else if (strcmp(o,"warning")==0) warning(buff);
//        else if (strcmp(o,"error")==0)   error(buff);
//        else if (strcmp(o,"none")!=0)    warning("options()$datatable.rbindlist.check=='%s' which is not 'message'|'warning'|'error'|'none'. See news item 5 in v1.12.2.", o);
//      }
    }
  }
  if (usenames==NA_LOGICAL) {
    usenames=FALSE;  // for backwards compatibility, see warning above which says this will change to TRUE in future
    ncol = length(VECTOR_ELT(l, first));  // ncol was increased as if fill=true, so reduce it back given fill=false (fill==false checked above)
  }

  int nprotect = 0;
  SEXP ans = PROTECT(allocVector(VECSXP, idcol + ncol)); nprotect++;
  SEXP ansNames;
  setAttrib(ans, R_NamesSymbol, ansNames=allocVector(STRSXP, idcol + ncol));
  if (idcol) {
    SET_STRING_ELT(ansNames, 0, STRING_ELT(idcolArg, 0));
    SEXP idval, listNames=getAttrib(l, R_NamesSymbol);
    if (length(listNames)) {
      SET_VECTOR_ELT(ans, 0, idval=allocVector(STRSXP, nrow));
      for (int i=0,ansloc=0; i<LENGTH(l); ++i) {
        SEXP li = VECTOR_ELT(l, i);
        if (!length(li)) continue;
        const int thisnrow = eachMax[i];
        SEXP thisname = STRING_ELT(listNames, i);
        for (int k=0; k<thisnrow; ++k) SET_STRING_ELT(idval, ansloc++, thisname);
      }
    } else {
      SET_VECTOR_ELT(ans, 0, idval=allocVector(INTSXP, nrow));
      int *idvald = INTEGER(idval);
      for (int i=0,ansloc=0; i<LENGTH(l); ++i) {
        SEXP li = VECTOR_ELT(l, i);
        if (!length(li)) continue;
        const int thisnrow = eachMax[i];
        for (int k=0; k<thisnrow; ++k) idvald[ansloc++] = i+1;
      }
    }
  }

  SEXP coercedForFactor = NULL;
  for(int j=0; j<ncol; ++j) {
    int maxType=LGLSXP;  // initialize with LGLSXP for test 2002.3 which has col x NULL in both lists to be filled with NA for #1871
    bool factor=false, orderedFactor=false;     // ordered factor is class c("ordered","factor"). isFactor() is true when isOrdered() is true.
    int longestLen=0, longestW=-1, longestI=-1; // just for ordered factor
    SEXP longestLevels=R_NilValue;              // just for ordered factor
    bool int64=false;
    const char *foundName=NULL;
    bool anyNotStringOrFactor=false;
    SEXP firstCol=R_NilValue;
    int firsti=-1, firstw=-1;
    for (int i=0; i<LENGTH(l); ++i) {
      SEXP li = VECTOR_ELT(l, i);
      if (!length(li)) continue;
      int w = usenames ? colMap[i*ncol + j] : j;  // colMap tells us which item to fetch for each of the final result columns, so we can stack column-by-column
      if (w==-1) continue;  // column j of final result has no input from this item (fill must be true)
      if (!foundName) {
        SEXP cn=PROTECT(getAttrib(li, R_NamesSymbol));
        if (length(cn)) { SEXP tt; SET_STRING_ELT(ansNames, idcol+j, tt=STRING_ELT(cn, w)); foundName=CHAR(tt); }
        UNPROTECT(1);
      }
      SEXP thisCol = VECTOR_ELT(li, w);
      int thisType = TYPEOF(thisCol);
      if (TYPEORDER(thisType)>TYPEORDER(maxType)) maxType=thisType;
      // return ScalarInteger(maxType);
      if (isFactor(thisCol)) {
        if (isNull(getAttrib(thisCol,R_LevelsSymbol))) error("Column %d of item %d has type 'factor' but has no levels; i.e. malformed.", w+1, i+1);
        factor = true;
        if (isOrdered(thisCol)) {
          orderedFactor = true;
          int thisLen = length(getAttrib(thisCol, R_LevelsSymbol));
          if (thisLen>longestLen) { longestLen=thisLen; longestLevels=getAttrib(thisCol, R_LevelsSymbol); /*for warnings later ...*/longestW=w; longestI=i; }
        }
      } else if (!isString(thisCol)) anyNotStringOrFactor=true;  // even for length 0 columns for consistency; test 2113.3
      if (INHERITS(thisCol, char_integer64)) {  // PRINTNAME(install("integer64"))
        if (firsti>=0 && !length(getAttrib(firstCol, R_ClassSymbol))) { firsti=i; firstw=w; firstCol=thisCol; } // so the integer64 attribute gets copied to target below
        int64=true;
      }
      if (firsti==-1) { firsti=i; firstw=w; firstCol=thisCol; }
      else {
        if (!factor && !int64) {
          if (!R_compute_identical(PROTECT(getAttrib(thisCol, R_ClassSymbol)),
                                   PROTECT(getAttrib(firstCol, R_ClassSymbol)),
                                   0)) {
            error("Class attribute on column %d of item %d does not match with column %d of item %d.", w+1, i+1, firstw+1, firsti+1);
          }
          UNPROTECT(2);
        }
      }
    }

    if (!foundName) { static char buff[12]; sprintf(buff,"V%d",j+1), SET_STRING_ELT(ansNames, idcol+j, mkChar(buff)); foundName=buff; }
    if (factor) maxType=INTSXP;  // if any items are factors then a factor is created (could be an option)
    if (int64 && maxType!=REALSXP)
      error("Internal error: column %d of result is determined to be integer64 but maxType=='%s' != REALSXP", j+1, type2char(maxType)); // # nocov
    SEXP target;
    SET_VECTOR_ELT(ans, idcol+j, target=allocVector(maxType, nrow));  // does not initialize logical & numerics, but does initialize character and list
    if (!factor) copyMostAttrib(firstCol, target); // all but names,dim and dimnames; mainly for class. And if so, we want a copy here, not keepattr's SET_ATTRIB.

    if (factor && anyNotStringOrFactor) {
      // in future warn, or use list column instead ... warning("Column %d contains a factor but not all items for the column are character or factor", idcol+j+1);
      // some coercing from (likely) integer/numeric to character will be needed. But this coerce can feasibly fail with out-of-memory, so we have to do it up-front
      // before the savetl_init() because we have no hook to clean up tl if coerceVector fails.
      if (coercedForFactor==NULL) { coercedForFactor=PROTECT(allocVector(VECSXP, LENGTH(l))); nprotect++; }
      for (int i=0; i<LENGTH(l); ++i) {
        int w = usenames ? colMap[i*ncol + j] : j;
        if (w==-1) continue;
        SEXP thisCol = VECTOR_ELT(VECTOR_ELT(l, i), w);
        if (!isFactor(thisCol) && !isString(thisCol)) {
          SET_VECTOR_ELT(coercedForFactor, i, coerceVector(thisCol, STRSXP));
        }
      }
    }
    int ansloc=0;
    if (factor) {
      char warnStr[1000] = "";
      savetl_init();  // no error from now (or warning given options(warn=2)) until savetl_end
      int nLevel=0, allocLevel=0;
      SEXP *levelsRaw = NULL;  // growing list of SEXP pointers. Raw since managed with raw realloc.
      if (orderedFactor) {
        // If all sets of ordered levels are compatible (no ambiguities or conflicts) then an ordered factor is created, otherwise regular factor.
        // Currently the longest set of ordered levels is taken and all other ordered levels must be a compatible subset of that.
        // e.g. c( a<c<b, z<a<c<b, a<b ) => z<a<c<b  [ the longest is the middle one, and the other two are ordered subsets of it ]
        //      c( a<c<b, z<c<a<b, a<b ) => regular factor because it contains an ambiguity: is a<c or c<a?
        //      c( a<c<b, c<b, 'c,b'   ) => a<c<b  because the regular factor/character items c and b exist in the ordered levels
        //      c( a<c<b, c<b, 'c,d'   ) => a<c<b<d  'd' from non-ordered item added on the end of longest ordered levels
        //      c( a<c<b, c<b<d<e )  => regular factor because this case isn't yet implemented. a<c<b<d<e would be possible in future (extending longest at the beginning or end)
        const SEXP *sd = STRING_PTR(longestLevels);
        nLevel = allocLevel = longestLen;
        levelsRaw = (SEXP *)malloc(nLevel * sizeof(SEXP));
        if (!levelsRaw) { savetl_end(); error("Failed to allocate working memory for %d ordered factor levels of result column %d", nLevel, idcol+j+1); }
        for (int k=0; k<longestLen; ++k) {
          SEXP s = sd[k];
          if (TRUELENGTH(s)>0) savetl(s);
          levelsRaw[k] = s;
          SET_TRUELENGTH(s,-k-1);
        }
        for (int i=0; i<LENGTH(l); ++i) {
          int w = usenames ? colMap[i*ncol + j] : j;
          if (w==-1) continue;
          SEXP thisCol = VECTOR_ELT(VECTOR_ELT(l, i), w);
          if (isOrdered(thisCol)) {
            SEXP levels = getAttrib(thisCol, R_LevelsSymbol);
            const SEXP *levelsD = STRING_PTR(levels);
            const int n = length(levels);
            for (int k=0, last=0; k<n; ++k) {
              SEXP s = levelsD[k];
              const int tl = TRUELENGTH(s);
              if (tl>=last) {  // if tl>=0 then also tl>=last because last<=0
                if (tl>=0) {
                  sprintf(warnStr,    // not direct warning as we're inside tl region
                  "Column %d of item %d is an ordered factor but level %d ['%s'] is missing from the ordered levels from column %d of item %d. " \
                  "Each set of ordered factor levels should be an ordered subset of the first longest. A regular factor will be created for this column.",
                  w+1, i+1, k+1, CHAR(s), longestW+1, longestI+1);
                } else {
                  sprintf(warnStr,
                  "Column %d of item %d is an ordered factor with '%s'<'%s' in its levels. But '%s'<'%s' in the ordered levels from column %d of item %d. " \
                  "A regular factor will be created for this column due to this ambiguity.",
                  w+1, i+1, CHAR(levelsD[k-1]), CHAR(s), CHAR(s), CHAR(levelsD[k-1]), longestW+1, longestI+1);
                  // k>=1 (so k-1 is ok) because when k==0 last==0 and this branch wouldn't happen
                }
                orderedFactor=false;
                i=LENGTH(l);  // break outer i loop
                break;        // break inner k loop
                // we leave the tl set for the longest levels; the regular factor will be created with the longest ordered levels first in case that useful for user
              }
              last = tl;  // negative ordinal; last should monotonically grow more negative if the levels are an ordered subset of the longest
            }
          }
        }
      }
      for (int i=0; i<LENGTH(l); ++i) {
        const int thisnrow = eachMax[i];
        SEXP li = VECTOR_ELT(l, i);
        if (!length(li)) continue;  // NULL items in the list() of DT/DF; not if thisnrow==0 because we need to retain (unused) factor levels (#3508)
        int w = usenames ? colMap[i*ncol + j] : j;
        if (w==-1) {
          writeNA(target, ansloc, thisnrow);
        } else {
          SEXP thisCol = VECTOR_ELT(li, w);
          SEXP thisColStr = isFactor(thisCol) ? getAttrib(thisCol, R_LevelsSymbol) : (isString(thisCol) ? thisCol : VECTOR_ELT(coercedForFactor, i));
          const int n = length(thisColStr);
          const SEXP *thisColStrD = STRING_PTR(thisColStr);  // D for data
          for (int k=0; k<n; ++k) {
            SEXP s = thisColStrD[k];
            if (s==NA_STRING ||             // remove NA from levels; test 1979 found by package emil when revdep testing 1.12.2 (#3473)
                TRUELENGTH(s)<0) continue;  // seen this level before; handles removing dups from levels as well as finding unique of character columns
            if (TRUELENGTH(s)>0) savetl(s);
            if (allocLevel==nLevel) {       // including initial time when allocLevel==nLevel==0
              SEXP *tt = NULL;
              if (allocLevel<INT_MAX) {
                int64_t new = (int64_t)allocLevel+n-k+1024; // if all remaining levels in this item haven't been seen before, plus 1024 margin in case of many very short levels
                allocLevel = (new>(int64_t)INT_MAX) ? INT_MAX : (int)new;
                tt = (SEXP *)realloc(levelsRaw, allocLevel*sizeof(SEXP));  // first time levelsRaw==NULL and realloc==malloc in that case
              }
              if (tt==NULL) {
                // # nocov start
                // C spec states that if realloc() fails (above) the original block (levelsRaw) is left untouched: it is not freed or moved. We ...
                for (int k=0; k<nLevel; k++) SET_TRUELENGTH(levelsRaw[k], 0);   // ... rely on that in this loop which uses levelsRaw.
                free(levelsRaw);
                savetl_end();
                error("Failed to allocate working memory for %d factor levels of result column %d when reading item %d of item %d", allocLevel, idcol+j+1, w+1, i+1);
                // # nocov end
              }
              levelsRaw = tt;
            }
            SET_TRUELENGTH(s,-(++nLevel));
            levelsRaw[nLevel-1] = s;
          }
          int *targetd = INTEGER(target);
          if (isFactor(thisCol)) {
            const int *id = INTEGER(thisCol);
            if (length(thisCol)<=1) {
              // recycle length-1, or NA-fill length-0
              SEXP lev;
              const int val = (length(thisCol)==1 && id[0]!=NA_INTEGER && (lev=thisColStrD[id[0]-1])!=NA_STRING) ? -TRUELENGTH(lev) : NA_INTEGER;
              //                                                                                    ^^ #3915 and tests 2015.2-5
              for (int r=0; r<thisnrow; ++r) targetd[ansloc+r] = val;
            } else {
              // length(thisCol)==thisnrow alreay checked before this truelength-clobber region
              // If all i==truelength(i) then just do a memcpy since hop is identity. Otherwise hop via the integer map.
              bool hop = false;
              if (orderedFactor) {
                // retain the position of NA level (if any) and the integer mappings to it
                for (int k=0; k<n; ++k) {
                  SEXP s = thisColStrD[k];
                  if (s!=NA_STRING && -TRUELENGTH(s)!=k+1) { hop=true; break; }
                }
              } else {
                for (int k=0; k<n; ++k) {
                  SEXP s = thisColStrD[k];
                  if (s==NA_STRING || -TRUELENGTH(s)!=k+1) { hop=true; break; }
                }
              }
              if (hop) {
                if (orderedFactor) {
                  for (int r=0; r<thisnrow; ++r)
                    targetd[ansloc+r] = id[r]==NA_INTEGER ? NA_INTEGER : -TRUELENGTH(thisColStrD[id[r]-1]);
                } else {
                  for (int r=0; r<thisnrow; ++r) {
                    SEXP lev;
                    targetd[ansloc+r] = id[r]==NA_INTEGER || (lev=thisColStrD[id[r]-1])==NA_STRING ? NA_INTEGER : -TRUELENGTH(lev);
                  }
                }
              } else {
                memcpy(targetd+ansloc, id, thisnrow*SIZEOF(thisCol));
              }
            }
          } else {
            const SEXP *sd = STRING_PTR(thisColStr);
            if (length(thisCol)<=1) {
              const int val = (length(thisCol)==1 && sd[0]!=NA_STRING) ? -TRUELENGTH(sd[0]) : NA_INTEGER;
              for (int r=0; r<thisnrow; ++r) targetd[ansloc+r] = val;
            } else {
              for (int r=0; r<thisnrow; ++r) targetd[ansloc+r] = sd[r]==NA_STRING ? NA_INTEGER : -TRUELENGTH(sd[r]);
            }
          }
        }
        ansloc += thisnrow;
      }
      for (int k=0; k<nLevel; ++k) SET_TRUELENGTH(levelsRaw[k], 0);
      savetl_end();
      if (warnStr[0]) warning(warnStr);  // now savetl_end() has happened it's safe to call warning (could error if options(warn=2))
      SEXP levelsSxp;
      setAttrib(target, R_LevelsSymbol, levelsSxp=allocVector(STRSXP, nLevel));
      for (int k=0; k<nLevel; ++k) SET_STRING_ELT(levelsSxp, k, levelsRaw[k]);
      free(levelsRaw);
      if (orderedFactor) {
        SEXP tt;
        setAttrib(target, R_ClassSymbol, tt=allocVector(STRSXP, 2));
        SET_STRING_ELT(tt, 0, char_ordered); // PRINTNAME(install("ordered"))
        SET_STRING_ELT(tt, 1, char_factor); //  PRINTNAME(install("factor"))
      } else {
        setAttrib(target, R_ClassSymbol, ScalarString(char_factor)); // "factor"
      }
    } else {  // factor==false
      for (int i=0; i<LENGTH(l); ++i) {
        const int thisnrow = eachMax[i];
        if (thisnrow==0) continue;
        SEXP li = VECTOR_ELT(l, i);
        int w = usenames ? colMap[i*ncol + j] : j;
        SEXP thisCol;
        if (w==-1 || !length(thisCol=VECTOR_ELT(li, w))) {  // !length for zeroCol warning above; #1871
          writeNA(target, ansloc, thisnrow);  // writeNA is integer64 aware and writes INT64_MIN
        } else {
          if (TYPEOF(target)==VECSXP && TYPEOF(thisCol)!=VECSXP) {
            // do an as.list() on the atomic column; #3528
            thisCol = PROTECT(coerceVector(thisCol, VECSXP)); nprotect++;
          }
          // else coerces if needed within memrecycle; with a no-alloc direct coerce from 1.12.4 (PR #3909)
          const char *ret = memrecycle(target, R_NilValue, ansloc, thisnrow, thisCol, idcol+j+1, foundName);
          if (ret) warning("Column %d of item %d: %s", w+1, i+1, ret);
          // e.g. when precision is lost like assigning 3.4 to integer64; test 2007.2
          // TODO: but maxType should handle that and this should never warn
        }
        ansloc += thisnrow;
      }
    }
  }
  UNPROTECT(nprotect);  // ans, coercedForFactor, thisCol
  return(ans);
}


