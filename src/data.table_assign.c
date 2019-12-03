#include "data.table.h"
#include <Rdefines.h>
#include <Rmath.h>


/* There are two reasons the finalizer doesn't restore the LENGTH to TRUELENGTH. i) The finalizer
happens too late after GC has already released the memory, and ii) copies by base R (e.g.
[<- in write.table just before test 894) allocate at length LENGTH but copy the TRUELENGTH over.
If the finalizer sets LENGTH to TRUELENGTH, that's a fail as it wasn't really allocated at
TRUELENGTH when R did the copy.
Karl Miller suggested an ENVSXP so that restoring LENGTH in finalizer should work. This is the
closest I got to getting it to pass all tests :

  SEXP env = PROTECT(allocSExp(ENVSXP));
  defineVar(SelfRefSymbol, x, env);
  defineVar(R_NamesSymbol, getAttrib(x, R_NamesSymbol), env);
  setAttrib(x, SelfRefSymbol, p = R_MakeExternalPtr(
    R_NilValue,         // for identical() to return TRUE. identical() doesn't look at tag and prot
    R_NilValue, //getAttrib(x, R_NamesSymbol), // to detect if names has been replaced and its tl lost, e.g. setattr(DT,"names",...)
    PROTECT(            // needed when --enable-strict-barrier it seems, iiuc. TO DO: test under that flag and remove if not needed.
      env               // wrap x in env to avoid an infinite loop in object.size() if prot=x were here
    )
  ));
  R_RegisterCFinalizerEx(p, finalizer, FALSE);
  UNPROTECT(2);

Then in finalizer:
  SETLENGTH(names, tl)
  SETLENGTH(dt, tl)

and that finalizer indeed now happens before the GC releases memory (thanks to the env wrapper).

However, we still have problem (ii) above and it didn't pass tests involving base R copies.

We really need R itself to start setting TRUELENGTH to be the allocated length and then
for GC to release TRUELENGTH not LENGTH.  Would really tidy this up.

Moved out of ?setkey Details section in 1.12.2 (Mar 2019). Revisit this w.r.t. to recent versions of R.
  The problem (for \code{data.table}) with the copy by \code{key<-} (other than
  being slower) is that \R doesn't maintain the over-allocated truelength, but it
  looks as though it has. Adding a column by reference using \code{:=} after a
  \code{key<-} was therefore a memory overwrite and eventually a segfault; the
  over-allocated memory wasn't really there after \code{key<-}'s copy. \code{data.table}s now have an attribute \code{.internal.selfref} to catch and warn about such copies.
  This attribute has been implemented in a way that is friendly with
  \code{identical()} and \code{object.size()}.
*/

// line 227:
// int checkOverAlloc(SEXP x)
// {
//  if (isNull(x))
//    error("Has getOption('datatable.alloccol') somehow become unset? It should be a number, by default 1024.");
//  if (!isInteger(x) && !isReal(x))
//    error("getOption('datatable.alloccol') should be a number, by default 1024. But its type is '%s'.", type2char(TYPEOF(x)));
//  if (LENGTH(x) != 1)
//    error("getOption('datatable.alloc') is a numeric vector ok but its length is %d. Its length should be 1.", LENGTH(x));
//  int ans = isInteger(x) ? INTEGER(x)[0] : (int)REAL(x)[0];
//  if (ans<0)
//    error("getOption('datatable.alloc')==%d.  It must be >=0 and not NA.", ans);
//  return ans;
// }


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


int *_Last_updated = NULL;

static bool anyNamed(SEXP x) {
  if (MAYBE_REFERENCED(x)) return true;
  if (isNewList(x)) for (int i=0; i<LENGTH(x); i++)
    if (anyNamed(VECTOR_ELT(x,i))) return true;
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
  } else if ((TYPEOF(target)!=TYPEOF(source) || targetIsI64!=sourceIsI64) && !isNewList(target)) {
    // The following checks are up front here, otherwise we'd need them twice in the two branches
    //   inside BODY that cater for 'where' or not. Maybe there's a way to merge the two macros in future.
    // The idea is to do these range checks without calling coerceVector() (which allocates)

#define CHECK_RANGE(STYPE, RFUN, COND, FMT, TO) {{                                                                                \
    const STYPE *sd = (const STYPE *)RFUN(source);                                                                                \
    for (int i=0; i<slen; ++i) {                                                                                                  \
      const STYPE val = sd[i];                                                                                                    \
      if (COND) {                                                                                                                 \
        const char *sType = sourceIsI64 ? "integer64" : type2char(TYPEOF(source));                                                \
        const char *tType = targetIsI64 ? "integer64" : type2char(TYPEOF(target));                                                \
        int n = snprintf(memrecycle_message, MSGSIZE,                                                                             \
                         FMT" (type '%s') at RHS position %d "TO" when assigning to type '%s'", val, sType, i+1, tType);          \
        if (colnum>0 && n>0 && n<MSGSIZE)                                                                                         \
          snprintf(memrecycle_message+n, MSGSIZE-n, " (column %d named '%s')", colnum, colname);                                  \
        /* string returned so that rbindlist/dogroups can prefix it with which item of its list this refers to  */                \
        break;                                                                                                                    \
      }                                                                                                                           \
    }                                                                                                                             \
    } break; }

  switch(TYPEOF(target)) {
  case LGLSXP:
    switch (TYPEOF(source)) {
    case RAWSXP:  CHECK_RANGE(Rbyte, RAW,      val!=0 && val!=1,                                        "%d",   "taken as TRUE")
    case INTSXP:  CHECK_RANGE(int, INTEGER,    val!=0 && val!=1 && val!=NA_INTEGER,                     "%d",   "taken as TRUE")
    case REALSXP: if (sourceIsI64)
      CHECK_RANGE(long long, REAL, val!=0 && val!=1 && val!=NA_INTEGER64,                   "%lld", "taken as TRUE")
      else  CHECK_RANGE(double, REAL,    !ISNAN(val) && val!=0.0 && val!=1.0,                     "%f",   "taken as TRUE")
    } break;
  case RAWSXP:
    switch (TYPEOF(source)) {
    case INTSXP:  CHECK_RANGE(int, INTEGER,    val<0 || val>255,                                        "%d",   "taken as 0")
    case REALSXP: if (sourceIsI64)
      CHECK_RANGE(long long, REAL, val<0 || val>255,                                        "%lld", "taken as 0")
      else  CHECK_RANGE(double, REAL,    !R_FINITE(val) || val<0.0 || val>256.0 || (int)val!=val, "%f",   "either truncated (precision lost) or taken as 0")
    } break;
  case INTSXP:
    if (TYPEOF(source)==REALSXP) {
      if (sourceIsI64)
        CHECK_RANGE(long long, REAL, val!=NA_INTEGER64 && (val<=NA_INTEGER || val>INT_MAX),   "%lld",  "out-of-range (NA)")
        else        CHECK_RANGE(double, REAL,    !ISNAN(val) && (!R_FINITE(val) || (int)val!=val),        "%f",    "truncated (precision lost)")
    } break;
  case REALSXP:
    if (targetIsI64 && isReal(source) && !sourceIsI64) {
      CHECK_RANGE(double, REAL,    !ISNAN(val) && (!R_FINITE(val) || (int)val!=val),        "%f",    "truncated (precision lost)")
    }
  }
  }

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
      BODY(SEXP, VECTOR_PTR, SEXP, val,  SET_VECTOR_ELT(target, off+i, cval))
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


SEXP setcolorder(SEXP x, SEXP o)
{
  SEXP names = getAttrib(x, R_NamesSymbol);
  const int *od = INTEGER(o), ncol=LENGTH(x);
  if (isNull(names)) error("dt passed to setcolorder has no names");
  if (ncol != LENGTH(names))
    error("Internal error: dt passed to setcolorder has %d columns but %d names", ncol, LENGTH(names));  // # nocov

  // Double-check here at C level that o[] is a strict permutation of 1:ncol. Reordering columns by reference makes no
  // difference to generations/refcnt so we can write behind barrier in this very special case of strict permutation.
  bool *seen = Calloc(ncol, bool);
  for (int i=0; i<ncol; ++i) {
    if (od[i]==NA_INTEGER || od[i]<1 || od[i]>ncol)
      error("Internal error: o passed to Csetcolorder contains an NA or out-of-bounds");  // # nocov
    if (seen[od[i]-1])
      error("Internal error: o passed to Csetcolorder contains a duplicate");             // # nocov
    seen[od[i]-1] = true;
  }
  Free(seen);

  SEXP *tmp = Calloc(ncol, SEXP);
  SEXP *xd = VECTOR_PTR(x), *namesd = STRING_PTR(names);
  for (int i=0; i<ncol; ++i) tmp[i] = xd[od[i]-1];
  memcpy(xd, tmp, ncol*sizeof(SEXP)); // sizeof is type size_t so no overflow here
  for (int i=0; i<ncol; ++i) tmp[i] = namesd[od[i]-1];
  memcpy(namesd, tmp, ncol*sizeof(SEXP));
  // No need to change key (if any); sorted attribute is column names not positions
  Free(tmp);
  return(R_NilValue);
}
