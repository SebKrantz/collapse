#ifdef _OPENMP
#include <omp.h>
#endif
#include "collapse_c.h"

SEXP Cna_rm(SEXP x) {
  const int n = LENGTH(x);
  if (n < 1) return x;
  int k = 0;
  switch(TYPEOF(x)) {
  case LGLSXP:
  case INTSXP: {
    const int *xd = INTEGER(x);
    for (int i = 0; i != n; ++i) if(xd[i] == NA_INTEGER) ++k;
    if(k == 0) return x;
    SEXP out = PROTECT(allocVector(TYPEOF(x), n - k));
    int *pout = INTEGER(out);
    k = 0;
    for (int i = 0; i != n; ++i) if(xd[i] != NA_INTEGER) pout[k++] = xd[i];
    copyMostAttrib(x, out);
    UNPROTECT(1);
    return out;
  }
  case REALSXP: { // What about integer64??
    const double *xd = REAL(x);
    for (int i = 0; i != n; ++i) if(ISNAN(xd[i])) ++k;
    if(k == 0) return x;
    SEXP out = PROTECT(allocVector(REALSXP, n - k));
    double *pout = REAL(out);
    k = 0;
    for (int i = 0; i != n; ++i) if(NISNAN(xd[i])) pout[k++] = xd[i]; // using xd[i] == xd[i] is not faster !!
    copyMostAttrib(x, out);
    UNPROTECT(1);
    return out;
  }
  case STRSXP: {
    const SEXP *xd = STRING_PTR(x);
    for (int i = 0; i != n; ++i) if(xd[i] == NA_STRING) ++k;
    if(k == 0) return x;
    SEXP out = PROTECT(allocVector(STRSXP, n - k));
    SEXP *pout = STRING_PTR(out);
    k = 0;
    for (int i = 0; i != n; ++i) if(xd[i] != NA_STRING) pout[k++] = xd[i];
    copyMostAttrib(x, out);
    UNPROTECT(1);
    return out;
  }
  case VECSXP: {
    const SEXP *xd = SEXPPTR(x);
    for (int i = 0; i != n; ++i) if(length(xd[i]) == 0) ++k;
    if(k == 0) return x;
    SEXP out = PROTECT(allocVector(VECSXP, n - k));
    SEXP *pout = SEXPPTR(out);
    k = 0;
    for (int i = 0; i != n; ++i) if(length(xd[i]) != 0) pout[k++] = xd[i];
    copyMostAttrib(x, out);
    UNPROTECT(1);
    return out;
  }
  }
  error("Unsupported type '%s' passed to na_rm()", type2char(TYPEOF(x)));
}

// Helper function to find a single sting in factor levels
int fchmatch(SEXP x, SEXP val, int nomatch) {
  const SEXP *px = STRING_PTR(x), v = asChar(val);
  for(int i = 0, l = length(x); i != l; ++i) if(px[i] == v) return i + 1;
  return nomatch;
}

SEXP whichv(SEXP x, SEXP val, SEXP Rinvert) {

  int j = 0, n = length(x), invert = asLogical(Rinvert);
  int *buf = (int *) R_alloc(n, sizeof(int));
  SEXP ans;

#define WHICHVLOOP                                               \
  if(invert) {                                                   \
    for(int i = 0; i != n; ++i) if(px[i] != v) buf[j++] = i+1;   \
  } else {                                                       \
    for(int i = 0; i != n; ++i) if(px[i] == v) buf[j++] = i+1;   \
  }

#define WHICHVLOOPLX                                               \
if(invert) {                                                       \
  for(int i = 0; i != n; ++i) if(px[i] != pv[i]) buf[j++] = i+1;   \
} else {                                                           \
  for(int i = 0; i != n; ++i) if(px[i] == pv[i]) buf[j++] = i+1;   \
}

if(length(val) == n && n > 1) {
  if(TYPEOF(val) != TYPEOF(x)) error("data types of x and value must be the same");
  switch(TYPEOF(x)) {
  case INTSXP:
  case LGLSXP:
  {
    const int *px = INTEGER(x);
    const int *pv = INTEGER(val);
    WHICHVLOOPLX
      break;
  }
  case REALSXP:
  {
    const double *px = REAL(x);
    const double *pv = REAL(val);
    WHICHVLOOPLX
      break;
  }
  case STRSXP:
  {
    const SEXP *px = STRING_PTR(x);
    const SEXP *pv = STRING_PTR(val);
    WHICHVLOOPLX
      break;
  }
  case RAWSXP :
  {
    const Rbyte *px = RAW(x);
    const Rbyte *pv = RAW(val);
    WHICHVLOOPLX
      break;
  }
  default: error("Unsupported type '%s' passed to whichv()", type2char(TYPEOF(x)));
  }
} else {
  if(length(val) != 1) error("length(value) needs to be length(x) or 1");
  switch(TYPEOF(x)) {
  case INTSXP:
  case LGLSXP:
  {
    const int *px = INTEGER(x);
    int v;
    if(TYPEOF(val) == STRSXP) {
      if(!isFactor(x)) error("Type mismatch: if value is character, x must be character or factor.");
      v = fchmatch(getAttrib(x, R_LevelsSymbol), val, 0);
    } else v = asInteger(val);
    WHICHVLOOP
      break;
  }
  case REALSXP:
  {
    const double *px = REAL(x);
    const double v = asReal(val);
    if(ISNAN(v)) {
      if(invert) {
        for(int i = 0; i != n; ++i) if(NISNAN(px[i])) buf[j++] = i+1;
      } else {
        for(int i = 0; i != n; ++i) if(ISNAN(px[i])) buf[j++] = i+1;
      }
    } else {
      WHICHVLOOP
    }
    break;
  }
  case STRSXP:
  {
    const SEXP *px = STRING_PTR(x);
    const SEXP v = asChar(val);
    WHICHVLOOP
      break;
  }
  case RAWSXP :
  {
    const Rbyte *px = RAW(x);
    const Rbyte v = RAW(val)[0];
    WHICHVLOOP
      break;
  }
  default: error("Unsupported type '%s' passed to whichv()", type2char(TYPEOF(x)));
  }
}
PROTECT(ans = allocVector(INTSXP, j));
if(j) memcpy(INTEGER(ans), buf, sizeof(int) * j);

UNPROTECT(1);
return(ans);
}

SEXP anyallv(SEXP x, SEXP val, SEXP Rall) {

  int n = length(x), all = asLogical(Rall);
  if(length(val) != 1) error("value needs to be length 1");

#define ALLANYVLOOP                                                      \
  if(all) {                                                              \
    for(int i = 0; i != n; ++i) if(px[i] != v) return ScalarLogical(0);  \
    return ScalarLogical(1);                                             \
  } else {                                                               \
    for(int i = 0; i != n; ++i) if(px[i] == v) return ScalarLogical(1);  \
    return ScalarLogical(0);                                             \
  }

switch(TYPEOF(x)) {
case INTSXP:
case LGLSXP:
{
  const int *px = INTEGER(x);
  int v;
  if(TYPEOF(val) == STRSXP) {
    if(!isFactor(x)) error("Type mismatch: if value is character, x must be character or factor.");
    v = fchmatch(getAttrib(x, R_LevelsSymbol), val, 0);
  } else v = asInteger(val);
  ALLANYVLOOP
    break;
}
case REALSXP:
{
  const double *px = REAL(x);
  const double v = asReal(val);
  if(ISNAN(v)) error("please use allNA()");
  ALLANYVLOOP
    break;
}
case STRSXP:
{
  const SEXP *px = STRING_PTR(x);
  const SEXP v = asChar(val);
  ALLANYVLOOP
    break;
}
case RAWSXP :
{
  const Rbyte *px = RAW(x);
  const Rbyte v = RAW(val)[0];
  ALLANYVLOOP
    break;
}
default: error("Unsupported type '%s' passed to allv() / anyv()", type2char(TYPEOF(x)));
}
  return(R_NilValue);
}

SEXP setcopyv(SEXP x, SEXP val, SEXP rep, SEXP Rinvert, SEXP Rset, SEXP Rind1) {

  const int n = length(x), lv = length(val), lr = length(rep),
    tx = TYPEOF(x), ind1 = asLogical(Rind1), invert = asLogical(Rinvert), set = asLogical(Rset);
  int nprotect = 0, tv = TYPEOF(val), tr = TYPEOF(rep);

  if(lv > 1 || ind1) {
    if(tv == LGLSXP) {
      if(lv != n) error("If v is a logical vector, length(v) needs to be equal to length(x)");
      if(lr != 1 && lr != n) error("If v is a logical vector, length(r) needs to be 1 or length(x)");
    } else if(tv == INTSXP || tv == REALSXP) {
      if(invert) error("invert = TRUE is only possible if v is a logical vector");
      if(lv == 0) return x; // integer(0) cannot cause error
      if(lv > n) error("length(v) must be <= length(x)");
      if(!(lr == 1 || lr == n || lr == lv)) error("length(r) must be either 1, length(v) or length(x)");
      if(tv == REALSXP) {
        if(lv == 1 && REAL_ELT(val, 0) == (int)REAL_ELT(val, 0)) {
          tv = INTSXP;
          val = PROTECT(coerceVector(val, INTSXP));
          ++nprotect;
        } else error("If length(v) > 1 or vind1 = TRUE, v must be an integer or logical vector");
      }
      // Just some heuristic checking as this is a programmers function
      const int v1 = INTEGER_ELT(val, 0), vn = INTEGER_ELT(val, lv-1);
      if(v1 < 1 || v1 > n || vn < 1 || vn > n) error("Detected index (v) outside of range [1, length(x)]");
    } else error("If length(v) > 1 or vind1 = TRUE, v must be an integer or logical vector");
  } else {
    if(lv == 0) return x; // empty replacement, good to return?
    if(lr != 1 && lr != n) error("If length(v) == 1, length(r) must be 1 or length(x)");
  }

  if(lr > 1 && tr != tx) { // lr == n &&
    if(!((tx == INTSXP && tr == LGLSXP) || (tx == LGLSXP && tr == INTSXP))) {
      // PROTECT_INDEX ipx;
      // PROTECT_WITH_INDEX(rep = coerceVector(rep, tx), &ipx);
      tr = tx;
      rep = PROTECT(coerceVector(rep, tx));
      ++nprotect;
    } // error("typeof(x) needs to match typeof(r)");
  }


  SEXP ans = R_NilValue;
  if(set == 0) {
    PROTECT(ans = shallow_duplicate(x)); // Fastest?? // copies attributes ?? -> Yes
    ++nprotect;
  }

  #define setcopyvLOOP(e)                                     \
  if(invert) {                                                \
    for(int i = 0; i != n; ++i) if(px[i] != v) px[i] = e;     \
  } else {                                                    \
    for(int i = 0; i != n; ++i) if(px[i] == v) px[i] = e;     \
  }

  #define setcopyvLOOPLVEC1                                     \
  if(tv == INTSXP) {                                            \
    for(int i = 0; i != lv; ++i) px[pv[i]-1] = r;               \
  } else if(invert == 0) {                                      \
    for(int i = 0; i != n; ++i) if(pv[i] > 0) px[i] = r;        \
  } else {                                                      \
    for(int i = 0; i != n; ++i) if(pv[i] == 0) px[i] = r;       \
  }

  #define setcopyvLOOPLVEC                                      \
  if(tv == INTSXP) {                                            \
    if(lr == n) {                                               \
      for(int i = 0; i != lv; ++i) px[pv[i]-1] = pr[pv[i]-1];   \
    } else {                                                    \
      for(int i = 0; i != lv; ++i) px[pv[i]-1] = pr[i];         \
    }                                                           \
  } else if(invert == 0) {                                      \
    for(int i = 0; i != n; ++i) if(pv[i] > 0) px[i] = pr[i];    \
  } else {                                                      \
    for(int i = 0; i != n; ++i) if(pv[i] == 0) px[i] = pr[i];   \
  }

  switch(tx) {
  case INTSXP:
  case LGLSXP:
  {
    int *restrict px = set ? INTEGER(x) : INTEGER(ans);
    if(lv == 1 && ind1 == 0) {
      int v;
      if(tv == STRSXP) {
        if(!isFactor(x)) error("Type mismatch: if v is character, x must be character or factor.");
        v = fchmatch(getAttrib(x, R_LevelsSymbol), val, 0);
      } else v = asInteger(val);
      if(lr == 1) {
        const int r = asInteger(rep);
        setcopyvLOOP(r)
      } else {
        const int *restrict pr = INTEGER(rep);
        setcopyvLOOP(pr[i])
      }
    } else {
      const int *restrict pv = INTEGER(val); // ALTREP(val) ? (const int *)ALTVEC_DATAPTR(val) :
      if(lr == 1) {
        const int r = asInteger(rep);
        setcopyvLOOPLVEC1
      } else {
        const int *restrict pr = INTEGER(rep);
        setcopyvLOOPLVEC
      }
    }
    break;
  }
  case REALSXP:
  {
    double *restrict px = set ? REAL(x) : REAL(ans);
    if(lv == 1 && ind1 == 0) {
      const double v = asReal(val);
      if(lr == 1) {
        const double r = asReal(rep);
        if(ISNAN(v)) {
          if(invert) {
            for(int i = 0; i != n; ++i) if(NISNAN(px[i])) px[i] = r;
          } else {
            for(int i = 0; i != n; ++i) if(ISNAN(px[i])) px[i] = r;
          }
        } else {
          setcopyvLOOP(r)
        }
      } else {
        const double *restrict pr = REAL(rep);
        if(ISNAN(v)) {
          if(invert) {
            for(int i = 0; i != n; ++i) if(NISNAN(px[i])) px[i] = pr[i];
          } else {
            for(int i = 0; i != n; ++i) if(ISNAN(px[i])) px[i] = pr[i];
          }
        } else {
          setcopyvLOOP(pr[i])
        }
      }
    } else {
      const int *restrict pv = INTEGER(val); // ALTREP(val) ? (const int *)ALTVEC_DATAPTR(val) :
      if(lr == 1) {
        const double r = asReal(rep);
        setcopyvLOOPLVEC1
      } else {
        const double *restrict pr = REAL(rep);
        setcopyvLOOPLVEC
      }
    }
    break;
  }
  case STRSXP:
  {
    SEXP *restrict px = set ? STRING_PTR(x) : STRING_PTR(ans);
    if(lv == 1 && ind1 == 0) {
      const SEXP v = PROTECT(asChar(val));
      if(lr == 1) {
        const SEXP r = asChar(rep);
        setcopyvLOOP(r)
      } else {
        const SEXP *restrict pr = STRING_PTR(rep);
        setcopyvLOOP(pr[i])
      }
      UNPROTECT(1);
    } else {
      const int *restrict pv = INTEGER(val); // ALTREP(val) ? (const int *)ALTVEC_DATAPTR(val) :
      if(lr == 1) {
        const SEXP r = asChar(rep);
        setcopyvLOOPLVEC1
      } else {
        const SEXP *restrict pr = STRING_PTR(rep);
        setcopyvLOOPLVEC
      }
    }
    break;
  }
  case VECSXP:
  {
    SEXP *restrict px = set ? SEXPPTR(x) : SEXPPTR(ans);
    if(lv == 1 && ind1 == 0) error("Cannot compare lists to a value");
    // if(tr != VECSXP) error("If X is a list and xlist = TRUE, R also needs to be a list");
    const int *restrict pv = INTEGER(val); // ALTREP(val) ? (const int *)ALTVEC_DATAPTR(val) :
    if(lr == 1) {
      const SEXP r = VECTOR_ELT(rep, 0);
      setcopyvLOOPLVEC1
    } else {
      const SEXP *restrict pr = SEXPPTR(rep);
      setcopyvLOOPLVEC
    }
    break;
  }
  case RAWSXP:
  {
    Rbyte *restrict px = set ? RAW(x) : RAW(ans);
    if(lv == 1 && ind1 == 0) {
      const Rbyte v = RAW(val)[0];
      if(lr == 1) {
        const Rbyte r = RAW(rep)[0];
        setcopyvLOOP(r)
      } else {
        const Rbyte *restrict pr = RAW(rep);
        setcopyvLOOP(pr[i])
      }
    } else {
      const int *restrict pv = INTEGER(val); // ALTREP(val) ? (const int *)ALTVEC_DATAPTR(val) :
      if(lr == 1) {
        const Rbyte r = RAW(rep)[0];
        setcopyvLOOPLVEC1
      } else {
        const Rbyte *restrict pr = RAW(rep);
        setcopyvLOOPLVEC
      }
    }
    break;
  }
  default: error("Unsupported type '%s' passed to setv() / copyv()", type2char(tx));
  }

  UNPROTECT(nprotect);
  if(set == 0) return(ans);
  return(x);
}

SEXP setop_core(SEXP x, SEXP val, SEXP op, SEXP roww) {

  int n = length(x), nv = length(val), o = asInteger(op), tx = TYPEOF(x);

#define OPSWITCH(e)                                  \
  switch(o) {                                        \
  case 1: for(int i = 0; i != n; ++i) px[i] += e;    \
    break;                                           \
  case 2: for(int i = 0; i != n; ++i) px[i] -= e;    \
    break;                                           \
  case 3: for(int i = 0; i != n; ++i) px[i] *= e;    \
    break;                                           \
  case 4: for(int i = 0; i != n; ++i) px[i] /= e;    \
    break;                                           \
  default: error("unsupported operation");           \
  }

if(nv == 1 || nv == n) {
  switch(tx) {
  case INTSXP:
  case LGLSXP:
  {
    int *px = INTEGER(x);
    if(nv == 1) {
      const int v = asInteger(val);
      OPSWITCH(v)
    } else {
      if(TYPEOF(val) == REALSXP) {
        // warning("adding real values to an integer: will truncate decimals");
        const double *v = REAL(val);
        OPSWITCH(v[i])
      } else {
        const int *v = INTEGER(val);
        OPSWITCH(v[i])
      }
    }
    break;
  }
  case REALSXP:
  {
    double *px = REAL(x);
    if(nv == 1) {
      const double v = asReal(val);
      OPSWITCH(v)
    } else {
      if(TYPEOF(val) == REALSXP) {
        const double *v = REAL(val);
        OPSWITCH(v[i])
      } else {
        const int *v = INTEGER(val);
        OPSWITCH(v[i])
      }
    }
    break;
  }
  default: error("Unsupported type '%s'", type2char(tx));
  }
} else {
  if(!isMatrix(x)) error("unequal argument lengths");
  int nr = nrows(x), nc = n / nr, rwl = asLogical(roww);
  if((rwl == 0 && nr != nv) || (rwl && nc != nv))
    error("length of vector must match matrix rows/columns or the size of the matrix itself");

#define OPSWITCHMAT(e)                               \
  switch(o) {                                        \
  case 1: for(int j = 0, cj; j != nc; ++j)  {        \
    cj = j * nr;                                     \
    for(int i = 0; i != nr; ++i) px[cj + i] += e;    \
  }                                                  \
  break;                                             \
  case 2: for(int j = 0, cj; j != nc; ++j)  {        \
    cj = j * nr;                                     \
    for(int i = 0; i != nr; ++i) px[cj + i] -= e;    \
  }                                                  \
  break;                                             \
  case 3: for(int j = 0, cj; j != nc; ++j)  {        \
    cj = j * nr;                                     \
    for(int i = 0; i != nr; ++i) px[cj + i] *= e;    \
  }                                                  \
  break;                                             \
  case 4: for(int j = 0, cj; j != nc; ++j)  {        \
    cj = j * nr;                                     \
    for(int i = 0; i != nr; ++i) px[cj + i] /= e;    \
  }                                                  \
  break;                                             \
  default: error("unsupported operation");           \
  }

switch(tx) {
case INTSXP:
case LGLSXP:
{
  int *px = INTEGER(x);
  if(TYPEOF(val) == REALSXP) {
    // warning("adding real values to an integer: will truncate decimals");
    const double *v = REAL(val);
    if(rwl) {
      OPSWITCHMAT(v[j])
    } else {
      OPSWITCHMAT(v[i])
    }
  } else {
    const int *v = INTEGER(val);
    if(rwl) {
      OPSWITCHMAT(v[j])
    } else {
      OPSWITCHMAT(v[i])
    }
  }
  break;
}
case REALSXP:
{
  double *px = REAL(x);
  if(TYPEOF(val) == REALSXP) {
    const double *v = REAL(val);
    if(rwl) {
      OPSWITCHMAT(v[j])
    } else {
      OPSWITCHMAT(v[i])
    }
  } else {
    const int *v = INTEGER(val);
    if(rwl) {
      OPSWITCHMAT(v[j])
    } else {
      OPSWITCHMAT(v[i])
    }
  }
  break;
}
default: error("Unsupported type '%s'", type2char(tx));
}

}
return(x);
}

SEXP setop(SEXP x, SEXP val, SEXP op, SEXP roww) {
  // IF x is a list, call function repeatedly..
  if(TYPEOF(x) == VECSXP) {
    SEXP *px = SEXPPTR(x);
    int lx = length(x);
    if(TYPEOF(val) == VECSXP) { // val is list: must match length(x)
      SEXP *pv = SEXPPTR(val);
      if(lx != length(val)) error("length(X) must match length(V)");
      for(int i = 0; i != lx; ++i) setop_core(px[i], pv[i], op, roww);
    } else if (length(val) == 1 || asLogical(roww) == 0) { // val is a scalar or vector but rowwise = FALSE
      for(int i = 0; i != lx; ++i) setop_core(px[i], val, op, roww);
    } else { // val is a numeric or logical vector to be applied rowwise
      if(lx != length(val)) error("length(X) must match length(V)");
      switch(TYPEOF(val)) {
      case REALSXP: {
        double *pv = REAL(val);
        for(int i = 0; i != lx; ++i) {
          setop_core(px[i], PROTECT(ScalarReal(pv[i])), op, roww); UNPROTECT(1);
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *pv = INTEGER(val);
        for(int i = 0; i != lx; ++i) {
          setop_core(px[i], PROTECT(ScalarInteger(pv[i])), op, roww); UNPROTECT(1);
        }
        break;
      }
      default: error("Unsupported type '%s'", type2char(TYPEOF(val)));
      }
    }
    return x;
  }
  return setop_core(x, val, op, roww);
}

SEXP vtypes(SEXP x, SEXP isnum) {
  int tx = TYPEOF(x);
  if(tx != VECSXP) return ScalarInteger(tx);
  SEXP *px = SEXPPTR(x); // This is ok, even if x contains ALTREP objects..
  int n = length(x);
  SEXP ans = PROTECT(allocVector(INTSXP, n));
  int *pans = INTEGER(ans);
  switch(asInteger(isnum)) {
  case 0:
    for(int i = 0; i != n; ++i) pans[i] = TYPEOF(px[i]) + 1;
    break;
  case 1: // Numeric variables: do_is with op = 100: https://github.com/wch/r-source/blob/2b0818a47199a0b64b6aa9b9f0e53a1e886e8e95/src/main/coerce.c
          // See also DispatchOrEval in https://github.com/wch/r-source/blob/trunk/src/main/eval.c
    {
    if(inherits(x, "indexed_frame")) { // NOT pdata.frame!! because columns in pdata.frame only become pseries when extracted from the frame
      for(int i = 0, tci, tnum, is_num; i != n; ++i) {
        is_num = 0;
        tci = TYPEOF(px[i]);
        tnum = tci == INTSXP || tci == REALSXP;
        if(tnum) is_num = inherits(px[i], "integer") || inherits(px[i], "numeric") || inherits(px[i], "ts") || inherits(px[i], "units");
        pans[i] = tnum && is_num;
      }
      // for(int i = 0; i != n; ++i) {
      //   int tci = TYPEOF(px[i]);
      //   pans[i] = (tci == INTSXP && inherits(px[i], "integer")) || (tci == REALSXP && inherits(px[i], "numeric")); // length(getAttrib(ci, R_ClassSymbol)) <= 2;
      // }
    } else {
      for(int i = 0, tci, tnum, is_num; i != n; ++i) {
        // pans[i] = isNumeric(px[i]) && !isLogical(px[i]); // Date is numeric, from: https://github.com/wch/r-source/blob/2b0818a47199a0b64b6aa9b9f0e53a1e886e8e95/src/main/coerce.c
        tci = TYPEOF(px[i]);
        tnum = tci == INTSXP || tci == REALSXP;
        is_num = tnum && OBJECT(px[i]) == 0;
        if(tnum && !is_num) is_num = inherits(px[i], "ts") || inherits(px[i], "units");
        pans[i] = is_num;
      }
    }
    SET_TYPEOF(ans, LGLSXP);
    break;
    }
  case 2: // is.factor
    for(int i = 0; i != n; ++i) pans[i] = (int)isFactor(px[i]);
    SET_TYPEOF(ans, LGLSXP);
    break;
  case 3: // is.list, needed for list processing functions
    for(int i = 0; i != n; ++i) pans[i] = TYPEOF(px[i]) == VECSXP;
    SET_TYPEOF(ans, LGLSXP);
    break;
  case 4: // is.sublist, needed for list processing functions
    for(int i = 0; i != n; ++i) pans[i] = TYPEOF(px[i]) == VECSXP && !isFrame(px[i]);
    SET_TYPEOF(ans, LGLSXP);
    break;
  case 7: // is.atomic(x), needed in atomic_elem()
    // is.atomic: do_is with op = 200:  https://github.com/wch/r-source/blob/9f9033e193071f256e21a181cb053cba983ed4a9/src/main/coerce.c
    for(int i = 0; i != n; ++i) {
      switch(TYPEOF(px[i])) {
      case NILSXP: /* NULL is atomic (S compatibly), but not in isVectorAtomic(.) */
      case CHARSXP:
      case LGLSXP:
      case INTSXP:
      case REALSXP:
      case CPLXSXP:
      case STRSXP:
      case RAWSXP:
        pans[i] = 1;
        break;
      default:
        pans[i] = 0;
      }
    }
    SET_TYPEOF(ans, LGLSXP);
    break;
  case 5: // is.atomic(x) || is.list(x), needed in reg_elem() and irreg_elem()
    for(int i = 0; i != n; ++i) {
      switch(TYPEOF(px[i])) {
      case VECSXP:
        pans[i] = 1;
        break;
      case NILSXP: /* NULL is atomic (S compatibly), but not in isVectorAtomic(.) */
      case CHARSXP:
      case LGLSXP:
      case INTSXP:
      case REALSXP:
      case CPLXSXP:
      case STRSXP:
      case RAWSXP:
        pans[i] = 1;
        break;
      default:
        pans[i] = 0;
      }
    }
    SET_TYPEOF(ans, LGLSXP);
    break;
  case 6:
    // Faster object type identification, needed in unlist2d:
    // idf <- function(x) if(inherits(x, "data.frame")) 2L else if (!length(x)) 1L else 3L*is.atomic(x)
    for(int i = 0; i != n; ++i) {
      if(length(px[i]) == 0)
        pans[i] = 1;
      else switch(TYPEOF(px[i])) {
           case VECSXP:
             pans[i] = isFrame(px[i]) ? 2 : 0;
             break;
           case NILSXP: /* NULL is atomic (S compatibly), but not in isVectorAtomic(.) */
           case CHARSXP:
           case LGLSXP:
           case INTSXP:
           case REALSXP:
           case CPLXSXP:
           case STRSXP:
           case RAWSXP:
            pans[i] = 3;
            break;
           default:
             pans[i] = 0;
           }
    }
    break;
  default:
    error("Unsupported vtypes option");
  }
  UNPROTECT(1);
  return ans;
}

SEXP vlengths(SEXP x, SEXP usenam) {
  if(TYPEOF(x) != VECSXP) return ScalarInteger(length(x));
  int n = length(x);
  SEXP ans = PROTECT(allocVector(INTSXP, n));
  int *pans = INTEGER(ans);
  if(ALTREP(x)) {
    for(int i = 0; i != n; ++i) pans[i] = length(VECTOR_ELT(x, i));
  } else {
    SEXP *px = SEXPPTR(x);
    for(int i = 0; i != n; ++i) pans[i] = length(px[i]);
  }
  if(asLogical(usenam)) {
    SEXP nam = getAttrib(x, R_NamesSymbol);
    if(TYPEOF(nam) != NILSXP) namesgets(ans, nam);
  }
  UNPROTECT(1);
  return ans;
}


// faster version of base::range, which calls both min() and max()
SEXP frange(SEXP x, SEXP Rnarm) {
  int l = length(x), narm = asLogical(Rnarm), tx = TYPEOF(x);

  SEXP out = PROTECT(allocVector(tx, 2));

  switch(tx) {
    case INTSXP:
    case LGLSXP:
    {
      if(l < 1) {
        INTEGER(out)[0] = INTEGER(out)[1] = NA_INTEGER;
        break;
      }
      int min, max, tmp, *px = INTEGER(x);
      if(narm) {
        int j = l-1;
        while(px[j] == NA_INTEGER && j!=0) --j;
        min = max = px[j];
        if(j != 0) for(int i = j; i--; ) {
          tmp = px[i];
          if(tmp == NA_INTEGER) continue;
          if(min > tmp) min = tmp;
          if(max < tmp) max = tmp;
        }
      } else {
        min = max = px[0];
        for(int i = 0; i != l; ++i) {
          tmp = px[i];
          if(tmp == NA_INTEGER) {
            min = max = tmp;
            break;
          } else {
            if(min > tmp) min = tmp;
            if(max < tmp) max = tmp;
          }
        }
      }
      INTEGER(out)[0] = min;
      INTEGER(out)[1] = max;
      break;
    }
    case REALSXP:
    {
      if(l < 1) {
        REAL(out)[0] = REAL(out)[1] = NA_REAL;
        break;
      }
      double min, max, tmp, *px = REAL(x);
      if(narm) {
        int j = l-1;
        while(ISNAN(px[j]) && j!=0) --j;
        min = max = px[j];
        if(j != 0) for(int i = j; i--; ) {
          tmp = px[i];
          if(min > tmp) min = tmp;
          if(max < tmp) max = tmp;
        }
      } else {
        min = max = px[0];
        for(int i = 0; i != l; ++i) {
          tmp = px[i];
          if(ISNAN(tmp)) {
            min = max = tmp;
            break;
          } else {
            if(min > tmp) min = tmp;
            if(max < tmp) max = tmp;
          }
        }
      }
      REAL(out)[0] = min;
      REAL(out)[1] = max;
      break;
    }
    default: error("Unsupported SEXP type: %s", type2char(tx));
  }

  copyMostAttrib(x, out);
  UNPROTECT(1);
  return out;
}


// faster distance matrices
// base R's version: https://github.com/wch/r-source/blob/79298c499218846d14500255efd622b5021c10ec/src/library/stats/src/distance.c
SEXP fdist(SEXP x, SEXP vec, SEXP Rret, SEXP Rnthreads) {

  SEXP dim = getAttrib(x, R_DimSymbol);
  int nrow, ncol, ret, nullv = isNull(vec), nthreads = asInteger(Rnthreads), nprotect = 1;
  if(nthreads > max_threads) nthreads = max_threads;
  if(TYPEOF(dim) != INTSXP) {
    nrow = 1;
    ncol= length(x);
  } else {
    nrow = INTEGER(dim)[0];
    ncol= INTEGER(dim)[1];
  }
  if(TYPEOF(x) != REALSXP) {
    x = PROTECT(coerceVector(x, REALSXP)); ++nprotect;
  }
  if(TYPEOF(Rret) == STRSXP) {
    const char *r = CHAR(STRING_ELT(Rret, 0));
    if(strcmp(r, "euclidean") == 0) ret = 1;
    else if(strcmp(r, "euclidean_squared") == 0) ret = 2;
    else error("Unsupported method: %s", r);
  } else {
    ret = asInteger(Rret);
    if(ret < 1 || ret > 2) error("method must be 1 ('euclidean') or 2 ('euclidean_squared')");
  }

  size_t l = nrow;
  if(nullv) { // Full distance matrix
    if(nrow <= 1) error("If v is left empty, x needs to be a matrix with at least 2 rows");
    l = ((double)nrow / 2) * (nrow - 1);
  } else if(length(vec) != ncol) error("length(v) must match ncol(x)");

  SEXP res = PROTECT(allocVector(REALSXP, l));
  double *px = REAL(x), *pres = REAL(res);
  memset(pres, 0, sizeof(double) * l);

  if(nullv) { // Full distance matrix
    if(nthreads > 1) {
      if(nthreads > nrow-1) nthreads = nrow-1;
      #pragma omp parallel for num_threads(nthreads)
      for (int k = 1; k < nrow; ++k) { // Row vectors to compute distances with
        int nmk = nrow - k;
        double *presk = pres + l - nmk*(nmk+1)/2, // https://en.wikipedia.org/wiki/1_%2B_2_%2B_3_%2B_4_%2B_%E2%8B%AF
               *pxj = px + k, v = pxj[-1], tmp;
        for (int j = 0; j < ncol; ++j) { // Elements of the row vector at hand
          for(int i = 0; i < nmk; ++i) { // All remaining rows to compute the distance to
            tmp = pxj[i] - v;
            presk[i] += tmp * tmp;
          }
          pxj += nrow; v = pxj[-1];
        }
      }
    } else {
      double *presk = pres, *pxj, v, tmp;
      for (int k = 1, nmk = nrow; k != nrow; ++k) { // Row vectors to compute distances with
        --nmk; pxj = px + k; v = pxj[-1];
        for (int j = 0; j != ncol; ++j) { // Elements of the row vector at hand
          for(int i = 0; i != nmk; ++i) { // All remaining rows to compute the distance to
            tmp = pxj[i] - v;
            presk[i] += tmp * tmp;
          }
          pxj += nrow; v = pxj[-1];
        }
        presk += nmk;
      }
    }
  } else { // Only a single vector
    if(TYPEOF(vec) != REALSXP) {
      vec = PROTECT(coerceVector(vec, REALSXP)); ++nprotect;
    }
    double *pv = REAL(vec);

    if(nrow > 1) { // x is a matrix
      if(nthreads > 1) {
        if(nthreads > nrow) nthreads = nrow;
        for (int j = 0; j < ncol; ++j) {
          double *pxj = px + j * nrow, v = pv[j];
          #pragma omp parallel for num_threads(nthreads)
          for (int i = 0; i < nrow; ++i) {
            double tmp = pxj[i] - v;
            pres[i] += tmp * tmp;
          }
        }
      } else {
        for (int j = 0; j != ncol; ++j) {
          double *pxj = px + j * nrow, v = pv[j], tmp;
          for (int i = 0; i != nrow; ++i) {
            tmp = pxj[i] - v;
            pres[i] += tmp * tmp;
          }
        }
      }
    } else { // x is a vector
      double dres = 0.0;
      if(nthreads > 1) {
        if(nthreads > ncol) nthreads = ncol;
        #pragma omp parallel for num_threads(nthreads) reduction(+:dres)
        for (int i = 0; i < ncol; ++i) {
          double tmp = px[i] - pv[i];
          dres += tmp * tmp;
        }
      } else {
        double tmp;
        for (int i = 0; i != ncol; ++i) {
          tmp = px[i] - pv[i];
          dres += tmp * tmp;
        }
      }
      pres[0] = ret == 1 ? sqrt(dres) : dres;
      ret = 2; // ensures we avoid the square root loop below
    }
  }

  // Square Root
  if(ret == 1) {
    if(nthreads > 1) {
      #pragma omp parallel for num_threads(nthreads)
      for (size_t i = 0; i < l; ++i) pres[i] = sqrt(pres[i]);
    } else {
      for (size_t i = 0; i != l; ++i) pres[i] = sqrt(pres[i]);
    }
  }

  if(nullv) { // Full distance matrix object
    // First creating symbols to avoid protect errors: https://blog.r-project.org/2019/04/18/common-protect-errors/
    SEXP sym_Size = install("Size"), sym_Labels = install("Labels"),
      sym_Diag = install("Diag"), sym_Upper = install("Upper"), sym_method = install("method");
    setAttrib(res, sym_Size, ScalarInteger(nrow));
    SEXP dn = getAttrib(x, R_DimNamesSymbol);
    if(TYPEOF(dn) == VECSXP && length(dn))
       setAttrib(res, sym_Labels, VECTOR_ELT(dn, 0));
    setAttrib(res, sym_Diag, ScalarLogical(0));
    setAttrib(res, sym_Upper, ScalarLogical(0));
    setAttrib(res, sym_method, mkString(ret == 1 ? "euclidean" : "euclidean_squared"));
    // Note: Missing "call" attribute
    classgets(res, mkString("dist"));
  }

  UNPROTECT(nprotect);
  return res;
}
