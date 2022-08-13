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

  const int n = length(x), lv = length(val), lr = length(rep), tv = TYPEOF(val),
    tx = TYPEOF(x), tr = TYPEOF(rep), ind1 = asLogical(Rind1), invert = asLogical(Rinvert), set = asLogical(Rset);
  int nprotect = 0;

  if(lv > 1 || ind1) {
    if(tv == LGLSXP) {
      if(lv != n) error("If v is a logical vector, length(v) needs to be equal to length(x)");
      if(lr != 1 && lr != n) error("If v is a logical vector, length(r) needs to be 1 or length(x)");
    } else if(tv == INTSXP) {
      if(invert) error("invert = TRUE is only possible if v is a logical vector");
      if(lv > n) error("length(v) must be <= length(x)");
      if(!(lr == 1 || lr == n || lr == lv)) error("length(r) must be either 1, length(v) or length(x)");
      // Just some heuristic checking as this is a programmers function
      const int v1 = INTEGER_ELT(val, 0), vn = INTEGER_ELT(val, lv-1);
      if(v1 < 1 || v1 > n || vn < 1 || vn > n) error("Detected index (v) outside of range [1, length(x)]");
    } else error("If length(v) > 1, v must be an integer or logical vector");
  } else if(lr != 1 && lr != n) error("If length(v) == 1, length(r) must be 1 or length(x)");

  if(lr > 1 && tr != tx) { // lr == n &&
    if(!((tx == INTSXP && tr == LGLSXP) || (tx == LGLSXP && tr == INTSXP))) {
      PROTECT_INDEX ipx;
      PROTECT_WITH_INDEX(rep = coerceVector(rep, tx), &ipx);
      ++nprotect;
    } // error("typeof(x) needs to match typeof(r)");
  }


  SEXP ans = R_NilValue;
  if(set == 0) {
    PROTECT(ans = duplicate(x)); // Fastest?? // copies attributes ?? -> Yes
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

  if(nprotect) UNPROTECT(nprotect);
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
  int n = length(x);
  SEXP ans = PROTECT(allocVector(INTSXP, n));
  int *pans = INTEGER(ans);
  switch(asInteger(isnum)) {
  case 0:
    for(int i = 0; i != n; ++i) pans[i] = TYPEOF(VECTOR_ELT(x, i)) + 1;
    break;
  case 1: // Numeric variables: do_is with op = 100: https://github.com/wch/r-source/blob/2b0818a47199a0b64b6aa9b9f0e53a1e886e8e95/src/main/coerce.c
    {
    if(inherits(x, "indexed_frame")) {
      for(int i = 0; i != n; ++i) {
        SEXP ci = VECTOR_ELT(x, i);
        int tci = TYPEOF(ci);
        pans[i] = (tci == INTSXP && inherits(ci, "integer")) || (tci == REALSXP && inherits(ci, "numeric")); // length(getAttrib(ci, R_ClassSymbol)) <= 2;
      }
    } else {
      for(int i = 0; i != n; ++i) {
        SEXP ci = VECTOR_ELT(x, i);
        int tci = TYPEOF(ci);
        pans[i] = (tci == INTSXP || tci == REALSXP) && OBJECT(ci) == 0;
      }
    }
    SET_TYPEOF(ans, LGLSXP);
    break;
    }
  case 2:
    for(int i = 0; i != n; ++i) pans[i] = (int)isFactor(VECTOR_ELT(x, i));
    SET_TYPEOF(ans, LGLSXP);
    break;
  default: error("Unsupported vtypes option");
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
    default: error("Unsupported SEXP type!");
  }

  copyMostAttrib(x, out);
  UNPROTECT(1);
  return out;
}
