#include "collapse_c.h" // Needs to be first because includes OpenMP, to avoid namespace conflicts.
#include "kit.h"


SEXP match_single(SEXP x, SEXP table, SEXP nomatch) {
  // Todo: optimizations for length 1 x or table???
  const int n = length(x), nt = length(table), nmv = asInteger(nomatch);
  if(n == 0) return allocVector(INTSXP, 0);
  if(nt == 0) return falloc(ScalarInteger(nmv), ScalarInteger(n), ScalarInteger(1));
  int nprotect = 1;
  // https://github.com/wch/r-source/blob/433b0c829018c7ad8cd6a585bf9c388f8aaae303/src/main/unique.c#L1356C4-L1356C4
  if(TYPEOF(x) > STRSXP || TYPEOF(table) > STRSXP) {
    if(TYPEOF(x) > STRSXP) {
      PROTECT(x = coerceVector(x, STRSXP)); ++nprotect;
    }
    if(TYPEOF(table) > STRSXP) {
      PROTECT(table = coerceVector(table, STRSXP)); ++nprotect;
    }
  }
  if(TYPEOF(x) != TYPEOF(table)) {
    if(TYPEOF(x) < TYPEOF(table)) {
      PROTECT(x	= coerceVector(x,	TYPEOF(table))); nprotect++;
    } else {
      PROTECT(table	= coerceVector(table,	TYPEOF(x))); nprotect++;
    }
  }
  if(isFactor(x)) {
    if(!R_compute_identical(getAttrib(x, R_LevelsSymbol), getAttrib(table, R_LevelsSymbol), 0)) {
      // TODO: does take care of levels ???
      PROTECT(x = coerceVector(x, STRSXP)); ++nprotect;
      PROTECT(table = coerceVector(table, STRSXP)); ++nprotect;
    }
  }

  int K = 0, tx = TYPEOF(x), anyNA = 0;
  size_t M;
  // if(n >= INT_MAX) error("Length of 'x' is too large. (Long vector not supported yet)"); // 1073741824
  if (tx == STRSXP || tx == REALSXP || tx == CPLXSXP || (tx == INTSXP && OBJECT(x) == 0)) {
    bigint:;
    const size_t n2 = 2U * (size_t) nt;
    M = 256;
    K = 8;
    while (M < n2) {
      M *= 2;
      K++;
    }
  } else if(tx == INTSXP) { // TODO: think about qG objects here...
    if(isFactor(x)) {
      tx = 1000;
      M = (size_t)nlevels(x) + 2;
    } else if(inherits(x, "qG")) {
      SEXP sym_ng = install("N.groups"), ngtab = getAttrib(table, sym_ng);
      if(isNull(ngtab)) goto bigint;
      int ng = asInteger(getAttrib(x, sym_ng)), ngt = asInteger(ngtab);
      if(ngt > ng) ng = ngt;
      M = (size_t)ng + 2;
      tx = 1000;
    } else goto bigint;
    anyNA = !(inherits(x, "na.included") && inherits(table, "na.included"));
  } else if (tx == LGLSXP) {
    M = 3;
  } else error("Type %s is not supported.", type2char(tx));

  int *restrict h = (int*)Calloc(M, int); // Table to save the hash values, table has size M
  SEXP ans = PROTECT(allocVector(INTSXP, n));
  int *restrict pans = INTEGER(ans);
  size_t id = 0;

  switch (tx) {
  case LGLSXP:
  case 1000: // This is for factors or logical vectors where the size of the table is known
  {
    const int *restrict px = INTEGER(x), *restrict pt = INTEGER(table);
    if(tx == 1000 && !anyNA) {
      // fill hash table with indices of 'table'
      for (int i = 0, j; i != nt; ++i) {
        j = pt[i];
        if(h[j]) continue;
        h[j] = i + 1;
      }
      // look up values of x in hash table
      for (int i = 0, j; i != n; ++i) {
        j = px[i];
        pans[i] = h[j] ? h[j] : nmv;
      }
    } else {
      // fill hash table with indices of 'table'
      for (int i = 0, j, k = (int)M-1; i != nt; ++i) {
        j = (pt[i] == NA_INTEGER) ? k : pt[i];
        if(h[j]) continue;
        h[j] = i + 1;
      }
      // look up values of x in hash table
      for (int i = 0, j, k = (int)M-1; i != n; ++i) {
        j = (px[i] == NA_INTEGER) ? k : px[i];
        pans[i] = h[j] ? h[j] : nmv;
      }
    }
  } break;
  case INTSXP: {
    const int *restrict px = INTEGER(x), *restrict pt = INTEGER(table);
    // fill hash table with indices of 'table'
    for (int i = 0; i != nt; ++i) {
      id = HASH(pt[i], K);
      while(h[id]) {
        if(pt[h[id]-1] == pt[i]) goto ibl;
        if(++id >= M) id = 0;
      }
      h[id] = i + 1;
      ibl:;
    }
    // look up values of x in hash table
    for (int i = 0; i != n; ++i) {
      id = HASH(px[i], K);
      while(h[id]) {
        if(pt[h[id]-1] == px[i]) {
          pans[i] = h[id];
          goto ibl2;
        }
        if(++id >= M) id = 0;
      }
      pans[i] = nmv;
      ibl2:;
    }
  } break;
  case REALSXP: {
    const double *restrict px = REAL(x), *restrict pt = REAL(table);
    union uno tpv;
    // fill hash table with indices of 'table'
    for (int i = 0; i != nt; ++i) {
      tpv.d = pt[i];
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        if(REQUAL(pt[h[id]-1], pt[i])) goto rbl;
        if(++id >= M) id = 0;
      }
      h[id] = i + 1;
      rbl:;
    }
    // look up values of x in hash table
    for (int i = 0; i != n; ++i) {
      tpv.d = px[i];
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        if(REQUAL(pt[h[id]-1], px[i])) {
          pans[i] = h[id];
          goto rbl2;
        }
        if(++id >= M) id = 0;
      }
      pans[i] = nmv;
      rbl2:;
    }
  } break;
  case CPLXSXP: {
    const Rcomplex *restrict px = COMPLEX(x), *restrict pt = COMPLEX(table);
    unsigned int u;
    union uno tpv;
    Rcomplex tmp;
    // fill hash table with indices of 'table'
    for (int i = 0; i != nt; ++i) {
      tmp = pt[i];
      if(C_IsNA(tmp)) {
        tmp.r = tmp.i = NA_REAL;
      } else if (C_IsNaN(tmp)) {
        tmp.r = tmp.i = R_NaN;
      }
      tpv.d = tmp.r;
      u = tpv.u[0] ^ tpv.u[1];
      tpv.d = tmp.i;
      u ^= tpv.u[0] ^ tpv.u[1];
      id = HASH(u, K);
      while(h[id]) {
        if(CEQUAL(pt[h[id]-1], pt[i])) goto cbl;
        if(++id >= M) id = 0;
      }
      h[id] = i + 1;
      cbl:;
    }
    // look up values of x in hash table
    for (int i = 0; i != n; ++i) {
      tmp = px[i];
      if(C_IsNA(tmp)) {
        tmp.r = tmp.i = NA_REAL;
      } else if (C_IsNaN(tmp)) {
        tmp.r = tmp.i = R_NaN;
      }
      tpv.d = tmp.r;
      u = tpv.u[0] ^ tpv.u[1];
      tpv.d = tmp.i;
      u ^= tpv.u[0] ^ tpv.u[1];
      id = HASH(u, K);
      while(h[id]) {
        if(CEQUAL(pt[h[id]-1], px[i])) {
          pans[i] = h[id];
          goto cbl2;
        }
        if(++id >= M) id = 0;
      }
      pans[i] = nmv;
      cbl2:;
    }
  } break;
  case STRSXP: {
    const SEXP *restrict px = STRING_PTR(x), *restrict pt = STRING_PTR(table);
    // fill hash table with indices of 'table'
    for (int i = 0; i != nt; ++i) {
      id = HASH(((intptr_t) pt[i] & 0xffffffff), K);
      while(h[id]) {
        if(pt[h[id]-1] == pt[i]) goto sbl;
        if(++id >= M) id = 0;
      }
      h[id] = i + 1;
      sbl:;
    }
    // look up values of x in hash table
    for (int i = 0; i != n; ++i) {
      id = HASH(((intptr_t) px[i] & 0xffffffff), K);
      while(h[id]) {
        if(pt[h[id]-1] == px[i]) {
          pans[i] = h[id];
          goto sbl2;
        }
        if(++id >= M) id = 0;
      }
      pans[i] = nmv;
      sbl2:;
    }
  } break;
  }
  Free(h);
  UNPROTECT(nprotect);
  return ans;
}


