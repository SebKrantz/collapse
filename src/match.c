#include "collapse_c.h" // Needs to be first because includes OpenMP, to avoid namespace conflicts.
#include "kit.h"


SEXP match_single(SEXP x, SEXP table, SEXP nomatch) {

    // Todo: optimizations for length 1 x or table???
  const int n = length(x), nt = length(table), nmv = asInteger(nomatch);
  if(n == 0) return allocVector(INTSXP, 0);
  if(nt == 0) return falloc(ScalarInteger(nmv), ScalarInteger(n), ScalarInteger(1));
  int nprotect = 1;

  // Allocating here. For factors there is a shorthand
  SEXP ans = PROTECT(allocVector(INTSXP, n));

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
    if(TYPEOF(x) < TYPEOF(table)) { // table could be double, complex, character....
      // TODO: What if x is logical and table is factor??
      if(isFactor(x)) { // For factors there is a shorthand: just match the levels against table...
        PROTECT(table = match_single(getAttrib(x, R_LevelsSymbol), table, ScalarInteger(nmv))); ++nprotect;
        int *pans = INTEGER(ans), *pt = INTEGER(table)-1, *px = INTEGER(x);
        if(inherits(x, "na.included")) {
          #pragma omp simd
          for(int i = 0; i < n; ++i) pans[i] = pt[px[i]];
        } else {
          #pragma omp simd
          for(int i = 0; i < n; ++i) pans[i] = px[i] == NA_INTEGER ? nmv : pt[px[i]];
        }
        UNPROTECT(nprotect);
        return ans;
      }
      PROTECT(x	= coerceVector(x,	TYPEOF(table))); ++nprotect; // Coercing to largest common type
    } else { // x has a larger type than table...
      if(isFactor(table)) { // There could be a complicated shorthand involving matching x against the levels and then replacing this by the first occurence index
        PROTECT(table = asCharacterFactor(table)); ++nprotect;
        if(TYPEOF(x) != STRSXP) { // Worst case: need to coerce x as well to make the match
          PROTECT(x = coerceVector(x, STRSXP)); ++nprotect;
        }
      } else {
        PROTECT(table = coerceVector(table,	TYPEOF(x))); ++nprotect;
      }
    }
  } else if(isFactor(x) && isFactor(table)) {
    if(!R_compute_identical(getAttrib(x, R_LevelsSymbol), getAttrib(table, R_LevelsSymbol), 0)) {
      // This is the inefficient way: coercing both to character
      // PROTECT(x = asCharacterFactor(x)); ++nprotect;
      // PROTECT(table = asCharacterFactor(table)); ++nprotect;

      // The efficient solution: matching the levels and regenerating table, taking zero as nomatch value here so that NA does not get matched against NA in x
      SEXP tab_ilev = PROTECT(match_single(getAttrib(table, R_LevelsSymbol), getAttrib(x, R_LevelsSymbol), ScalarInteger(0))); ++nprotect;
      SEXP table_new = PROTECT(duplicate(table)); ++nprotect;
      subsetVectorRaw(table_new, tab_ilev, table, /*anyNA=*/!inherits(table, "na.included"));
      table = table_new;
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



SEXP match_two_vectors(SEXP x, SEXP table, SEXP nomatch) {

  if(TYPEOF(x) != VECSXP || TYPEOF(table) != VECSXP) error("both x and table need to be atomic vectors or lists");
  const int l = length(x), lt = length(table), nmv = asInteger(nomatch);
  if(l == 0) return allocVector(INTSXP, 0);
  if(lt == 0) return falloc(ScalarInteger(nmv), ScalarInteger(length(VECTOR_ELT(x, 0))), ScalarInteger(1));

  if(l != lt) error("length(n) must match length(nt)");
  if(l != 2) error("Internal function match_two_vectors() only supports lists of length 2");

  // Shallow copy in case we need to coerce
  int nprotect = 1;
  SEXP clist = PROTECT(allocVector(VECSXP, 4)), *pc = SEXPPTR(clist);
  for(int j = 0; j < 2; ++j) {
    pc[j] = VECTOR_ELT(x, j);
    pc[j+2] = VECTOR_ELT(table, j);
  }
  const int n = length(pc[0]), nt = length(pc[2]);
  if(n != length(pc[1])) error("both vectors in x must have the same length");
  if(nt != length(pc[3])) error("both vectors in table must have the same length");

  // Now coercing as necessary
  for(int j = 0; j < 2; ++j) {
    if(TYPEOF(pc[j]) == CPLXSXP || TYPEOF(pc[j]) > STRSXP) SET_VECTOR_ELT(clist, j, coerceVector(pc[j], STRSXP));
    if(TYPEOF(pc[j+2]) == CPLXSXP || TYPEOF(pc[j+2]) > STRSXP) SET_VECTOR_ELT(clist, j+2, coerceVector(pc[j+2], STRSXP));
    x = pc[j]; table = pc[j+2];
    if(TYPEOF(x) != TYPEOF(table)) {
      if(TYPEOF(x) > TYPEOF(table)) {
        SEXP tmp = table; table = x; x = tmp;
      }
      if(isFactor(x)) { // TODO: could implement as in single case.. What if x is logical and table is factor??
        SET_VECTOR_ELT(clist, j, asCharacterFactor(x));
        if(TYPEOF(table) != STRSXP) SET_VECTOR_ELT(clist, j+2, coerceVector(table, STRSXP));
      } else SET_VECTOR_ELT(clist, j, coerceVector(x, TYPEOF(table)));
    } else if(isFactor(x) && isFactor(table)) {
      if(!R_compute_identical(getAttrib(x, R_LevelsSymbol), getAttrib(table, R_LevelsSymbol), 0)) {
        SEXP tab_ilev = PROTECT(match_single(getAttrib(table, R_LevelsSymbol), getAttrib(x, R_LevelsSymbol), ScalarInteger(0))); ++nprotect;
        SEXP table_new;
        SET_VECTOR_ELT(clist, j+2, table_new = duplicate(table));
        subsetVectorRaw(table_new, tab_ilev, table, /*anyNA=*/!inherits(table, "na.included")); // TODO: check this !!
      }
    }
  }

  int K = 0;
  size_t M;
  const size_t n2 = 2U * (size_t) nt;
  M = 256;
  K = 8;
  while (M < n2) {
    M *= 2;
    K++;
  }

  int *restrict h = (int*)Calloc(M, int); // Table to save the hash values, table has size M
  SEXP ans = PROTECT(allocVector(INTSXP, n)); ++nprotect;
  int *restrict pans = INTEGER(ans);
  size_t id = 0;

  const int t1 = TYPEOF(pc[0]), t2 = TYPEOF(pc[1]);

  // 6 cases: 3 same type and 3 different types
  if(t1 == t2) { // same type
    switch(t1) {
      case INTSXP:
      case LGLSXP: {
        const int *restrict px1 = INTEGER(pc[0]), *restrict px2 = INTEGER(pc[1]),
                  *restrict pt1 = INTEGER(pc[2]), *restrict pt2 = INTEGER(pc[3]);
        // fill hash table with indices of 'table'
        for (int i = 0; i != nt; ++i) {
          id = HASH((unsigned)pt1[i] * (unsigned)pt2[i], K); // TODO: bitwise? combine hash values?
          while(h[id]) {
            if(pt1[h[id]-1] == pt1[i] && pt2[h[id]-1] == pt2[i]) goto ibl;
            if(++id >= M) id = 0;
          }
          h[id] = i + 1;
          ibl:;
        }
        // look up values of x in hash table
        for (int i = 0; i != n; ++i) {
          id = HASH((unsigned)px1[i] * (unsigned)px2[i], K); // TODO: bitwise? combine hash values?
          while(h[id]) {
            if(pt1[h[id]-1] == px1[i] && pt2[h[id]-1] == px2[i]) {
              pans[i] = h[id];
              goto ibl2;
            }
            if(++id >= M) id = 0;
          }
          pans[i] = nmv;
          ibl2:;
        }
      } break;
      case STRSXP: {
        const SEXP *restrict px1 = STRING_PTR(pc[0]), *restrict px2 = STRING_PTR(pc[1]),
                   *restrict pt1 = STRING_PTR(pc[2]), *restrict pt2 = STRING_PTR(pc[3]);
        // fill hash table with indices of 'table'
        for (int i = 0; i != nt; ++i) {
          id = HASH(((intptr_t) pt1[i] ^ (intptr_t) pt2[i]) & 0xffffffff, K);
          while(h[id]) {
            if(pt1[h[id]-1] == pt1[i] && pt2[h[id]-1] == pt2[i]) goto sbl;
            if(++id >= M) id = 0;
          }
          h[id] = i + 1;
          sbl:;
        }
        // look up values of x in hash table
        for (int i = 0; i != n; ++i) {
          id = HASH(((intptr_t) px1[i] ^ (intptr_t) px2[i]) & 0xffffffff, K);
          while(h[id]) {
            if(pt1[h[id]-1] == px1[i] && pt2[h[id]-1] == px2[i]) {
              pans[i] = h[id];
              goto sbl2;
            }
            if(++id >= M) id = 0;
          }
          pans[i] = nmv;
          sbl2:;
        }
      } break;
      case REALSXP: {
        const double *restrict px1 = REAL(pc[0]), *restrict px2 = REAL(pc[1]),
                     *restrict pt1 = REAL(pc[2]), *restrict pt2 = REAL(pc[3]);
        union uno tpv1, tpv2;
        // fill hash table with indices of 'table'
        for (int i = 0; i != nt; ++i) {
          tpv1.d = pt1[i]; tpv2.d = pt2[i];
          id = HASH(tpv1.u[0] + tpv1.u[1] + tpv2.u[0] + tpv2.u[1], K); // adding all a good idea??
          while(h[id]) {
            if(REQUAL(pt1[h[id]-1], pt1[i]) && REQUAL(pt2[h[id]-1], pt2[i])) goto rbl;
            if(++id >= M) id = 0;
          }
          h[id] = i + 1;
          rbl:;
        }
        // look up values of x in hash table
        for (int i = 0; i != n; ++i) {
          tpv1.d = px1[i]; tpv2.d = px2[i];
          id = HASH(tpv1.u[0] + tpv1.u[1] + tpv2.u[0] + tpv2.u[1], K); // adding all a good idea??
          while(h[id]) {
            if(REQUAL(pt1[h[id]-1], px1[i]) && REQUAL(pt2[h[id]-1], px2[i])) {
              pans[i] = h[id];
              goto rbl2;
            }
            if(++id >= M) id = 0;
          }
          pans[i] = nmv;
          rbl2:;
        }
      } break;
      default: error("Type %s is not supported.", type2char(t1)); // Should never be reached
    }
  } else { // different types
    // First case: integer and real
    if(((t1 == INTSXP || t1 == LGLSXP) && t2 == REALSXP) || (t1 == REALSXP && (t2 == INTSXP || t2 == LGLSXP))) {
      const int rev = t1 == REALSXP;
      const int *restrict pxi = INTEGER(pc[0+rev]), *restrict pti = INTEGER(pc[2+rev]);
      const double *restrict pxr = REAL(pc[1-rev]), *restrict ptr = REAL(pc[3-rev]);
      union uno tpv;
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        tpv.d = ptr[i];
        id = HASH(pti[i] + tpv.u[0] + tpv.u[1], K); // TODO: bitwise? combine hash values?
        while(h[id]) {
          if(pti[h[id]-1] == pti[i] && REQUAL(ptr[h[id]-1], ptr[i])) goto irbl;
          if(++id >= M) id = 0;
        }
        h[id] = i + 1;
        irbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        tpv.d = pxr[i];
        id = HASH(pxi[i] + tpv.u[0] + tpv.u[1], K); // TODO: bitwise? combine hash values?
        while(h[id]) {
          if(pti[h[id]-1] == pxi[i] && REQUAL(ptr[h[id]-1], pxr[i])) {
            pans[i] = h[id];
            goto irbl2;
          }
          if(++id >= M) id = 0;
        }
        pans[i] = nmv;
        irbl2:;
      }

    // Second case: real and string
    } else if ((t1 == REALSXP && t2 == STRSXP) || (t1 == STRSXP && t2 == REALSXP)) {
      const int rev = t1 == STRSXP;
      const double *restrict pxr = REAL(pc[0+rev]), *restrict ptr = REAL(pc[2+rev]);
      const SEXP *restrict pxs = STRING_PTR(pc[1-rev]), *restrict pts = STRING_PTR(pc[3-rev]);
      union uno tpv;
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        tpv.d = ptr[i];
        id = HASH((intptr_t)pts[i] ^ tpv.u[0] ^ tpv.u[1], K); // TODO: bitwise best??
        while(h[id]) {
          if(pts[h[id]-1] == pts[i] && REQUAL(ptr[h[id]-1], ptr[i])) goto rsbl;
          if(++id >= M) id = 0;
        }
        h[id] = i + 1;
        rsbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        tpv.d = pxr[i];
        id = HASH((intptr_t)pxs[i] ^ tpv.u[0] ^ tpv.u[1], K); // TODO: bitwise best??
        while(h[id]) {
          if(pts[h[id]-1] == pxs[i] && REQUAL(ptr[h[id]-1], pxr[i])) {
            pans[i] = h[id];
            goto rsbl2;
          }
          if(++id >= M) id = 0;
        }
        pans[i] = nmv;
        rsbl2:;
      }
    // Third case: integer and string
    } else if(((t1 == INTSXP || t1 == LGLSXP) && t2 == STRSXP) || (t1 == STRSXP && (t2 == INTSXP || t2 == LGLSXP))) {
      const int rev = t1 == STRSXP;
      const int *restrict pxi = INTEGER(pc[0+rev]), *restrict pti = INTEGER(pc[2+rev]);
      const SEXP *restrict pxs = STRING_PTR(pc[1-rev]), *restrict pts = STRING_PTR(pc[3-rev]);

      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        id = HASH((intptr_t)pts[i] ^ pti[i], K); // TODO: bitwise best??
        while(h[id]) {
          if(pts[h[id]-1] == pts[i] && pti[h[id]-1] == pti[i]) goto isbl;
          if(++id >= M) id = 0;
        }
        h[id] = i + 1;
        isbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        id = HASH((intptr_t)pxs[i] ^ pxi[i], K); // TODO: bitwise best??
        while(h[id]) {
          if(pts[h[id]-1] == pxs[i] && pti[h[id]-1] == pxi[i]) {
            pans[i] = h[id];
            goto isbl2;
          }
          if(++id >= M) id = 0;
        }
        pans[i] = nmv;
        isbl2:;
      }
    } else error("Unsupported types: %s and %s", type2char(t1), type2char(t2));
  }

  UNPROTECT(nprotect);
  return ans;
}

// TODO: create match_multiple_vectors: a generalization of match_two_vectors that works for multiple vectors
// This will have to involve bucketing and subgroup matching
// Also idea: combine matches using the maximum before the next largest value?

// Function for export
SEXP fmatchC(SEXP x, SEXP table, SEXP nomatch) {
  if(TYPEOF(x) == VECSXP) return match_two_vectors(x, table, nomatch);
  return match_single(x, table, nomatch);
}
