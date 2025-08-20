#include "collapse_c.h" // Needs to be first because includes OpenMP, to avoid namespace conflicts.
#include "data.table.h"
#include "kit.h"


SEXP match_single(SEXP x, SEXP table, SEXP nomatch) {

    // Todo: optimizations for length 1 x or table???
  const int n = length(x), nt = length(table), nmv = asInteger(nomatch);
  if(n == 0) return allocVector(INTSXP, 0);
  if(nt == 0) {
    SEXP nmvint = PROTECT(ScalarInteger(nmv));
    SEXP nint = PROTECT(ScalarInteger(n));
    SEXP sint1 = PROTECT(ScalarInteger(1));
    SEXP res = falloc(nmvint, nint, sint1);
    UNPROTECT(3);
    return res;
  }

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
  int tx = TYPEOF(x), tt = TYPEOF(table);
  // factor is between logical and integer
  if(tx == LGLSXP) tx = INTSXP;
  else if(tx == INTSXP && isFactor(x)) tx -= 1;
  else if(tx == REALSXP && isObject(x) && INHERITS(x, char_integer64) && !INHERITS(table, char_integer64)) {
    PROTECT(x = integer64toREAL(x)); ++nprotect;
  }
  if(tt == LGLSXP) tt = INTSXP;
  else if(tt == INTSXP && isFactor(table)) tt -= 1;
  else if(tt == REALSXP && isObject(table) && INHERITS(table, char_integer64) && !INHERITS(x, char_integer64)) {
    PROTECT(table = integer64toREAL(table)); ++nprotect;
  }

  if(tx != tt) {
    if(tx < tt) { // table could be integer, double, complex, character....
      if(tx == INTSXP-1) { // For factors there is a shorthand: just match the levels against table...
        SEXP nmvint = PROTECT(ScalarInteger(nmv)); ++nprotect;
        SEXP tab = PROTECT(match_single(getAttrib(x, R_LevelsSymbol), table, nmvint)); ++nprotect;
        int *pans = INTEGER(ans), *pt = INTEGER(tab), *px = INTEGER(x);
        if(inherits(x, "na.included")) {
          #pragma omp simd
          for(int i = 0; i < n; ++i) pans[i] = pt[px[i]-1];
        } else {
          int na_ind = 0;
          // Need to take care of possible NA matches in table..
          switch(tt) {
            case INTSXP: {
              const int *ptt = INTEGER_RO(table);
              for(int i = 0; i != nt; ++i) {
                if(ptt[i] == NA_INTEGER) {
                  na_ind = i+1; break;
                }
              }
            } break;
            case REALSXP: {
              const double *ptt = REAL_RO(table);
              for(int i = 0; i != nt; ++i) {
                if(ISNAN(ptt[i])) {
                  na_ind = i+1; break;
                }
              }
            } break;
            case STRSXP: {
              const SEXP *ptt = STRING_PTR_RO(table);
              for(int i = 0; i != nt; ++i) {
                if(ptt[i] == NA_STRING) {
                  na_ind = i+1; break;
                }
              }
            } break;
            case CPLXSXP: {
              const Rcomplex *ptt = COMPLEX_RO(table);
              for(int i = 0; i != nt; ++i) {
                if(C_IsNA(ptt[i]) || C_IsNaN(ptt[i])) {
                  na_ind = i+1; break;
                }
              }
            } break;
            default: error("Type %s for 'table' is not supported.", type2char(tt));
          }
          if(na_ind == 0) na_ind = nmv;
          #pragma omp simd
          for(int i = 0; i < n; ++i) pans[i] = px[i] == NA_INTEGER ? na_ind : pt[px[i]-1];
        }
        UNPROTECT(nprotect);
        return ans;
      }
      PROTECT(x	= coerceVector(x,	tt)); ++nprotect; // Coercing to largest common type
    } else { // x has a larger type than table...
      if(tt == INTSXP-1) { // There could be a complicated shorthand involving matching x against the levels and then replacing this by the first occurrence index
        PROTECT(table = asCharacterFactor(table)); ++nprotect;
        if(tx != STRSXP) { // Worst case: need to coerce x as well to make the match
          PROTECT(x = coerceVector(x, STRSXP)); ++nprotect;
        }
      } else {
        PROTECT(table = coerceVector(table,	tx)); ++nprotect;
      }
    }
  } else if(tx == INTSXP-1 && tt == INTSXP-1) { // Both factors
    SEXP x_lev = PROTECT(getAttrib(x, R_LevelsSymbol)); ++nprotect; // Unnecessary but appeases RCHK
    if(!R_compute_identical(x_lev, getAttrib(table, R_LevelsSymbol), 0)) {
      // This is the inefficient way: coercing both to character
      // PROTECT(x = asCharacterFactor(x)); ++nprotect;
      // PROTECT(table = asCharacterFactor(table)); ++nprotect;

      // The efficient solution: matching the levels and regenerating table, taking zero as nomatch value here so that NA does not get matched against NA in x
      SEXP sint0 = PROTECT(ScalarInteger(0)); ++nprotect;
      SEXP tab_ilev = PROTECT(match_single(getAttrib(table, R_LevelsSymbol), x_lev, sint0)); ++nprotect;
      SEXP table_new = PROTECT(duplicate(table)); ++nprotect;
      subsetVectorRaw(table_new, tab_ilev, table, /*anyNA=*/!inherits(table, "na.included"));
      table = table_new;
    }
  }

  tx = TYPEOF(x);
  int K = 0, anyNA = 0;
  size_t M;
  // if(n >= INT_MAX) error("Length of 'x' is too large. (Long vector not supported yet)"); // 1073741824
  if (tx == STRSXP || tx == REALSXP || tx == CPLXSXP || (tx == INTSXP && !isObject(x))) {
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
      SEXP ngtab = getAttrib(table, sym_n_groups);
      if(isNull(ngtab)) goto bigint;
      int ng = asInteger(getAttrib(x, sym_n_groups)), ngt = asInteger(ngtab);
      if(ngt > ng) ng = ngt;
      M = (size_t)ng + 2;
      tx = 1000;
    } else goto bigint;
    anyNA = !(inherits(x, "na.included") && inherits(table, "na.included"));
  } else if (tx == LGLSXP) {
    M = 3;
  } else error("Type %s is not supported.", type2char(tx));

  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values, table has size M
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
    if (need2utf8(x)) {
      PROTECT(x = coerceUtf8IfNeeded(x)); ++nprotect;
    }
    if (need2utf8(table)) {
      PROTECT(table = coerceUtf8IfNeeded(table)); ++nprotect;
    }
    const SEXP *restrict px = SEXPPTR_RO(x), *restrict pt = SEXPPTR_RO(table);
    // fill hash table with indices of 'table'
    for (int i = 0; i != nt; ++i) {
      id = HASH(((uintptr_t) pt[i] & 0xffffffff), K);
      while(h[id]) {
        if(pt[h[id]-1] == pt[i]) goto sbl;
        if(++id >= M) id = 0;
      }
      h[id] = i + 1;
      sbl:;
    }
    // look up values of x in hash table
    for (int i = 0; i != n; ++i) {
      id = HASH(((uintptr_t) px[i] & 0xffffffff), K);
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
  R_Free(h);
  UNPROTECT(nprotect);
  return ans;
}


// Outsourcing the conversions to a central function

SEXP coerce_single_to_equal_types(SEXP x, SEXP table) {

  int nprotect = 1;
  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, x);
  SET_VECTOR_ELT(out, 1, table);

  // https://github.com/wch/r-source/blob/433b0c829018c7ad8cd6a585bf9c388f8aaae303/src/main/unique.c#L1356C4-L1356C4
  if(TYPEOF(x) == CPLXSXP || TYPEOF(x) > STRSXP) SET_VECTOR_ELT(out, 0, coerceVector(x, STRSXP));
  if(TYPEOF(table) == CPLXSXP || TYPEOF(table) > STRSXP) SET_VECTOR_ELT(out, 1, coerceVector(table, STRSXP));
  x = VECTOR_ELT(out, 0);
  table = VECTOR_ELT(out, 1);
  int tx = TYPEOF(x), tt = TYPEOF(table);
  if(tx == LGLSXP) tx = INTSXP;
  else if(tx == INTSXP && isFactor(x)) tx -= 1;
  else if(tx == REALSXP && isObject(x) && INHERITS(x, char_integer64) && !INHERITS(table, char_integer64)) {
    SET_VECTOR_ELT(out, 0, integer64toREAL(x)); x = VECTOR_ELT(out, 0);
  }
  if(tt == LGLSXP) tt = INTSXP;
  else if(tt == INTSXP && isFactor(table)) tt -= 1;
  else if(tt == REALSXP && isObject(table) && INHERITS(table, char_integer64) && !INHERITS(x, char_integer64)) {
    SET_VECTOR_ELT(out, 1, integer64toREAL(table)); table = VECTOR_ELT(out, 1);
  }

  if(tx != tt) {
    if(tx > tt) {
      if(tt == INTSXP-1) { // TODO: could implement as in single case..
        SET_VECTOR_ELT(out, 1, asCharacterFactor(table));
        if(tx != STRSXP) SET_VECTOR_ELT(out, 0, coerceVector(x, STRSXP));
      } else SET_VECTOR_ELT(out, 1, coerceVector(table, tx));
    } else {
      if(tx == INTSXP-1) { // TODO: could implement as in single case..
        SET_VECTOR_ELT(out, 0, asCharacterFactor(x));
        if(tt != STRSXP) SET_VECTOR_ELT(out, 1, coerceVector(table, STRSXP));
      } else SET_VECTOR_ELT(out, 0, coerceVector(x, tt));
    }
  } else if(tx == INTSXP-1 && tt == INTSXP-1) { // Both factors
    SEXP x_lev = PROTECT(getAttrib(x, R_LevelsSymbol)); ++nprotect; // Unnecessary but appeases RCHK
    if(!R_compute_identical(x_lev, getAttrib(table, R_LevelsSymbol), 0)) {
      SEXP sint0 = PROTECT(ScalarInteger(0)); ++nprotect;
      SEXP tab_ilev = PROTECT(match_single(getAttrib(table, R_LevelsSymbol), x_lev, sint0)); ++nprotect;
      SEXP table_new;
      SET_VECTOR_ELT(out, 1, table_new = duplicate(table));
      subsetVectorRaw(table_new, tab_ilev, table, /*anyNA=*/!inherits(table, "na.included")); // TODO: check this !!
    }
  }

  UNPROTECT(nprotect);
  return out;
}


SEXP coerce_to_equal_types(SEXP x, SEXP table) {

  if(TYPEOF(x) == VECSXP || TYPEOF(table) == VECSXP) {
    if(TYPEOF(x) != TYPEOF(table)) error("x and table must both be lists when one is a list");
    int l = length(x);
    if(length(table) != l) error("lengths of x and table must be equal of both are lists");
    SEXP out = PROTECT(allocVector(VECSXP, l));
    for(int i = 0; i < l; i++) {
      SEXP xi = VECTOR_ELT(x, i);
      SEXP ti = VECTOR_ELT(table, i);
      SET_VECTOR_ELT(out, i, coerce_single_to_equal_types(xi, ti));
    }
    UNPROTECT(1);
    return out;
  }

  return coerce_single_to_equal_types(x, table);
}


// Still See: https://www.cockroachlabs.com/blog/vectorized-hash-joiner/

SEXP match_two_vectors(SEXP x, SEXP table, SEXP nomatch) {

  if(TYPEOF(x) != VECSXP || TYPEOF(table) != VECSXP) error("both x and table need to be atomic vectors or lists");
  const int l = length(x), lt = length(table), nmv = asInteger(nomatch);
  if(l == 0) return allocVector(INTSXP, 0);
  if(lt == 0) {
    SEXP nmvint = PROTECT(ScalarInteger(nmv));
    SEXP lx0 = PROTECT(ScalarInteger(length(VECTOR_ELT(x, 0))));
    SEXP sint1 = PROTECT(ScalarInteger(1));
    SEXP res = falloc(nmvint, lx0, sint1);
    UNPROTECT(3);
    return res;
  }

  if(l != lt) error("length(n) must match length(nt)");
  if(l != 2) error("Internal function match_two_vectors() only supports lists of length 2");

  // Shallow copy and coercing as necessary
  int nprotect = 1;
  SEXP clist = PROTECT(coerce_to_equal_types(x, table));
  const SEXP *pc = SEXPPTR_RO(clist), *pc1 = SEXPPTR_RO(pc[0]), *pc2 = SEXPPTR_RO(pc[1]);
  const int n = length(pc1[0]), nt = length(pc1[1]);
  if(n != length(pc2[0])) error("both vectors in x must have the same length");
  if(nt != length(pc2[1])) error("both vectors in table must have the same length");

  int K = 0;
  size_t M;
  const size_t n2 = 2U * (size_t) nt;
  M = 256;
  K = 8;
  while (M < n2) {
    M *= 2;
    K++;
  }

  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values, table has size M
  SEXP ans = PROTECT(allocVector(INTSXP, n)); ++nprotect;
  int *restrict pans = INTEGER(ans);
  size_t id = 0;

  int t1 = TYPEOF(pc1[0]), t2 = TYPEOF(pc2[0]);
  if(t1 == LGLSXP) t1 = INTSXP;
  if(t2 == LGLSXP) t2 = INTSXP;

  // 6 cases: 3 same type and 3 different types
  if(t1 == t2) { // same type
    switch(t1) {
      case INTSXP: {
        const int *restrict px1 = INTEGER(pc1[0]), *restrict px2 = INTEGER(pc2[0]),
                  *restrict pt1 = INTEGER(pc1[1]), *restrict pt2 = INTEGER(pc2[1]);
        // fill hash table with indices of 'table'
        for (int i = 0; i != nt; ++i) {
          id = HASH(pt1[i] + (64988430769U * pt2[i]), K);
          while(h[id]) {
            if(pt1[h[id]-1] == pt1[i] && pt2[h[id]-1] == pt2[i]) goto ibl;
            if(++id >= M) id = 0;
          }
          h[id] = i + 1;
          ibl:;
        }
        // look up values of x in hash table
        for (int i = 0; i != n; ++i) {
          id = HASH(px1[i] + (64988430769U * px2[i]), K);
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
        for(int i = 0; i < 2; ++i) {
          if(need2utf8(pc1[i])) SET_VECTOR_ELT(pc[0], i, coerceUtf8IfNeeded(pc1[i]));
          if(need2utf8(pc2[i])) SET_VECTOR_ELT(pc[1], i, coerceUtf8IfNeeded(pc2[i]));
        }
        const SEXP *restrict px1 = SEXPPTR_RO(pc1[0]), *restrict px2 = SEXPPTR_RO(pc2[0]),
                   *restrict pt1 = SEXPPTR_RO(pc1[1]), *restrict pt2 = SEXPPTR_RO(pc2[1]);
        // fill hash table with indices of 'table'
        for (int i = 0; i != nt; ++i) {
          id = HASH(64988430769U * ((uintptr_t)pt1[i] & 0xffffffff) + ((uintptr_t)pt2[i] & 0xffffffff), K);
          while(h[id]) {
            if(pt1[h[id]-1] == pt1[i] && pt2[h[id]-1] == pt2[i]) goto sbl;
            if(++id >= M) id = 0;
          }
          h[id] = i + 1;
          sbl:;
        }
        // look up values of x in hash table
        for (int i = 0; i != n; ++i) {
          id = HASH(64988430769U * ((uintptr_t)px1[i] & 0xffffffff) + ((uintptr_t)px2[i] & 0xffffffff), K);
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
        const double *restrict px1 = REAL(pc1[0]), *restrict px2 = REAL(pc2[0]),
                     *restrict pt1 = REAL(pc1[1]), *restrict pt2 = REAL(pc2[1]);
        union uno tpv1, tpv2;
        // fill hash table with indices of 'table'
        for (int i = 0; i != nt; ++i) {
          tpv1.d = pt1[i]; tpv2.d = pt2[i];
          id = HASH((64988430769U * (tpv1.u[0] + tpv1.u[1])) + tpv2.u[0] + tpv2.u[1], K);
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
          id = HASH((64988430769U * (tpv1.u[0] + tpv1.u[1])) + tpv2.u[0] + tpv2.u[1], K);
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
    if((t1 == INTSXP && t2 == REALSXP) || (t1 == REALSXP && t2 == INTSXP)) {
      const int rev = t1 == REALSXP;
      const int *restrict pxi = INTEGER(VECTOR_ELT(pc[rev], 0)), *restrict pti = INTEGER(VECTOR_ELT(pc[rev], 1));
      const double *restrict pxr = REAL(VECTOR_ELT(pc[1-rev], 0)), *restrict ptr = REAL(VECTOR_ELT(pc[1-rev], 1));
      union uno tpv;
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        tpv.d = ptr[i];
        id = HASH((64988430769U * pti[i]) + tpv.u[0] + tpv.u[1], K); // TODO: improve!
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
        id = HASH((64988430769U * pxi[i]) + tpv.u[0] + tpv.u[1], K); // TODO: improve!
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
      const double *restrict pxr = REAL(VECTOR_ELT(pc[rev], 0)), *restrict ptr = REAL(VECTOR_ELT(pc[rev], 1));
      for(int i = 0; i < 2; ++i) {
        if(need2utf8(VECTOR_ELT(pc[1-rev], i))) SET_VECTOR_ELT(pc[1-rev], i, coerceUtf8IfNeeded(VECTOR_ELT(pc[1-rev], i)));
      }
      const SEXP *restrict pxs = SEXPPTR_RO(VECTOR_ELT(pc[1-rev], 0)), *restrict pts = SEXPPTR_RO(VECTOR_ELT(pc[1-rev], 1));
      union uno tpv;
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        tpv.d = ptr[i];
        id = HASH((tpv.u[0] + tpv.u[1]) * ((uintptr_t)pts[i] & 0xffffffff), K);
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
        id = HASH((tpv.u[0] + tpv.u[1]) * ((uintptr_t)pxs[i] & 0xffffffff), K);
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
    } else if((t1 == INTSXP && t2 == STRSXP) || (t1 == STRSXP && t2 == INTSXP)) {
      const int rev = t1 == STRSXP;
      const int *restrict pxi = INTEGER(VECTOR_ELT(pc[rev], 0)), *restrict pti = INTEGER(VECTOR_ELT(pc[rev], 1));
      for(int i = 0; i < 2; ++i) {
        if(need2utf8(VECTOR_ELT(pc[1-rev], i))) SET_VECTOR_ELT(pc[1-rev], i, coerceUtf8IfNeeded(VECTOR_ELT(pc[1-rev], i)));
      }
      const SEXP *restrict pxs = SEXPPTR_RO(VECTOR_ELT(pc[1-rev], 0)), *restrict pts = SEXPPTR_RO(VECTOR_ELT(pc[1-rev], 1));

      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        id = HASH(pti[i] * ((uintptr_t)pts[i] & 0xffffffff), K); // TODO: improve!
        while(h[id]) {
          if(pts[h[id]-1] == pts[i] && pti[h[id]-1] == pti[i]) goto isbl;
          if(++id >= M) id = 0;
        }
        h[id] = i + 1;
        isbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        id = HASH(pxi[i] * ((uintptr_t)pxs[i] & 0xffffffff), K);
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

  R_Free(h);
  UNPROTECT(nprotect);
  return ans;
}

// TODO: create match_multiple_vectors: a generalization of match_two_vectors that works for multiple vectors
// This will have to involve bucketing and subgroup matching
// Also idea: combine matches using the maximum before the next largest value?


// This is a workhorse function for matching more than 2 vectors: it matches the first two vectors and also
// saves the unique value count and a group-id for the table which is used to match further columns using the same logic
void match_two_vectors_extend(const SEXP *pc, const int nmv, const int n, const int nt, const size_t M, const int K,
                              int *ng, int *pans, int *ptab)  {

  const SEXP *pc1 = SEXPPTR_RO(pc[0]), *pc2 = SEXPPTR_RO(pc[1]);
  if(n != length(pc2[0])) error("both vectors in x must have the same length");
  if(nt != length(pc2[1])) error("both vectors in table must have the same length");

  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values, table has size M
  size_t id = 0;
  int ngt = 0;

  int t1 = TYPEOF(pc1[0]), t2 = TYPEOF(pc2[0]);
  if(t1 == LGLSXP) t1 = INTSXP;
  if(t2 == LGLSXP) t2 = INTSXP;

  // 6 cases: 3 same type and 3 different types
  if(t1 == t2) { // same type
    switch(t1) {
    case INTSXP: {
      const int *restrict px1 = INTEGER(pc1[0]), *restrict px2 = INTEGER(pc2[0]),
                *restrict pt1 = INTEGER(pc1[1]), *restrict pt2 = INTEGER(pc2[1]);
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        id = HASH(pt1[i] + (64988430769U * pt2[i]), K);
        while(h[id]) {
          if(pt1[h[id]-1] == pt1[i] && pt2[h[id]-1] == pt2[i]) {
            ptab[i] = ptab[h[id]-1];
            goto ibl;
          }
          if(++id >= M) id = 0;
        }
        ptab[i] = h[id] = i + 1; ++ngt;
        ibl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        id = HASH(px1[i] + (64988430769U * px2[i]), K);
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
      for(int i = 0; i < 2; ++i) {
        if(need2utf8(pc1[i])) SET_VECTOR_ELT(pc[0], i, coerceUtf8IfNeeded(pc1[i]));
        if(need2utf8(pc2[i])) SET_VECTOR_ELT(pc[1], i, coerceUtf8IfNeeded(pc2[i]));
      }
      const SEXP *restrict px1 = SEXPPTR_RO(pc1[0]), *restrict px2 = SEXPPTR_RO(pc2[0]),
                 *restrict pt1 = SEXPPTR_RO(pc1[1]), *restrict pt2 = SEXPPTR_RO(pc2[1]);
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        id = HASH(64988430769U * ((uintptr_t)pt1[i] & 0xffffffff) + ((uintptr_t)pt2[i] & 0xffffffff), K);
        while(h[id]) {
          if(pt1[h[id]-1] == pt1[i] && pt2[h[id]-1] == pt2[i]) {
            ptab[i] = ptab[h[id]-1];
            goto sbl;
          }
          if(++id >= M) id = 0;
        }
        ptab[i] = h[id] = i + 1; ++ngt;
        sbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        id = HASH(64988430769U * ((uintptr_t)px1[i] & 0xffffffff) + ((uintptr_t)px2[i] & 0xffffffff), K);
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
      const double *restrict px1 = REAL(pc1[0]), *restrict px2 = REAL(pc2[0]),
                   *restrict pt1 = REAL(pc1[1]), *restrict pt2 = REAL(pc2[1]);
      union uno tpv1, tpv2;
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        tpv1.d = pt1[i]; tpv2.d = pt2[i];
        id = HASH((64988430769U * (tpv1.u[0] + tpv1.u[1])) + tpv2.u[0] + tpv2.u[1], K);
        while(h[id]) {
          if(REQUAL(pt1[h[id]-1], pt1[i]) && REQUAL(pt2[h[id]-1], pt2[i])) {
            ptab[i] = ptab[h[id]-1];
            goto rbl;
          }
          if(++id >= M) id = 0;
        }
        ptab[i] = h[id] = i + 1; ++ngt;
        rbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        tpv1.d = px1[i]; tpv2.d = px2[i];
        id = HASH((64988430769U * (tpv1.u[0] + tpv1.u[1])) + tpv2.u[0] + tpv2.u[1], K);
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
    if((t1 == INTSXP && t2 == REALSXP) || (t1 == REALSXP && t2 == INTSXP)) {
      const int rev = t1 == REALSXP;
      const int *restrict pxi = INTEGER(VECTOR_ELT(pc[rev], 0)), *restrict pti = INTEGER(VECTOR_ELT(pc[rev], 1));
      const double *restrict pxr = REAL(VECTOR_ELT(pc[1-rev], 0)), *restrict ptr = REAL(VECTOR_ELT(pc[1-rev], 1));
      union uno tpv;
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        tpv.d = ptr[i];
        id = HASH((64988430769U * pti[i]) + tpv.u[0] + tpv.u[1], K);
        while(h[id]) {
          if(pti[h[id]-1] == pti[i] && REQUAL(ptr[h[id]-1], ptr[i])) {
            ptab[i] = ptab[h[id]-1];
            goto irbl;
          }
          if(++id >= M) id = 0;
        }
        ptab[i] = h[id] = i + 1; ++ngt;
        irbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        tpv.d = pxr[i];
        id = HASH((64988430769U * pxi[i]) + tpv.u[0] + tpv.u[1], K);
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
      const double *restrict pxr = REAL(VECTOR_ELT(pc[rev], 0)), *restrict ptr = REAL(VECTOR_ELT(pc[rev], 1));
      for(int i = 0; i < 2; ++i) {
        if(need2utf8(VECTOR_ELT(pc[1-rev], i))) SET_VECTOR_ELT(pc[1-rev], i, coerceUtf8IfNeeded(VECTOR_ELT(pc[1-rev], i)));
      }
      const SEXP *restrict pxs = SEXPPTR_RO(VECTOR_ELT(pc[1-rev], 0)), *restrict pts = SEXPPTR_RO(VECTOR_ELT(pc[1-rev], 1));
      union uno tpv;
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        tpv.d = ptr[i];
        id = HASH((tpv.u[0] + tpv.u[1]) * ((uintptr_t)pts[i] & 0xffffffff), K);
        while(h[id]) {
          if(pts[h[id]-1] == pts[i] && REQUAL(ptr[h[id]-1], ptr[i])) { // TODO: which comparison is more expensive?
            ptab[i] = ptab[h[id]-1];
            goto rsbl;
          }
          if(++id >= M) id = 0;
        }
        ptab[i] = h[id] = i + 1; ++ngt;
        rsbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        tpv.d = pxr[i];
        id = HASH((tpv.u[0] + tpv.u[1]) * ((uintptr_t)pxs[i] & 0xffffffff), K);
        while(h[id]) {
          if(pts[h[id]-1] == pxs[i] && REQUAL(ptr[h[id]-1], pxr[i])) { // TODO: which comparison is more expensive?
            pans[i] = h[id];
            goto rsbl2;
          }
          if(++id >= M) id = 0;
        }
        pans[i] = nmv;
        rsbl2:;
      }
      // Third case: integer and string
    } else if((t1 == INTSXP && t2 == STRSXP) || (t1 == STRSXP && t2 == INTSXP)) {
      const int rev = t1 == STRSXP;
      const int *restrict pxi = INTEGER(VECTOR_ELT(pc[rev], 0)), *restrict pti = INTEGER(VECTOR_ELT(pc[rev], 1));
      for(int i = 0; i < 2; ++i) {
        if(need2utf8(VECTOR_ELT(pc[1-rev], i))) SET_VECTOR_ELT(pc[1-rev], i, coerceUtf8IfNeeded(VECTOR_ELT(pc[1-rev], i)));
      }
      const SEXP *restrict pxs = SEXPPTR_RO(VECTOR_ELT(pc[1-rev], 0)), *restrict pts = SEXPPTR_RO(VECTOR_ELT(pc[1-rev], 1));

      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        id = HASH(pti[i] * ((uintptr_t)pts[i] & 0xffffffff), K);
        while(h[id]) {
          if(pti[h[id]-1] == pti[i] && pts[h[id]-1] == pts[i]) {
            ptab[i] = ptab[h[id]-1];
            goto isbl;
          }
          if(++id >= M) id = 0;
        }
        ptab[i] = h[id] = i + 1; ++ngt;
        isbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        id = HASH(pxi[i] * ((uintptr_t)pxs[i] & 0xffffffff), K);
        while(h[id]) {
          if(pti[h[id]-1] == pxi[i] && pts[h[id]-1] == pxs[i]) {
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

  *ng = ngt;
  R_Free(h); // Free hash table
}

// Helper function to match an additional vector
void match_additional(const SEXP *pcj, const int nmv, const int n, const int nt, const size_t M, const int K,
                      int *ng, int *pans_copy, int *pans, int *ptab_copy, int *ptab) {

  if(length(pcj[0]) != n) error("all vectors in x must have the same length");
  if(length(pcj[1]) != nt) error("all vectors in table must have the same length");

  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values, table has size M
  size_t id = 0;

  const unsigned int mult = (M-1) / nt; // TODO: This faster? or better hash ans ? -> Seems faster ! but possible failures ?
  int ngt = 0;

  // Copies really needed??
  memcpy(pans_copy, pans, n * sizeof(int));
  memcpy(ptab_copy, ptab, nt * sizeof(int));

  // TODO: Special case for factors !!!!
  switch(TYPEOF(pcj[0])) {
    case INTSXP:
    case LGLSXP: {
      const int *restrict px = INTEGER(pcj[0]), *restrict pt = INTEGER(pcj[1]);
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        if(ptab_copy[i] == nmv) {
          ++ngt;
          continue;
        }
        id = (ptab_copy[i]*mult) ^ HASH(pt[i], K); // HASH(ptab_copy[i], K)
        while(h[id]) {
          if(ptab_copy[h[id]-1] == ptab_copy[i] && pt[h[id]-1] == pt[i]) {
            ptab[i] = ptab[h[id]-1];
            goto itbl;
          }
          if(++id >= M) id = 0;
        }
        ptab[i] = h[id] = i + 1; ++ngt;
        itbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        if(pans_copy[i] == nmv) continue;
        id = (pans_copy[i]*mult) ^ HASH(px[i], K); // HASH(pans_copy[i], K)
        while(h[id]) {
          if(ptab_copy[h[id]-1] == pans_copy[i] && pt[h[id]-1] == px[i]) {
            pans[i] = h[id];
            goto itbl2;
          }
          if(++id >= M) id = 0;
        }
        pans[i] = nmv;
        itbl2:;
      }
    } break;
    case STRSXP: {
      const SEXP *restrict px = SEXPPTR_RO(PROTECT(coerceUtf8IfNeeded(pcj[0]))),
                 *restrict pt = SEXPPTR_RO(PROTECT(coerceUtf8IfNeeded(pcj[1])));
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        if(ptab_copy[i] == nmv) {
          ++ngt;
          continue;
        }
        id = (ptab_copy[i]*mult) ^ HASH(((uintptr_t) pt[i] & 0xffffffff), K); // HASH(ptab_copy[i], K)
        while(h[id]) {
          if(ptab_copy[h[id]-1] == ptab_copy[i] && pt[h[id]-1] == pt[i]) {
            ptab[i] = ptab[h[id]-1];
            goto stbl;
          }
          if(++id >= M) id = 0;
        }
        ptab[i] = h[id] = i + 1; ++ngt;
        stbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        if(pans_copy[i] == nmv) continue;
        id = (pans_copy[i]*mult) ^ HASH(((uintptr_t) px[i] & 0xffffffff), K); // HASH(pans_copy[i], K)
        while(h[id]) {
          if(ptab_copy[h[id]-1] == pans_copy[i] && pt[h[id]-1] == px[i]) {
            pans[i] = h[id];
            goto stbl2;
          }
          if(++id >= M) id = 0;
        }
        pans[i] = nmv;
        stbl2:;
      }
      UNPROTECT(2);
    } break;
    case REALSXP: {
      const double *restrict px = REAL(pcj[0]), *restrict pt = REAL(pcj[1]);
      union uno tpv;
      // fill hash table with indices of 'table'
      for (int i = 0; i != nt; ++i) {
        if(ptab_copy[i] == nmv) {
          ++ngt;
          continue;
        }
        tpv.d = pt[i];
        id = (ptab_copy[i]*mult) ^ HASH(tpv.u[0] + tpv.u[1], K); // HASH(ptab_copy[i], K)
        while(h[id]) {
          if(ptab_copy[h[id]-1] == ptab_copy[i] && REQUAL(pt[h[id]-1], pt[i])) {
            ptab[i] = ptab[h[id]-1];
            goto rtbl;
          }
          if(++id >= M) id = 0;
        }
        ptab[i] = h[id] = i + 1; ++ngt;
        rtbl:;
      }
      // look up values of x in hash table
      for (int i = 0; i != n; ++i) {
        if(pans_copy[i] == nmv) continue;
        tpv.d = px[i];
        id = (pans_copy[i]*mult) ^ HASH(tpv.u[0] + tpv.u[1], K); // HASH(pans_copy[i], K)
        while(h[id]) {
          if(ptab_copy[h[id]-1] == pans_copy[i] && REQUAL(pt[h[id]-1], px[i])) {
            pans[i] = h[id];
            goto rtbl2;
          }
          if(++id >= M) id = 0;
        }
        pans[i] = nmv;
        rtbl2:;
      }
    } break;
    default: error("Type %s is not supported.", type2char(TYPEOF(pcj[0]))); // Should never be reached
  }

  *ng = ngt;
  R_Free(h); // Free hash table
}

// This is after unique table rows have already been found, we simply need to check if the remaining columns are equal...
void match_rest(const SEXP *pcj, const int nmv, const int n, const int nt, int *pans) {

  if(length(pcj[0]) != n) error("all vectors in x must have the same length");
  if(length(pcj[1]) != nt) error("all vectors in table must have the same length");

  switch(TYPEOF(pcj[0])) {
    case INTSXP:
    case LGLSXP: {
      const int *restrict px = INTEGER(pcj[0]), *restrict pt = INTEGER(pcj[1])-1;
      for (int i = 0; i != n; ++i) {
        if(pans[i] == nmv) continue;
        if(px[i] != pt[pans[i]]) pans[i] = nmv;
      }
    } break;
    case STRSXP: {
      const SEXP *restrict px = SEXPPTR_RO(PROTECT(coerceUtf8IfNeeded(pcj[0]))), *restrict pt = SEXPPTR_RO(PROTECT(coerceUtf8IfNeeded(pcj[1])))-1;
      for (int i = 0; i != n; ++i) {
        if(pans[i] == nmv) continue;
        if(px[i] != pt[pans[i]]) pans[i] = nmv;
      }
      UNPROTECT(2);
    } break;
    case REALSXP: {
      const double *restrict px = REAL(pcj[0]), *restrict pt = REAL(pcj[1])-1;
      for (int i = 0; i != n; ++i) {
        if(pans[i] == nmv) continue;
        if(!REQUAL(px[i], pt[pans[i]])) pans[i] = nmv;
      }
    } break;
    default: error("Type %s is not supported.", type2char(TYPEOF(pcj[0]))); // Should never be reached
  }
}


SEXP match_multiple(SEXP x, SEXP table, SEXP nomatch, SEXP overid)  {

  if(TYPEOF(x) != VECSXP || TYPEOF(table) != VECSXP) error("both x and table need to be atomic vectors or lists");
  const int l = length(x), lt = length(table), nmv = asInteger(nomatch);
  if(l == 0) return allocVector(INTSXP, 0);
  if(lt == 0) {
    SEXP nmvint = PROTECT(ScalarInteger(nmv));
    SEXP lx0 = PROTECT(ScalarInteger(length(VECTOR_ELT(x, 0))));
    SEXP sint1 = PROTECT(ScalarInteger(1));
    SEXP res = falloc(nmvint, lx0, sint1);
    UNPROTECT(3);
    return res;
  }
  if(l != lt) error("length(n) must match length(nt)");

  // Shallow copy and coercing as necessary
  SEXP clist = PROTECT(coerce_to_equal_types(x, table));
  const SEXP *pc = SEXPPTR_RO(clist);
  const int n = length(VECTOR_ELT(pc[0], 0)), nt = length(VECTOR_ELT(pc[0], 1));

  // Determining size of hash table
  const size_t n2 = 2U * (size_t) nt;
  size_t M = 256;
  int K = 8;
  while (M < n2) {
    M *= 2;
    K++;
  }

  int *restrict ptab = (int*)R_alloc(nt, sizeof(int)); // Table to contain the group-id of table
  int ng = 0; // Number of groups
  SEXP ans = PROTECT(allocVector(INTSXP, n));
  int *restrict pans = INTEGER(ans);

  // Initial matching two vectors
  match_two_vectors_extend(pc, nmv, n, nt, M, K, &ng, pans, ptab);

  // Early termination if table is already unique or we only have 2 vectors (should use match_two_vectors() directly)
  if(l > 2) {
    int oid = asInteger(overid); // 0 = early termination, 1 = proceed with warning, 2 = proceed without warning
    if(oid > 0 || ng != nt) {
      // Need to copy table and ans: enters as first vector
      int *restrict ptab_copy = (int*)R_alloc(nt, sizeof(int));
      int *restrict pans_copy = (int*)R_alloc(n, sizeof(int));
      for (int j = 2; j < l; ++j) {
        if(ng != nt) match_additional(SEXPPTR_RO(pc[j]), nmv, n, nt, M, K, &ng, pans_copy, pans, ptab_copy, ptab);
        else {
          if(oid == 1) warning("Overidentified match/join: the first %d of %d columns uniquely match the records. With overid > 0, fmatch() continues to match columns. Consider removing columns or setting overid = 0 to terminate the algorithm after %d columns (the results may differ, see ?fmatch). Alternatively set overid = 2 to silence this warning.", j, l/oid++, j);
          if(oid <= 0) break;
          match_rest(SEXPPTR_RO(pc[j]), nmv, n, nt, pans);
        }
      }
    }
  }

  UNPROTECT(2);
  return ans;
}


SEXP fmatch_internal(SEXP x, SEXP table, SEXP nomatch, SEXP overid) {
  if(TYPEOF(x) == VECSXP) {
    if(length(x) == 2) return match_two_vectors(x, table, nomatch);
    if(length(x) == 1) return match_single(VECTOR_ELT(x, 0), VECTOR_ELT(table, 0), nomatch);
    return match_multiple(x, table, nomatch, overid);
  }
  return match_single(x, table, nomatch);
}

void count_match(SEXP res, int nt, int nmv) {
  const int *restrict pres = INTEGER(res);
  int n = length(res), nd = 0, nnm = 0;
  int *restrict cnt = (int*)R_Calloc(nt+1, int);
  for (int i = 0; i != n; ++i) {
    if(pres[i] == nmv) ++nnm;
    else if(cnt[pres[i]] == 0) {
      cnt[pres[i]] = 1;
      ++nd;
    }
  }
  R_Free(cnt);
  SEXP sym_nomatch = install("N.nomatch");
  SEXP sym_distinct = install("N.distinct");
  setAttrib(res, sym_nomatch, ScalarInteger(nnm));
  setAttrib(res, sym_n_groups, ScalarInteger(nt));
  setAttrib(res, sym_distinct, ScalarInteger(nd));
  classgets(res, mkString("qG"));
}

// This is for export
SEXP fmatchC(SEXP x, SEXP table, SEXP nomatch, SEXP count, SEXP overid) {
  if(asLogical(count) <= 0) return fmatch_internal(x, table, nomatch, overid);
  SEXP res = PROTECT(fmatch_internal(x, table, nomatch, overid));
  int nt = isNewList(table) ? length(VECTOR_ELT(table, 0)) : length(table);
  count_match(res, nt, asInteger(nomatch));
  UNPROTECT(1);
  return res;
}
