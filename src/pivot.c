#include "collapse_c.h"


// Needed ?? rbindlist() is already pretty fast...
// SEXP pivot_long_replicate_id_columns(SEXP data, SEXP times) {
//
// }

// Helper for pivot_long
void writeValueByIndex(SEXP target, SEXP source, const int from, SEXP index) {

  const int tt = TYPEOF(target), coerce = TYPEOF(source) != tt, li = length(index);
  if(coerce) source = PROTECT(coerceVector(source, tt));
  if(length(source) < li) error("Attempting to write %d elements to a vector of length %d", li, length(source));
  if(TYPEOF(index) != INTSXP) error("Indices must be integers");
  const int *restrict pi = INTEGER(index);

  // TODO: SIMD??
  switch(tt) {
  case INTSXP:
  case LGLSXP: {
    const int *restrict ps = INTEGER_RO(source)-1;
    int *restrict pt = INTEGER(target)+from;
    for(int i = 0; i != li; ++i) pt[i] = ps[pi[i]];
    break;
  }
  case REALSXP: {
    const double *restrict ps = REAL_RO(source)-1;
    double *restrict pt = REAL(target)+from;
    for(int i = 0; i != li; ++i) pt[i] = ps[pi[i]];
    break;
  }
  case CPLXSXP: {
    const Rcomplex *restrict ps = COMPLEX_RO(source)-1;
    Rcomplex *restrict pt = COMPLEX(target)+from;
    for(int i = 0; i != li; ++i) pt[i] = ps[pi[i]];
    break;
  }
  case RAWSXP: {
    const Rbyte *restrict ps = RAW_RO(source)-1;
    Rbyte *restrict pt = RAW(target)+from;
    for(int i = 0; i != li; ++i) pt[i] = ps[pi[i]];
    break;
  }
  case STRSXP:
  case VECSXP:
  case EXPRSXP: {
    const SEXP *restrict ps = SEXPPTR_RO(source)-1;
    SEXP *restrict pt = SEXPPTR(target)+from;
    for(int i = 0; i != li; ++i) pt[i] = ps[pi[i]];
    break;
  }
  default: error("Unsupported SEXP type: '%s'", type2char(tt));
  }
  if(coerce == 0) return;
  UNPROTECT(1);
}


SEXP pivot_long(SEXP data, SEXP ind, SEXP idcol) {
  if(TYPEOF(data) != VECSXP) error("pivot_long: input data is of type '%s', but needs to be a list", type2char(TYPEOF(data)));
  const int l = length(data);
  if(l == 1) return VECTOR_ELT(data, 0);
  if(l == 0) error("pivot_long: input data needs to have 1 or more columns. Current number of columns: 0");

  const SEXP *pd = SEXPPTR_RO(data), *pind = pd;

  if(!isNull(ind)) {
    if(TYPEOF(ind) != VECSXP) error("pivot_long with missing value removal: list of indices of type '%s', but needs to be a list", type2char(TYPEOF(ind)));
    if(length(ind) != l) error("length(data) must match lenth(indlist)");
    pind = SEXPPTR_RO(ind);
  }

  int max_type = 0, distinct_types = 0, len = 0;
  for (int j = 0, tj, tj_first = TYPEOF(pd[0]), oj, oj_first = OBJECT(pd[0]); j != l; ++j) {
    tj = TYPEOF(pd[j]);
    oj = OBJECT(pd[j]);
    len += length(pind[j]);
    if(tj > max_type) max_type = tj;
    if(tj != tj_first || oj != oj_first) distinct_types = 1;
  }

  SEXP res;
  // Case 1: no indices, which means we simply melt a single column: same as rbindlist()
  if(isNull(ind)) {
    res = PROTECT(allocVector(max_type, len)); len = 0;
    for (int j = 0; j != l; ++j) {
      int tmp = length(pd[j]);
      writeValue(res, pd[j], len, tmp); // from data.table_rbindlist.c
      len += tmp;
    }
  } else {
  // Now the more interesting case: we have a list of indices for the non-missing cases of each column.
    res = PROTECT(allocVector(max_type, len)); len = 0;
    for (int j = 0; j != l; ++j) {
      writeValueByIndex(res, pd[j], len, pind[j]); // See above
      len += length(pind[j]);
    }
  }

  if(distinct_types == 0) {
    copyMostAttrib(pd[0], res);
    // setAttrib(res, install("label"), R_NilValue); // better to keep, this is also used for id-columns if na.rm = TRUE
  }

  // Add ID column
  if(asLogical(idcol)) {
    SEXP names = getAttrib(data, R_NamesSymbol);
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SEXP id_column;
    SET_VECTOR_ELT(result, 0, id_column = allocVector(isNull(names) ? INTSXP : STRSXP, length(res)));
    SET_VECTOR_ELT(result, 1, res);
    if(isNull(names)) {
      int *restrict pid = INTEGER(id_column);
      for (int j = 0, end = 0, v = 1; j != l; ++j) {
        end = length(pind[j]); // SIMD??
        for (int i = 0; i != end; ++i) pid[i] = v;
        pid += end; ++v;
      }
    } else {
      SEXP *restrict pid = STRING_PTR(id_column), *pnam = STRING_PTR(names);
      for (int j = 0, end = 0; j != l; ++j) {
        SEXP namj = pnam[j];
        end = length(pind[j]); // SIMD??
        for (int i = 0; i != end; ++i) pid[i] = namj;
        pid += end;
      }
    }
    UNPROTECT(2);
    return result;
  }

  UNPROTECT(1);
  return res;
}

// TODO: How to check for duplicate rows?
SEXP pivot_wide(SEXP index, SEXP id, SEXP column, SEXP fill, SEXP Rnthreads) {

  SEXP sym_ng = install("N.groups");
  const int *restrict pix = INTEGER_RO(index), *restrict pid = INTEGER_RO(id), l = length(index),
    nr = asInteger(getAttrib(index, sym_ng)),
    nc = asInteger(getAttrib(id, sym_ng)), tx = TYPEOF(column);
  if(l != length(id)) error("Internal error: length(index) must match length(id)");
  if(l != length(column)) error("Internal error: length(index) must match length(column)");
  if(nr < 1 || nc < 1) error("Resulting data frame after pivoting needs to have at least one row and column");

  int nthreads = asInteger(Rnthreads);
  if(l < 100000) nthreads = 1; // No improvements from multithreading on small data.
  if(nthreads > max_threads) nthreads = max_threads;

  SEXP out = PROTECT(allocVector(VECSXP, nc)), *restrict pout = SEXPPTR(out)-1;

  SEXP fill_val;
  if(fill == R_NilValue) {
    fill_val = tx == REALSXP ? ScalarReal(NA_REAL) : tx == INTSXP ? ScalarInteger(NA_INTEGER) :
    tx == LGLSXP ? ScalarLogical(NA_LOGICAL) : tx == STRSXP ? ScalarString(NA_STRING) :
    tx == CPLXSXP ? ScalarComplex(asComplex(ScalarReal(NA_REAL))) : tx == RAWSXP ? ScalarRaw(0) : R_NilValue;
  } else if(TYPEOF(fill) == tx) {
    fill_val = fill;
  } else fill_val = coerceVector(fill, tx);
  PROTECT(fill_val);
  SEXP out1;
  SET_VECTOR_ELT(out, 0, out1 = falloc(fill_val, ScalarInteger(nr), ScalarLogical(1)));
  copyMostAttrib(column, out1); // TODO: Check that this works!!
  // TODO: can multithread?? -> NOPE!, as expected
  for (int j = 1; j < nc; ++j) SET_VECTOR_ELT(out, j, duplicate(out1));


  // TODO: SIMD: doesn't vectorize on clang 16. Also multithreading gives only minor performance improvements..
  switch(tx) {
    case INTSXP:
    case LGLSXP: {
      const int *restrict pc = INTEGER_RO(column);
      if(nthreads <= 1) {
        for(int i = 0; i != l; ++i) INTEGER(pout[pid[i]])[pix[i]-1] = pc[i];
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int i = 0; i < l; ++i) INTEGER(pout[pid[i]])[pix[i]-1] = pc[i];
      }
      break;
    }
    case REALSXP: {
      const double *restrict pc = REAL_RO(column);
      if(nthreads <= 1) {
        for(int i = 0; i < l; ++i) REAL(pout[pid[i]])[pix[i]-1] = pc[i];
        // // cool idea but not really faster...
        // double *restrict pout_i = REAL(pout[pid[0]])-1;
        // for(int i = 0, prev = pid[0]; i != l; ++i) {
        //   if(pid[i] != prev) pout_i = REAL(pout[pid[i]])-1;
        //   pout_i[pix[i]] = pc[i];
        // }
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int i = 0; i < l; ++i) REAL(pout[pid[i]])[pix[i]-1] = pc[i];
      }
      break;
    }
    case CPLXSXP: {
      const Rcomplex *restrict pc = COMPLEX_RO(column);
      if(nthreads <= 1) {
        for(int i = 0; i != l; ++i) COMPLEX(pout[pid[i]])[pix[i]-1] = pc[i];
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int i = 0; i < l; ++i) COMPLEX(pout[pid[i]])[pix[i]-1] = pc[i];
      }
      break;
    }
    case RAWSXP: {
      const Rbyte *restrict pc = RAW_RO(column);
      if(nthreads <= 1) {
        for(int i = 0; i != l; ++i) RAW(pout[pid[i]])[pix[i]-1] = pc[i];
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int i = 0; i < l; ++i) RAW(pout[pid[i]])[pix[i]-1] = pc[i];
      }
      break;
    }
    case STRSXP:
    case VECSXP:
    case EXPRSXP: {
      const SEXP *restrict pc = SEXPPTR_RO(column);
      if(nthreads <= 1 || tx != STRSXP) {
        for(int i = 0; i != l; ++i) SEXPPTR(pout[pid[i]])[pix[i]-1] = pc[i];
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int i = 0; i < l; ++i) SEXPPTR(pout[pid[i]])[pix[i]-1] = pc[i];
      }
      break;
    }
    default: error("Unsupported SEXP type: '%s'", type2char(tx));
  }
  UNPROTECT(2);
  return out;
}
