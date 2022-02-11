#include "collapse_c.h"
// #include "data.table.h"
// #ifndef USE_RINTERNALS
// #define USE_RINTERNALS
// #endif

void matCopyAttr(SEXP out, SEXP x, SEXP Rdrop, int ng) {
  SEXP dn = getAttrib(x, R_DimNamesSymbol);
  SEXP cn = isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1); // PROTECT ??
  if(ng == 0 && asLogical(Rdrop)) {
    if(length(cn)) setAttrib(out, R_NamesSymbol, cn);
  } else {
    int nprotect = 1;
    SEXP dim = PROTECT(duplicate(getAttrib(x, R_DimSymbol)));
    INTEGER(dim)[0] = ng == 0 ? 1 : ng;
    dimgets(out, dim);
    if(length(cn)) {
      ++nprotect;
      SEXP dn = PROTECT(allocVector(VECSXP, 2));
      SET_VECTOR_ELT(dn, 0, R_NilValue);
      SET_VECTOR_ELT(dn, 1, cn);
      dimnamesgets(out, dn);
    }
    if(!isObject(x)) copyMostAttrib(x, out);
    UNPROTECT(nprotect);
  }
}

void DFcopyAttr(SEXP out, SEXP x, int ng) {
  SHALLOW_DUPLICATE_ATTRIB(out, x);
  if(ng == 0) {
    setAttrib(out, R_RowNamesSymbol, ScalarInteger(1));
  } else {
    SEXP rn = PROTECT(allocVector(INTSXP, 2)); // Needed here, now unsafe to pass uninitialized vectors to R_RowNamesSymbol.
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -ng;
    setAttrib(out, R_RowNamesSymbol, rn);
    UNPROTECT(1);
  }
}

SEXP geteptr(SEXP x) {
  return R_ExternalPtrProtected(x);
}

// Faster than rep_len(value, n) and slightly faster than matrix(value, n) (which in turn is faster than rep_len)...
SEXP falloc(SEXP value, SEXP n) {
  int l = asInteger(n), tval = TYPEOF(value);
  if(length(value) > 1) error("Must supply a single value to alloc()");
  SEXP out = PROTECT(allocVector(tval, l));
  switch(tval) {
    case INTSXP:
    case LGLSXP: {
      int val = asInteger(value), *pout = INTEGER(out);
      if(val == 0) memset(pout, 0, l*sizeof(int));
      else for(int i = 0; i != l; ++i) pout[i] = val;
      break;
    }
    case REALSXP: {
      double val = asReal(value), *pout = REAL(out);
      if(val == 0.0) memset(pout, 0.0, l*sizeof(double));
      else for(int i = 0; i != l; ++i) pout[i] = val;
      break;
    }
    case STRSXP: {
      SEXP val = asChar(value), *pout = STRING_PTR(out);
      for(int i = 0; i != l; ++i) pout[i] = val;
      break;
    }
    case VECSXP: {
      SEXP *pout = SEXPPTR(out);
      for(int i = 0; i != l; ++i) pout[i] = value;
      break;
    }
    default: error("Not supportd SEXP Type in alloc()");
  }
  copyMostAttrib(value, out);
  UNPROTECT(1);
  return out;
}


SEXP groups2GRP(SEXP x, SEXP lx, SEXP gs) {
  int l = length(x);
  SEXP out = PROTECT(allocVector(INTSXP, asInteger(lx)));
  int *pout = INTEGER(out)-1, *pgs = INTEGER(gs);
  // SEXP *px = VECTOR_PTR(x); // -> Depreciated interface: https://github.com/hadley/r-internals/blob/ea892fa79bbffe961e78dbe9c90ce4ca3bf2d9bc/vectors.md
  // Matt Dowle Commented:
  // VECTOR_PTR does exist but returns 'not safe to return vector pointer' when USE_RINTERNALS is not defined.
  // VECTOR_DATA and LIST_POINTER exist too but call VECTOR_PTR. All are clearly not intended to be used by packages.
  // The concern is overhead inside VECTOR_ELT() biting when called repetitively in a loop like we do here. That's why
  // we take the R API (INTEGER()[i], REAL()[i], etc) outside loops for the simple types even when not parallel. For this
  // type list case (VECSXP) it might be that some items are ALTREP for example, so we really should use the heavier
  // _ELT accessor (VECTOR_ELT) inside the loop in this case.
  SEXP *px = SEXPPTR(x);

  for(int j = l; j--; ) { // This can go in any direction..
    // SEXP column = VECTOR_ELT(x, j);
    int *pcolumn = INTEGER(px[j]), jp = j+1;
    for(int i = pgs[j]; i--; ) pout[pcolumn[i]] = jp; // This can go in any direction...
  }
  UNPROTECT(1);
  return out;
}

// Faster version of base R's spit based on grouping objects..
// TODO: Support DF's !! And check attribute preservation !!
// -> works for factors, Date and POSIXct, but not for POSIXlt
SEXP gsplit(SEXP x, SEXP gobj, SEXP toint) {
  if(TYPEOF(gobj) != VECSXP || !inherits(gobj, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
  const SEXP g = VECTOR_ELT(gobj, 1), gs = VECTOR_ELT(gobj, 2),
    ord = VECTOR_ELT(gobj, 5), order = VECTOR_ELT(gobj, 6);
  const int ng = length(gs), *pgs = INTEGER(gs), tx = TYPEOF(x), l = length(g);
  if(ng != INTEGER(VECTOR_ELT(gobj, 0))[0]) error("'GRP' object needs to have valid vector of group-sizes");
  SEXP res = PROTECT(allocVector(VECSXP, ng));
  // Output as integer or not
  if(asLogical(toint)) {
    for(int i = 0; i != ng; ++i) SET_VECTOR_ELT(res, i, allocVector(INTSXP, pgs[i]));
  } else { // Allocate split vectors and copy attributes and object bits
    SEXP x1 = PROTECT(allocVector(tx, 1));
    copyMostAttrib(x, x1);
    SEXP ax = ATTRIB(x1);
    if(length(ax) == 1 && TAG(ax) == install("label")) ax = R_NilValue;
    int ox = OBJECT(x);
    // FAZIT: Need to use SET_VECTOR_ELT!! pres[i] = allocVector() doesn't work!!
    if(TYPEOF(ax) != NILSXP && ox != 0) {
      for(int i = 0, s4o = IS_S4_OBJECT(x); i != ng; ++i) {
        SEXP resi;
        SET_VECTOR_ELT(res, i, resi = allocVector(tx, pgs[i]));
        SET_ATTRIB(resi, ax);
        SET_OBJECT(resi, ox);
        if(s4o) SET_S4_OBJECT(resi);
      }
    } else if(TYPEOF(ax) != NILSXP) {
      for(int i = 0; i != ng; ++i) {
        SEXP resi;
        SET_VECTOR_ELT(res, i, resi = allocVector(tx, pgs[i])); // SET_ATTRIB(pres[i] = allocVector(tx, pgs[i]), ax);
        SET_ATTRIB(resi, ax);
      }
    } else if(ox != 0) { // Is this even possible? Object bits but no attributes?
      for(int i = 0, s4o = IS_S4_OBJECT(x); i != ng; ++i) {
        SEXP resi;
        SET_VECTOR_ELT(res, i, resi = allocVector(tx, pgs[i]));
        SET_OBJECT(resi, ox);
        if(s4o) SET_S4_OBJECT(resi);
      }
    } else {
      for(int i = 0; i != ng; ++i) SET_VECTOR_ELT(res, i, allocVector(tx, pgs[i]));
    }
    UNPROTECT(1);
  }

  SEXP *pres = SEXPPTR(res);
  // If grouping is sorted
  if(LOGICAL(ord)[1] == 1) { // This only works if data is already ordered in order of the groups
    int count = 0;
    if(asLogical(toint)) {
      for(int j = 0; j != ng; ++j) {
        int *pgj = INTEGER(pres[j]), gsj = pgs[j];
        for(int i = 0; i != gsj; ++i) pgj[i] = ++count;
      }
    } else {
      if(length(x) != l) error("length(x) must match length(g)");
      switch(tx) {
      case INTSXP:
      case LGLSXP: {
        const int *px = INTEGER(x);
        for(int j = 0; j != ng; ++j) {
          int *pgj = INTEGER(pres[j]), gsj = pgs[j];
          for(int i = 0; i != gsj; ++i) pgj[i] = px[count++];
        }
        break;
      }
      case REALSXP: {
        const double *px = REAL(x);
        for(int j = 0, gsj; j != ng; ++j) {
          double *pgj = REAL(pres[j]);
          gsj = pgs[j];
          for(int i = 0; i != gsj; ++i) pgj[i] = px[count++];
        }
        break;
      }
      case CPLXSXP: {
        const Rcomplex *px = COMPLEX(x);
        for(int j = 0, gsj; j != ng; ++j) {
          Rcomplex *pgj = COMPLEX(pres[j]);
          gsj = pgs[j];
          for(int i = 0; i != gsj; ++i) pgj[i] = px[count++];
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = STRING_PTR(x);
        for(int j = 0, gsj; j != ng; ++j) {
          SEXP *pgj = STRING_PTR(pres[j]);
          gsj = pgs[j];
          for(int i = 0; i != gsj; ++i) pgj[i] = px[count++];
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR(x);
        for(int j = 0, gsj; j != ng; ++j) {
          SEXP *pgj = SEXPPTR(pres[j]);
          gsj = pgs[j];
          for(int i = 0; i != gsj; ++i) pgj[i] = px[count++];
        }
        break;
      }
      case RAWSXP: {
        const Rbyte *px = RAW(x);
        for(int j = 0, gsj; j != ng; ++j) {
          Rbyte *pgj = RAW(pres[j]);
          gsj = pgs[j];
          for(int i = 0; i != gsj; ++i) pgj[i] = px[count++];
        }
        break;
      }
      default: error("Unsupported type '%s' passed to gsplit", type2char(tx));
      }
    }
  } else if(length(order) == l) { // Grouping not sorted but we have the ordering..
    SEXP sym_starts = PROTECT(install("starts"));
    const SEXP starts = getAttrib(order, sym_starts);
    UNPROTECT(1);
    if(length(starts) != ng) goto unsno;
    const int *po = INTEGER(order), *ps = INTEGER(starts);

    if(asLogical(toint)) {
      for(int i = 0; i != ng; ++i) {
        int *pri = INTEGER(pres[i]);
        for(int j = ps[i]-1, end = ps[i]+pgs[i]-1, k = 0; j < end; j++) pri[k++] = po[j];
      }
    } else {
      if(length(x) != l) error("length(x) must match length(g)");
      switch(tx) {
      case INTSXP:
      case LGLSXP: {
        const int *px = INTEGER(x);
        for(int i = 0; i != ng; ++i) {
          int *pri = INTEGER(pres[i]);
          for(int j = ps[i]-1, end = ps[i]+pgs[i]-1, k = 0; j < end; ++j) pri[k++] = px[po[j]-1];
        }
        break;
      }
      case REALSXP: {
        double *px = REAL(x);
        for(int i = 0; i != ng; ++i) {
          double *pri = REAL(pres[i]);
          for(int j = ps[i]-1, end = ps[i]+pgs[i]-1, k = 0; j < end; ++j) pri[k++] = px[po[j]-1];
        }
        break;
      }
      case CPLXSXP: {
        Rcomplex *px = COMPLEX(x);
        for(int i = 0; i != ng; ++i) {
          Rcomplex *pri = COMPLEX(pres[i]);
          for(int j = ps[i]-1, end = ps[i]+pgs[i]-1, k = 0; j < end; ++j) pri[k++] = px[po[j]-1];
        }
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x);
        for(int i = 0; i != ng; ++i) {
          SEXP *pri = STRING_PTR(pres[i]);
          for(int j = ps[i]-1, end = ps[i]+pgs[i]-1, k = 0; j < end; ++j) pri[k++] = px[po[j]-1];
        }
        break;
      }
      case VECSXP: {
        SEXP *px = SEXPPTR(x);
        for(int i = 0; i != ng; ++i) {
          SEXP *pri = SEXPPTR(pres[i]);
          for(int j = ps[i]-1, end = ps[i]+pgs[i]-1, k = 0; j < end; ++j) pri[k++] = px[po[j]-1];
        }
        break;
      }
      case RAWSXP: {
        const Rbyte *px = RAW(x);
        for(int i = 0; i != ng; ++i) {
          Rbyte *pri = RAW(pres[i]);
          for(int j = ps[i]-1, end = ps[i]+pgs[i]-1, k = 0; j < end; ++j) pri[k++] = px[po[j]-1];
        }
        break;
      }
      default: error("Unsupported type '%s' passed to gsplit", type2char(tx));
      }
    }

  } else { // Unsorted, without ordering
    unsno:;
    int *count = (int*)Calloc(ng, int);
    // memset(count, 0, sizeof(int)*(ng+1)); // Needed here ??
    // int *count = (int *) R_alloc(ng+1, sizeof(int));

    const int *pg = INTEGER(g);
    // --pres;
    if(asLogical(toint)) {
      for(int i = 0, gi; i != l; ++i) {
        gi = pg[i]-1;
        INTEGER(pres[gi])[count[gi]++] = i+1;
      }
    } else {
      if(length(x) != l) error("length(x) must match length(g)");
      switch(tx) {
      case INTSXP:
      case LGLSXP: {
        const int *px = INTEGER(x);
        for(int i = 0, gi; i != l; ++i) {
          gi = pg[i]-1;
          INTEGER(pres[gi])[count[gi]++] = px[i];
        }
        break;
      }
      case REALSXP: {
        const double *px = REAL(x);
        for(int i = 0, gi; i != l; ++i) {
          gi = pg[i]-1;
          REAL(pres[gi])[count[gi]++] = px[i];
        }
        break;
      }
      case CPLXSXP: {
        const Rcomplex *px = COMPLEX(x);
        for(int i = 0, gi; i != l; ++i) {
          gi = pg[i]-1;
          COMPLEX(pres[gi])[count[gi]++] = px[i];
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = STRING_PTR(x);
        for(int i = 0, gi; i != l; ++i) {
          gi = pg[i]-1;
          STRING_PTR(pres[gi])[count[gi]++] = px[i];
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR(x);
        for(int i = 0, gi; i != l; ++i) {
          gi = pg[i]-1;
          SEXPPTR(pres[gi])[count[gi]++] = px[i];
        }
        break;
      }
      case RAWSXP: {
        const Rbyte *px = RAW(x);
        for(int i = 0, gi; i != l; ++i) {
          gi = pg[i]-1;
          RAW(pres[gi])[count[gi]++] = px[i];
        }
        break;
      }
      default: error("Unsupported type '%s' passed to gsplit", type2char(tx));
      }
    }
    Free(count);
  }
  UNPROTECT(1);
  return res;
}

// This is for fmutate, to reorder the result of grouped data if the result has the same length as x
SEXP greorder(SEXP x, SEXP gobj) {
  if(TYPEOF(gobj) != VECSXP || !inherits(gobj, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
  if(LOGICAL(VECTOR_ELT(gobj, 5))[1] == 1) return x;
  const SEXP g = VECTOR_ELT(gobj, 1), gs = VECTOR_ELT(gobj, 2), order = VECTOR_ELT(gobj, 6);
  const int ng = length(gs), l = length(g), tx = TYPEOF(x),
            *pgs = INTEGER(gs), *pg = INTEGER(g);
  if(ng != INTEGER(VECTOR_ELT(gobj, 0))[0]) error("'GRP' object needs to have valid vector of group-sizes");
  if(l != length(x)) error("length(x) must match length(g)");

  SEXP res = PROTECT(allocVector(tx, l));

  // Note: This is only faster for a large number of groups...
  if(length(order) == l) { // Grouping not sorted but we have the ordering..
    SEXP sym_starts = PROTECT(install("starts"));
    const SEXP starts = getAttrib(order, sym_starts);
    UNPROTECT(1);
    if(length(starts) != ng) goto unsno2;
    const int *po = INTEGER(order), *ps = INTEGER(starts);

    switch(tx) {
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x), *pr = INTEGER(res);
      for(int i = 0, k = 0; i != ng; ++i) {
        for(int j = ps[i]-1, end = ps[i]+pgs[i]-1; j < end; ++j) pr[po[j]-1] = px[k++];
      }
      break;
    }
    case REALSXP: {
      double *px = REAL(x), *pr = REAL(res);
      for(int i = 0, k = 0; i != ng; ++i) {
        for(int j = ps[i]-1, end = ps[i]+pgs[i]-1; j < end; ++j) pr[po[j]-1] = px[k++];
      }
      break;
    }
    case CPLXSXP: {
      Rcomplex *px = COMPLEX(x), *pr = COMPLEX(res);
      for(int i = 0, k = 0; i != ng; ++i) {
        for(int j = ps[i]-1, end = ps[i]+pgs[i]-1; j < end; ++j) pr[po[j]-1] = px[k++];
      }
      break;
    }
    case STRSXP: {
      SEXP *px = STRING_PTR(x), *pr = STRING_PTR(res);
      for(int i = 0, k = 0; i != ng; ++i) {
        for(int j = ps[i]-1, end = ps[i]+pgs[i]-1; j < end; ++j) pr[po[j]-1] = px[k++];
      }
      break;
    }
    case VECSXP: {
      SEXP *px = SEXPPTR(x), *pr = SEXPPTR(res);
      for(int i = 0, k = 0; i != ng; ++i) {
        for(int j = ps[i]-1, end = ps[i]+pgs[i]-1; j < end; ++j) pr[po[j]-1] = px[k++];
      }
      break;
    }
    case RAWSXP: {
      Rbyte *px = RAW(x), *pr = RAW(res);
      for(int i = 0, k = 0; i != ng; ++i) {
        for(int j = ps[i]-1, end = ps[i]+pgs[i]-1; j < end; ++j) pr[po[j]-1] = px[k++];
      }
      break;
    }
    default: error("Unsupported type '%s' passed to gsplit", type2char(tx));
    }

  } else { // Unsorted, without ordering
    unsno2:;
    int *count = (int *) R_alloc(ng+1, sizeof(int));
    int *cgs = (int *) R_alloc(ng+2, sizeof(int)); cgs[1] = 0;
    for(int i = 0; i != ng; ++i) {
      count[i+1] = 0;
      cgs[i+2] = cgs[i+1] + pgs[i];
    }
    switch(tx) {
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x), *pr = INTEGER(res);
      for(int i = 0; i != l; ++i) pr[i] = px[cgs[pg[i]]+count[pg[i]]++];
      break;
    }
    case REALSXP: {
      double *px = REAL(x), *pr = REAL(res);
      for(int i = 0; i != l; ++i) pr[i] = px[cgs[pg[i]]+count[pg[i]]++];
      break;
    }
    case CPLXSXP: {
      Rcomplex *px = COMPLEX(x), *pr = COMPLEX(res);
      for(int i = 0; i != l; ++i) pr[i] = px[cgs[pg[i]]+count[pg[i]]++];
      break;
    }
    case STRSXP: {
      SEXP *px = STRING_PTR(x), *pr = STRING_PTR(res);
      for(int i = 0; i != l; ++i) pr[i] = px[cgs[pg[i]]+count[pg[i]]++];
      break;
    }
    case VECSXP: {
      SEXP *px = SEXPPTR(x), *pr = SEXPPTR(res);
      for(int i = 0; i != l; ++i) pr[i] = px[cgs[pg[i]]+count[pg[i]]++];
      break;
    }
    case RAWSXP: {
      Rbyte *px = RAW(x), *pr = RAW(res);
      for(int i = 0; i != l; ++i) pr[i] = px[cgs[pg[i]]+count[pg[i]]++];
      break;
    }
    default: error("Unsupported type '%s' passed to gsplit", type2char(tx));
    }
  }
  SHALLOW_DUPLICATE_ATTRIB(res, x);
  UNPROTECT(1);
  return res;
}

// Note: Only supports numeric data!!!!
SEXP lassign(SEXP x, SEXP s, SEXP rows, SEXP fill) {
  int l = length(x), tr = TYPEOF(rows), ss = asInteger(s), rs = LENGTH(rows);
  SEXP out = PROTECT(allocVector(VECSXP, l));
  // SEXP *px = VECTOR_PTR(x); // -> Depreciated interface: https://github.com/hadley/r-internals/blob/ea892fa79bbffe961e78dbe9c90ce4ca3bf2d9bc/vectors.md
  SEXP *px = SEXPPTR(x);
  double dfill = asReal(fill);

  if(tr == INTSXP) {
    int *rowsv = INTEGER(rows); //, vs = ss * sizeof(double);
    for(int j = l; j--; ) {
      SEXP column = px[j]; // VECTOR_ELT(x, j);
      if(length(column) != rs) error("length(rows) must match nrow(x)");
      SEXP outj;
      SET_VECTOR_ELT(out, j, outj = allocVector(REALSXP, ss));
      double *pcolumn = REAL(column), *poutj = REAL(outj);
      // memset(poutj, dfill, vs); // cannot memset missing values... can only memset 0
      for(int i = ss; i--; ) poutj[i] = dfill;
      for(int i = 0; i != rs; ++i) poutj[rowsv[i]-1] = pcolumn[i];
      SHALLOW_DUPLICATE_ATTRIB(outj, column);
    }
  } else if(tr == LGLSXP) {
    int *rowsv = LOGICAL(rows);
    if(ss != rs) error("length(rows) must match length(s) if rows is a logical vector");
    for(int j = l; j--; ) {
      SEXP column = px[j]; // VECTOR_ELT(x, j);
      SEXP outj;
      SET_VECTOR_ELT(out, j, outj = allocVector(REALSXP, ss));
      double *pcolumn = REAL(column), *poutj = REAL(outj);
      for(int i = 0, k = 0; i != rs; ++i) poutj[i] = rowsv[i] ? pcolumn[k++] : dfill;
      SHALLOW_DUPLICATE_ATTRIB(outj, column);
    }
  } else error("rows must be positive integers or a logical vector");
  SHALLOW_DUPLICATE_ATTRIB(out, x);
  UNPROTECT(1);
  return out;
}


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
  if(length(val) != 1) error("value needs to be length 1");
  int *buf = (int *) R_alloc(n, sizeof(int));
  SEXP ans;

  #define WHICHVLOOP                                             \
  if(invert) {                                                   \
    for(int i = 0; i != n; ++i) if(px[i] != v) buf[j++] = i+1;   \
  } else {                                                       \
    for(int i = 0; i != n; ++i) if(px[i] == v) buf[j++] = i+1;   \
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
  PROTECT(ans = allocVector(INTSXP, j));
  if(j) memcpy(INTEGER(ans), buf, sizeof(int) * j);

  UNPROTECT(1);
  return(ans);
}

SEXP anyallv(SEXP x, SEXP val, SEXP Rall) {

  int n = length(x), all = asLogical(Rall);
  if(length(val) != 1) error("value needs to be length 1");

  #define ALLANYVLOOP                                                    \
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

  int n = length(x), lv = length(val), lr = length(rep),
    ind1 = asLogical(Rind1), invert = asLogical(Rinvert), set = asLogical(Rset);

  if(lv > 1 || ind1) {
    if(lr != n) error("If length(v) > 1, length(r) must match length(x). Note that x[v] <- r is efficient in base R, only x[v] <- r[v] is optimized here.");
    if(TYPEOF(val) == LGLSXP) {
      if(lv != n) error("If v is a logical vector, length(v) needs to be equal to length(x). Note that x[v] <- r is efficient in base R, only x[v] <- r[v] is optimized here.");
    } else if(TYPEOF(val) == INTSXP) {
      if(invert) error("invert = TRUE is only possible if v is a logical vector");
    } else error("If length(v) > 1, v must be an integer or logical vector used to subset both x and r");
  } else if(lr != 1 && lr != n) error("If length(v) == 1, length(r) must be 1 or length(x)");


  SEXP ans = R_NilValue;
  if(set == 0) PROTECT(ans = duplicate(x)); // Fastest?? // copies attributes ?? -> Yes

  #define setcopyvLOOP(e)                                   \
  if(invert) {                                              \
    for(int i = 0; i != n; ++i) if(px[i] != v) px[i] = e;   \
  } else {                                                  \
    for(int i = 0; i != n; ++i) if(px[i] == v) px[i] = e;   \
  }

  #define setcopyvLOOPLVEC                                      \
  if(tv == INTSXP) {                                            \
    for(int i = 0; i != lv; ++i) px[pv[i]-1] = pr[pv[i]-1];     \
  } else if(invert == 0) {                                      \
    for(int i = 0; i != n; ++i) if(pv[i]) px[i] = pr[i];        \
  } else {                                                      \
    for(int i = 0; i != n; ++i) if(pv[i] == 0) px[i] = pr[i];   \
  }

  switch(TYPEOF(x)) {
  case INTSXP:
  case LGLSXP:
  {
    int *px = set ? INTEGER(x) : INTEGER(ans);
    if(lv == 1 && ind1 == 0) {
      int v;
      if(TYPEOF(val) == STRSXP) {
        if(!isFactor(x)) error("Type mismatch: if value is character, x must be character or factor.");
        v = fchmatch(getAttrib(x, R_LevelsSymbol), val, 0);
      } else v = asInteger(val);
      if(lr == 1) {
        const int r = asInteger(rep);
        setcopyvLOOP(r)
      } else {
        const int *pr = INTEGER(rep);
        setcopyvLOOP(pr[i])
      }
    } else {
      const int tv = TYPEOF(val), *pv = tv == INTSXP ? INTEGER(val) : LOGICAL(val), *pr = INTEGER(rep);
      setcopyvLOOPLVEC
    }
    break;
  }
  case REALSXP:
  {
    double *px = set ? REAL(x) : REAL(ans);
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
        const double *pr = REAL(rep);
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
      const int tv = TYPEOF(val), *pv = tv == INTSXP ? INTEGER(val) : LOGICAL(val);
      const double *pr = REAL(rep);
      setcopyvLOOPLVEC
    }
    break;
  }
  case STRSXP:
  {
    SEXP *px = set ? STRING_PTR(x) : STRING_PTR(ans);
    if(lv == 1 && ind1 == 0) {
      const SEXP v = asChar(val);
      if(lr == 1) {
        const SEXP r = asChar(rep);
        setcopyvLOOP(r)
      } else {
        const SEXP *pr = STRING_PTR(rep);
        setcopyvLOOP(pr[i])
      }
    } else {
      const int tv = TYPEOF(val), *pv = tv == INTSXP ? INTEGER(val) : LOGICAL(val);
      const SEXP *pr = STRING_PTR(rep);
      setcopyvLOOPLVEC
    }
    break;
  }
  case VECSXP:
  {
    SEXP *px = set ? SEXPPTR(x) : SEXPPTR(ans);
    if(lv == 1 && ind1 == 0) error("Cannot compare lists to a value");
    const int tv = TYPEOF(val), *pv = tv == INTSXP ? INTEGER(val) : LOGICAL(val);
    const SEXP *pr = SEXPPTR(rep);
    setcopyvLOOPLVEC
    break;
  }
  // case RAWSXP:
  // {
  //   Rbyte *px = set ? RAW(x) : RAW(ans);
  //   const Rbyte v = RAW(val)[0], r = RAW(rep)[0];
  //   setcopyvLOOP
  //   break;
  // }
  default: error("Unsupported type '%s' passed to setv() / copyv()", type2char(TYPEOF(x)));
  }
  if(set == 0) {
    UNPROTECT(1);
    return(ans);
  }
  return(x);
}

SEXP setop_core(SEXP x, SEXP val, SEXP op, SEXP roww) {

  int n = length(x), nv = length(val), o = asInteger(op), tx = TYPEOF(x);

    #define OPSWITCH(e)                                \
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

    #define OPSWITCHMAT(e)                             \
    switch(o) {                                        \
    case 1: for(int j = 0, cj; j != nc; ++j)  {        \
      cj = j * nr;                                     \
      for(int i = 0; i != nr; ++i) px[cj + i] += e;    \
      }                                                \
      break;                                           \
    case 2: for(int j = 0, cj; j != nc; ++j)  {        \
      cj = j * nr;                                     \
      for(int i = 0; i != nr; ++i) px[cj + i] -= e;    \
      }                                                \
      break;                                           \
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
        for(int i = 0; i != lx; ++i) setop_core(px[i], ScalarReal(pv[i]), op, roww);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *pv = INTEGER(val);
        for(int i = 0; i != lx; ++i) setop_core(px[i], ScalarInteger(pv[i]), op, roww);
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
    for(int i = 0; i != n; ++i) {
      SEXP ci = VECTOR_ELT(x, i);
      int tci = TYPEOF(ci);
      pans[i] = (tci == INTSXP || tci == REALSXP) && OBJECT(ci) == 0;
    }
    SET_TYPEOF(ans, LGLSXP);
    break;
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

// SEXP CasChar(SEXP x) {
//  return coerceVector(x, STRSXP);
// }

/* Inspired by:
 * do_list2env : .Internal(list2env(x, envir))
 */
SEXP multiassign(SEXP lhs, SEXP rhs, SEXP envir) {
  if(TYPEOF(lhs) != STRSXP) error("lhs needs to be character");
  int n = length(lhs);
  if(n == 1) { // lazy_duplicate appears not necessary (copy-on modify is automatically implemented, and <- also does not use it).
    defineVar(installChar(STRING_ELT(lhs, 0)), rhs, envir);
    return R_NilValue;
  }
  if(length(rhs) != n) error("length(lhs) must be equal to length(rhs)");
  SEXP *plhs = STRING_PTR(lhs);
  switch(TYPEOF(rhs)) { // installTrChar translates to native encoding, installChar does the same now, but also is available on older systems.
    case REALSXP: {
      double *prhs = REAL(rhs);
      for(int i = 0; i < n; ++i) defineVar(installChar(plhs[i]), ScalarReal(prhs[i]), envir);
      break;
    }
    case INTSXP: {
      int *prhs = INTEGER(rhs);
      for(int i = 0; i < n; ++i) defineVar(installChar(plhs[i]), ScalarInteger(prhs[i]), envir);
      break;
    }
    case STRSXP: {
      SEXP *prhs = STRING_PTR(rhs);
      for(int i = 0; i < n; ++i) defineVar(installChar(plhs[i]), ScalarString(prhs[i]), envir);
      break;
    }
    case LGLSXP: {
      int *prhs = LOGICAL(rhs);
      for(int i = 0; i < n; ++i) defineVar(installChar(plhs[i]), ScalarLogical(prhs[i]), envir);
      break;
    }
    case VECSXP: { // lazy_duplicate appears not necessary (copy-on modify is automatically implemented, and <- also does not use it).
      for(int i = 0; i < n; ++i) defineVar(installChar(plhs[i]), VECTOR_ELT(rhs, i), envir);
      break;
    }
    default: {
      SEXP rhsl = PROTECT(coerceVector(rhs, VECSXP));
      for(int i = 0; i < n; ++i) defineVar(installChar(plhs[i]), VECTOR_ELT(rhsl, i), envir);
      UNPROTECT(1);
    }
  }
  return R_NilValue;
}


SEXP vlabels(SEXP x, SEXP attrn, SEXP usenam) {
  if(!isString(attrn)) error("'attrn' must be of mode character");
  if(length(attrn) != 1) error("exactly one attribute 'attrn' must be given");
  SEXP sym_attrn = PROTECT(installChar(STRING_ELT(attrn, 0)));
  int l = length(x);
  if(TYPEOF(x) != VECSXP) {
    SEXP labx = getAttrib(x, sym_attrn);
    UNPROTECT(1);
    if(labx == R_NilValue) return ScalarString(NA_STRING);
    return labx;
  }
  SEXP res = PROTECT(allocVector(STRSXP, l));
  SEXP *pres = STRING_PTR(res), *px = SEXPPTR(x);
  for(int i = 0; i < l; ++i) {
    SEXP labxi = getAttrib(px[i], sym_attrn);
    pres[i] = labxi == R_NilValue ? NA_STRING : STRING_ELT(labxi, 0);
  }
  if(asLogical(usenam)) {
    SEXP nam = getAttrib(x, R_NamesSymbol);
    if(TYPEOF(nam) != NILSXP) namesgets(res, nam);
  }
  UNPROTECT(2);
  return res;
}

// Note: ind can be NULL...
SEXP setvlabels(SEXP x, SEXP attrn, SEXP value, SEXP ind) { // , SEXP sc
 if(!isString(attrn)) error("'attrn' must be of mode character");
 if(length(attrn) != 1) error("exactly one attribute 'attrn' must be given");
 if(TYPEOF(x) != VECSXP) error("X must be a list");
 int nprotect = 1, l = length(x), tv = TYPEOF(value); // , scl = asLogical(sc);
 SEXP *px = SEXPPTR(x); // , xsc;
 // if(scl) { // Create shallow copy
 //   if(INHERITS(x, char_datatable)) {
 //     xsc = PROTECT(Calloccol(x));
 //   } else {
 //     xsc = PROTECT(shallow_duplicate(x));
 //   }
 //   ++nprotect;
 //   px = SEXPPTR(xsc);
 // }
 SEXP *pv = px;
 if(tv != NILSXP) {
   if(tv == VECSXP || tv == STRSXP) {
    pv = SEXPPTR(value);
   } else {
    SEXP vl = PROTECT(coerceVector(value, VECSXP));
    pv = SEXPPTR(vl); ++nprotect;
   }
 }
 SEXP sym_attrn = PROTECT(installChar(STRING_ELT(attrn, 0)));
 if(length(ind) == 0) {
   if(tv != NILSXP && l != length(value)) error("length(x) must match length(value)");
   if(tv == NILSXP) {
     for(int i = 0; i < l; ++i) setAttrib(px[i], sym_attrn, R_NilValue);
   } else if(tv == STRSXP) {
     for(int i = 0; i < l; ++i) setAttrib(px[i], sym_attrn, ScalarString(pv[i]));
   } else {
     for(int i = 0; i < l; ++i) setAttrib(px[i], sym_attrn, pv[i]);
   }
 } else {
   if(TYPEOF(ind) != INTSXP) error("vlabels<-: ind must be of type integer");
   int li = length(ind), *pind = INTEGER(ind), ii;
   if(tv != NILSXP && li != length(value)) error("length(ind) must match length(value)");
   if(li == 0 || li > l) error("vlabels<-: length(ind) must be > 0 and <= length(x)");
   if(tv == NILSXP) {
     for(int i = 0; i < li; ++i) {
       ii = pind[i]-1;
       if(ii < 0 || ii >= l) error("vlabels<-: ind must be between 1 and length(x)");
       setAttrib(px[ii], sym_attrn, R_NilValue);
     }
   } else if(tv == STRSXP) {
     for(int i = 0; i < li; ++i) {
       ii = pind[i]-1;
       if(ii < 0 || ii >= l) error("vlabels<-: ind must be between 1 and length(x)");
       setAttrib(px[ii], sym_attrn, ScalarString(pv[i]));
     }
   } else {
     for(int i = 0; i < li; ++i) {
       ii = pind[i]-1;
       if(ii < 0 || ii >= l) error("vlabels<-: ind must be between 1 and length(x)");
       setAttrib(px[ii], sym_attrn, pv[i]);
     }
   }
 }
 UNPROTECT(nprotect);
 // return scl ? xsc : x;
 return x;
}


SEXP setnames(SEXP x, SEXP nam) {
  setAttrib(x, R_NamesSymbol, nam);
  return x;
}

