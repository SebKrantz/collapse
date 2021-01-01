#include <R.h>
#include <Rinternals.h>

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

  for(int j = l; j--; ) { // This can go in any direction..
    SEXP column = VECTOR_ELT(x, j);
    int *pcolumn = INTEGER(column), jp = j+1;
    for(int i = pgs[j]; i--; ) pout[pcolumn[i]] = jp; // This can go in any direction...
  }
  UNPROTECT(1);
  return out;
}

// Note: Only supports numeric data!!!!
SEXP lassign(SEXP x, SEXP s, SEXP rows, SEXP fill) {
  int l = length(x), tr = TYPEOF(rows), ss = asInteger(s), rs = length(rows);
  SEXP out = PROTECT(allocVector(VECSXP, l));
  // SEXP *px = VECTOR_PTR(x); // -> Depreciated interface: https://github.com/hadley/r-internals/blob/ea892fa79bbffe961e78dbe9c90ce4ca3bf2d9bc/vectors.md
  double dfill = asReal(fill);

  if(tr == INTSXP) {
    int *rowsv = INTEGER(rows); //, vs = ss * sizeof(double);
    for(int j = l; j--; ) {
      SEXP column = VECTOR_ELT(x, j);
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
      SEXP column = VECTOR_ELT(x, j);
      SEXP outj;
      SET_VECTOR_ELT(out, j, outj = allocVector(REALSXP, ss));
      double *pcolumn = REAL(column), *poutj = REAL(outj);
      for(int i = 0, k = 0; i != rs; ++i) poutj[i] = rowsv[i] ? pcolumn[k++] : dfill;
      SHALLOW_DUPLICATE_ATTRIB(outj, column);
    }
  } else error("rows must be positive integers or a logical vector");
  DUPLICATE_ATTRIB(out, x);
  UNPROTECT(1);
  return out;
}

