
// greorder by reference...

SEXP greorder(SEXP x, SEXP gobj) {
  if(TYPEOF(gobj) != VECSXP || !inherits(gobj, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
  if(LOGICAL(VECTOR_ELT(gobj, 5))[1] == 1) return x;
  const SEXP g = VECTOR_ELT(gobj, 1), gs = VECTOR_ELT(gobj, 2), order = VECTOR_ELT(gobj, 6);
  const int ng = length(gs), l = length(g), tx = TYPEOF(x),
            *pgs = INTEGER(gs), *pg = INTEGER(g);
  if(ng != INTEGER(VECTOR_ELT(gobj, 0))[0]) error("'GRP' object needs to have valid vector of group-sizes");
  if(l != length(x)) error("length(x) must match length(g)");

  // Note: This is only faster for a large number of groups...
  if(length(order) == l) { // Grouping not sorted but we have the ordering..
    const SEXP starts = getAttrib(order, install("starts"));
    if(length(starts) != ng) goto unsno2;
    const int *po = INTEGER(order), *ps = INTEGER(starts);

    #define GRORDLOOP_O                                               \
    for(int i = 0, k = 0, o; i != ng; ++i) {                          \
      for(int j = ps[i]-1, end = ps[i]+pgs[i]-1; j < end; ++j) {      \
        o = po[j]-1;                                                  \
        temp = px[o];                                                 \
        px[o] = px[k];                                                \
        px[k++] = temp;                                               \
      }                                                               \
    }
  // Old: not by reference
  // for(int i = 0, k = 0; i != ng; ++i) {
  //   for(int j = ps[i]-1, end = ps[i]+pgs[i]-1; j < end; ++j) pr[po[j]-1] = px[k++];
  // }

    switch(tx) {
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x), temp;
      GRORDLOOP_O
      break;
    }
    case REALSXP: {
      double *px = REAL(x), temp;
      GRORDLOOP_O
      break;
    }
    case CPLXSXP: {
      Rcomplex *px = COMPLEX(x), temp;
      GRORDLOOP_O
      break;
    }
    case STRSXP: {
      SEXP *px = STRING_PTR(x), temp;
      GRORDLOOP_O
      break;
    }
    case VECSXP: {
      SEXP *px = SEXPPTR(x), temp;
      GRORDLOOP_O
      break;
    }
    case RAWSXP: {
      Rbyte *px = RAW(x), temp;
      GRORDLOOP_O
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

    #define GRORDLOOP_UO              \
    for(int i = 0, o; i != l; ++i) {  \
      o = cgs[pg[i]]+count[pg[i]]++;  \
      temp = px[i];                   \
      px[i] = px[o];                  \
      px[o] = temp;                   \
    }
  // Old: not by reference
  // for(int i = 0; i != l; ++i) pr[i] = px[cgs[pg[i]]+count[pg[i]]++];

    switch(tx) {
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x), temp;
      GRORDLOOP_UO
      break;
    }
    case REALSXP: {
      double *px = REAL(x), temp;
      GRORDLOOP_UO
      break;
    }
    case CPLXSXP: {
      Rcomplex *px = COMPLEX(x), temp;
      GRORDLOOP_UO
      break;
    }
    case STRSXP: {
      SEXP *px = STRING_PTR(x), temp;
      GRORDLOOP_UO
      break;
    }
    case VECSXP: {
      SEXP *px = SEXPPTR(x), temp;
      GRORDLOOP_UO
      break;
    }
    case RAWSXP: {
      Rbyte *px = RAW(x), temp;
      GRORDLOOP_UO
      break;
    }
    default: error("Unsupported type '%s' passed to gsplit", type2char(tx));
    }
  }
  return x;
}
