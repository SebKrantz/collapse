#include <R.h>
#include <Rinternals.h>

#define SEXPPTR(x) ((SEXP *)DATAPTR(x))  // to avoid overhead of looped VECTOR_ELT

void fcumsum_double_impl(double *pout, double *px, int ng, int *pg, int narm, int fill, int l) {
  if(ng == 0) {
    if(narm <= 0) {
      pout[0] = px[0];
      for(int i = 1; i != l; ++i) pout[i] = pout[i-1] + px[i];
    } else if(fill) {
      pout[0] = ISNAN(px[0]) ? 0.0 : px[0];
      for(int i = 1; i != l; ++i) pout[i] = pout[i-1] + (ISNAN(px[i]) ? 0.0 : px[i]);
    } else {
      double last = 0;
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i])) pout[i] = px[i];
        else last = pout[i] = px[i] + last;
      }
    }
  } else {
    double last[ng+1]; // Also pass pointer to function ??
    memset(last, 0.0, sizeof(double) * (ng+1));
    if(narm <= 0) {
      for(int i = 0; i != l; ++i) last[pg[i]] = pout[i] = last[pg[i]] + px[i];
    } else if(fill) {
      for(int i = 0; i != l; ++i) last[pg[i]] = pout[i] = last[pg[i]] + (ISNAN(px[i]) ? 0.0 : px[i]);
    } else {
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i])) pout[i] = px[i];
        else last[pg[i]] = pout[i] = last[pg[i]] + px[i];
      }
    }
  }
}

void fcumsum_double_impl_order(double *pout, double *px, int ng, int *pg, int *po, int narm, int fill, int l) {
  if(ng == 0) {
    if(narm <= 0) {
      pout[po[0]] = px[po[0]];
      for(int i = 1; i != l; ++i) pout[po[i]] = pout[po[i-1]] + px[po[i]];
    } else if(fill) {
      pout[po[0]] = ISNAN(px[po[0]]) ? 0.0 : px[po[0]];
      for(int i = 1; i != l; ++i) pout[po[i]] = pout[po[i-1]] + (ISNAN(px[po[i]]) ? 0.0 : px[po[i]]);
    } else {
      double last = 0;
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i];
        if(ISNAN(px[poi])) pout[poi] = px[poi];
        else last = pout[poi] = px[poi] + last;
      }
    }
  } else {
    double last[ng+1]; // Also pass pointer to function ??
    memset(last, 0.0, sizeof(double) * (ng+1));
    if(narm <= 0) {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i];
        last[pg[poi]] = pout[poi] = last[pg[poi]] + px[poi];
      }
    } else if(fill) {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i];
        last[pg[poi]] = pout[poi] = last[pg[poi]] + (ISNAN(px[poi]) ? 0.0 : px[poi]);
      }
    } else {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i];
        if(ISNAN(px[poi])) pout[poi] = px[poi];
        else last[pg[poi]] = pout[poi] = last[pg[poi]] + px[poi];
      }
    }
  }
}


void fcumsum_int_impl(int *pout, int *px, int ng, int *pg, int narm, int fill, int l) {
  if(ng == 0) {
    if(narm <= 0) {
      pout[0] = px[0];
      for(int i = 1; i != l; ++i) pout[i] = pout[i-1] + px[i];
    } else if(fill) {
      pout[0] = (px[0] == NA_INTEGER) ? 0 : px[0];
      for(int i = 1; i != l; ++i) pout[i] = pout[i-1] + (px[i] == NA_INTEGER ? 0 : px[i]);
    } else {
      for(int i = 0, last = 0; i != l; ++i) {
        if(px[i] == NA_INTEGER) pout[i] = NA_INTEGER;
        else last = pout[i] = px[i] + last;
      }
    }
  } else {
    int last[ng+1]; // Also pass pointer to function ??
    memset(last, 0, sizeof(int) * (ng+1));
    if(narm <= 0) {
      for(int i = 0; i != l; ++i) last[pg[i]] = pout[i] = last[pg[i]] + px[i];
    } else if(fill) {
      for(int i = 0; i != l; ++i) last[pg[i]] = pout[i] = last[pg[i]] + (px[i] == NA_INTEGER ? 0 : px[i]);
    } else {
      for(int i = 0; i != l; ++i) {
        if(px[i] == NA_INTEGER) pout[i] = NA_INTEGER;
        else last[pg[i]] = pout[i] = last[pg[i]] + px[i];
      }
    }
  }
}

void fcumsum_int_impl_order(int *pout, int *px, int ng, int *pg, int *po, int narm, int fill, int l) {
  if(ng == 0) {
    if(narm <= 0) {
      pout[po[0]] = px[po[0]];
      for(int i = 1; i != l; ++i) pout[po[i]] = pout[po[i-1]] + px[po[i]];
    } else if(fill) {
      pout[po[0]] = (px[po[0]] == NA_INTEGER) ? 0 : px[po[0]];
      for(int i = 1; i != l; ++i) pout[po[i]] = pout[po[i-1]] + (px[po[i]] == NA_INTEGER ? 0 : px[po[i]]);
    } else {
      for(int i = 0, last = 0, poi; i != l; ++i) {
        poi = po[i];
        if(px[poi] == NA_INTEGER) pout[poi] = NA_INTEGER;
        else last = pout[poi] = px[poi] + last;
      }
    }
  } else {
    int last[ng+1]; // Also pass pointer to function ??
    memset(last, 0, sizeof(int) * (ng+1));
    if(narm <= 0) {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i];
        last[pg[poi]] = pout[poi] = last[pg[poi]] + px[poi];
      }
    } else if(fill) {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i];
        last[pg[poi]] = pout[poi] = last[pg[poi]] + (px[poi] == NA_INTEGER ? 0 : px[poi]);
      }
    } else {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i];
        if(px[poi] == NA_INTEGER) pout[poi] = NA_INTEGER;
        else last[pg[poi]] = pout[poi] = last[pg[poi]] + px[poi];
      }
    }
  }
}

SEXP fcumsumC(SEXP x, SEXP Rng, SEXP g, SEXP o, SEXP Rnarm, SEXP Rfill) {
  int l = length(x), tx = TYPEOF(x), ng = asInteger(Rng),
    narm = asInteger(Rnarm), fill = asInteger(Rfill), *pg = INTEGER(g),
    ord  = length(o) > 1, *po = ord ? INTEGER(o) : pg;
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng > 0 && l != length(g)) error("length(g) must match length(x)");
  if(ord && l != length(o)) error("length(o) must match length(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(tx, l));
  switch(tx) {
  case REALSXP: {
    double *px = REAL(x), *pout = REAL(out);
    if(ord) fcumsum_double_impl_order(pout, px, ng, pg, po, narm, fill, l);
    else fcumsum_double_impl(pout, px, ng, pg, narm, fill, l);
    break;
  }
  case INTSXP: {
    int *px = INTEGER(x), *pout = INTEGER(out);
    if(ord) fcumsum_int_impl_order(pout, px, ng, pg, po, narm, fill, l);
    else fcumsum_int_impl(pout, px, ng, pg, narm, fill, l);
    break;
  }
  default: error("Unsupported SEXP type");
  }
  DUPLICATE_ATTRIB(out, x);
  UNPROTECT(1);
  return out;
}

SEXP fcumsummC(SEXP x, SEXP Rng, SEXP g, SEXP o, SEXP Rnarm, SEXP Rfill) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1],
     ng = asInteger(Rng), narm = asInteger(Rnarm), fill = asInteger(Rfill), *pg = INTEGER(g),
     ord  = length(o) > 1, *po = ord ? INTEGER(o) : pg;
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng > 0 && l != length(g)) error("length(g) must match nrow(x)");
  if(ord && l != length(o)) error("length(o) must match nrow(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(tx, l * col));
  switch(tx) {
  case REALSXP: {
    double *px = REAL(x), *pout = REAL(out);
    if(ord) for(int j = 0; j != col; ++j) fcumsum_double_impl_order(pout + j*l, px + j*l, ng, pg, po, narm, fill, l);
    else for(int j = 0; j != col; ++j) fcumsum_double_impl(pout + j*l, px + j*l, ng, pg, narm, fill, l);
    break;
  }
  case INTSXP: {
    int *px = INTEGER(x), *pout = INTEGER(out);
    if(ord) for(int j = 0; j != col; ++j) fcumsum_int_impl_order(pout + j*l, px + j*l, ng, pg, po, narm, fill, l);
    else for(int j = 0; j != col; ++j) fcumsum_int_impl(pout + j*l, px + j*l, ng, pg, narm, fill, l);
    break;
  }
  default: error("Unsupported SEXP type");
  }
  DUPLICATE_ATTRIB(out, x);
  UNPROTECT(1);
  return out;
}

SEXP fcumsumlC(SEXP x, SEXP Rng, SEXP g, SEXP o, SEXP Rnarm, SEXP Rfill) {
  int l = length(x);
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101
  SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out), *px = SEXPPTR(x);
  for(int j = 0; j != l; ++j) pout[j] = fcumsumC(px[j], Rng, g, o, Rnarm, Rfill);
  DUPLICATE_ATTRIB(out, x);
  UNPROTECT(1);
  return out;
}
