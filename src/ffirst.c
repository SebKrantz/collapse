#include "collapse_c.h"
// #include <stdint.h>
// #include <stdbool.h>

// TODO: Implemented smarter copy names ?!
// About Pointers
// https://www.tutorialspoint.com/cprogramming/c_pointers.htm
// https://www.tutorialspoint.com/cprogramming/c_pointer_arithmetic.htm


// Use const ?
SEXP ffirst_impl(SEXP x, int ng, SEXP g, int narm, int *gl) {

  int l = length(x), tx = TYPEOF(x), end = l-1;
  if (l < 2) return x; // Prevents seqfault for numeric(0) #101
  if (ng == 0) {
    SEXP out = PROTECT(allocVector(tx, 1));
    int j = 0;
    if(narm) {
      switch(tx) {
        case REALSXP: {
          double *px = REAL(x);
          while(ISNAN(px[j]) && j != end) ++j;
          REAL(out)[0] = px[j];
          break;
        }
        case STRSXP: {
          const SEXP *px = SEXPPTR_RO(x);
          while(px[j] == NA_STRING && j != end) ++j;
          SET_STRING_ELT(out, 0, px[j]);
          break;
        }
        case INTSXP:
        case LGLSXP: {
          int *px = INTEGER(x);
          while(px[j] == NA_INTEGER && j != end) ++j;
          INTEGER(out)[0] = px[j];
          break;
        }
        case VECSXP: {
          const SEXP *px = SEXPPTR_RO(x);
          while(length(px[j]) == 0 && j != end) ++j;
          SET_VECTOR_ELT(out, 0, px[j]);
          break;
        }
        default: error("Unsupported SEXP type!");
      }
    } else {
      switch(tx) {
      case REALSXP: REAL(out)[0] = REAL(x)[0];
        break;
      case STRSXP: SET_STRING_ELT(out, 0, STRING_ELT(x, 0));
        break;
      case INTSXP:
      case LGLSXP: INTEGER(out)[0] = INTEGER(x)[0];
        break;
      case VECSXP: SET_VECTOR_ELT(out, 0, VECTOR_ELT(x, 0));
        break;
      default: error("Unsupported SEXP type!");
      }
    }
    if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")))
       copyMostAttrib(x, out);
    if(!isNull(getAttrib(x, R_NamesSymbol)))
      namesgets(out, ScalarString(STRING_ELT(getAttrib(x, R_NamesSymbol), j)));
    UNPROTECT(1);
    return out;
  } else { // with groups
    if(length(g) != l) error("length(g) must match nrow(X)");
    SEXP out = PROTECT(allocVector(tx, ng));
    if(narm) {
      int ngs = 0, *pg = INTEGER(g);
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int i = ng; i--; ) pout[i] = NA_REAL;
        --pout;
        for(int i = 0; i != l; ++i) {
          if(NISNAN(px[i])) { // Fastest ???
            if(ISNAN(pout[pg[i]])) {
              pout[pg[i]] = px[i];
              if(++ngs == ng) break;
            }
          }
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int i = ng; i--; ) pout[i] = NA_STRING;
        --pout;
        for(int i = 0; i != l; ++i) {
          if(px[i] != NA_STRING) {
            if(pout[pg[i]] == NA_STRING) {
              pout[pg[i]] = px[i];
              if(++ngs == ng) break;
            }
          }
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int i = ng; i--; ) pout[i] = NA_INTEGER;
        --pout;
        for(int i = 0; i != l; ++i) {
          if(px[i] != NA_INTEGER) {
            if(pout[pg[i]] == NA_INTEGER) {
              pout[pg[i]] = px[i];
              if(++ngs == ng) break;
            }
          }
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int i = ng; i--; ) pout[i] = R_NilValue; // R_NilValue or just leave empty ??
        --pout;
        for(int i = 0; i != l; ++i) {
          if(length(px[i])) {
            if(pout[pg[i]] == R_NilValue) {
              pout[pg[i]] = px[i];
              if(++ngs == ng) break;
            }
          }
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    } else {
      // Old Implementation: With boolean array
      //      bool gl[ng+1];
      //      memset(gl, 1, sizeof(bool) * (ng+1));
      //        for(int i = 0; i != l; ++i) {
      //          if(gl[pg[i]]) {
      //            gl[pg[i]] = false;
      //            pout[pg[i]] = px[i];
      //            ++ngs;
      //            if(ngs == ng) break;
      //          }
      //        }
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x)-1, *pout = REAL(out);
        for(int i = ng; i--; ) pout[i] = gl[i] == NA_INTEGER ? NA_REAL : px[gl[i]];
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x)-1, *pout = INTEGER(out);
        for(int i = ng; i--; ) pout[i] = gl[i] == NA_INTEGER ? NA_INTEGER : px[gl[i]];
        break;
      }
      case STRSXP:{
        const SEXP *px = SEXPPTR_RO(x)-1;
        SEXP *pout = SEXPPTR(out);
        for(int i = ng; i--; ) pout[i] = gl[i] == NA_INTEGER ? NA_STRING : px[gl[i]];
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x)-1;
        SEXP *pout = SEXPPTR(out);
        for(int i = ng; i--; ) pout[i] = gl[i] == NA_INTEGER ? R_NilValue : px[gl[i]];
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")))
       copyMostAttrib(x, out); // SHALLOW_DUPLICATE_ATTRIB(out, x);
    UNPROTECT(1);
    return out;
  }
}

SEXP ffirstC(SEXP x, SEXP Rng, SEXP g, SEXP gst, SEXP Rnarm) {
  int *pgl, ng = asInteger(Rng), narm = asLogical(Rnarm);
  if(ng == 0 || narm) {
    pgl = &ng; // TO avoid Wmaybe uninitialized
    return ffirst_impl(x, ng, g, narm, pgl);
  }
  if(length(gst) != ng) {
  // Using C-Array -> Not a good idea, variable length arrays give note on gcc11
  SEXP gl = PROTECT(allocVector(INTSXP, ng));
  int *pg = INTEGER(g), lg = length(g);
  pgl = INTEGER(gl);
  for(int i = ng; i--; ) pgl[i] = NA_INTEGER;
  --pgl; // &gl[0]-1 Or gl-1; // Pointer to -1 array element (since g starts from 1): https://beginnersbook.com/2014/01/c-pointer-to-array-example/
         // Above gives gcc11 issue !! (works with R INTEGER() pointer, not plain C array)
  for(int i = 0; i != lg; ++i) if(pgl[pg[i]] == NA_INTEGER) pgl[pg[i]] = i+1;

  //  SEXP gl = PROTECT(allocVector(INTSXP, ng));
  //  memset(gl, 0, sizeof(int)*ng); //
  //  int *pg = INTEGER(g);
  //  pgl = INTEGER(gl)-1; // Pointer to -1 array element (since g starts from 1): https://beginnersbook.com/2014/01/c-pointer-to-array-example/
  //  for(int i = length(g); i--; ) if(!pgl[pg[i]]) pgl[pg[i]] = i; // Correct? even for first value ?

  // SEXP out = PROTECT(allocVector(INTSXP, ng));
  // int *pout = INTEGER(out);
  // for(int i = ng; i--; ) pout[i] = pgl[i+1];
  // UNPROTECT(1);
  // return out; // Checking pointer: appears to be correct...
  // UNPROTECT(1);
  // return gl;
  SEXP res = ffirst_impl(x, ng, g, narm, ++pgl);
  UNPROTECT(1);
  return res;
  } else return ffirst_impl(x, ng, g, narm, INTEGER(gst));
}

SEXP ffirstlC(SEXP x, SEXP Rng, SEXP g, SEXP gst, SEXP Rnarm) {
  int l = length(x), *pgl, ng = asInteger(Rng), narm = asLogical(Rnarm), nprotect = 1;
  if(ng > 0 && !narm) {
    if(length(gst) != ng) {
    // Can't use integer array here because apparently it is removed by the garbage collector when passed to a new function
    SEXP gl = PROTECT(allocVector(INTSXP, ng)); ++nprotect;
    int *pg = INTEGER(g), lg = length(g); // gl[ng],
    pgl = INTEGER(gl); // pgl = &gl[0];
    for(int i = ng; i--; ) pgl[i] = NA_INTEGER;
    --pgl;
    for(int i = 0; i != lg; ++i) if(pgl[pg[i]] == NA_INTEGER) pgl[pg[i]] = i+1;
    ++pgl;
    } else pgl = INTEGER(gst);
  } else pgl = &l; // To avoid Wmaybe uninitialized..
  // return ffirst_impl(VECTOR_ELT(x, 0), ng, g, narm, pgl);
  SEXP out = PROTECT(allocVector(VECSXP, l));
  const SEXP *px = SEXPPTR_RO(x);
  for(int j = 0; j != l; ++j) SET_VECTOR_ELT(out, j, ffirst_impl(px[j], ng, g, narm, pgl));
  DFcopyAttr(out, x, ng);
  UNPROTECT(nprotect);
  return out;
}

// For matrix writing a separate function to increase efficiency.
SEXP ffirstmC(SEXP x, SEXP Rng, SEXP g, SEXP gst, SEXP Rnarm, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), ng = asInteger(Rng), narm = asLogical(Rnarm),
    l = INTEGER(dim)[0], col = INTEGER(dim)[1], end = l-1;
  if (l < 2) return x;
  if (ng == 0) {
    SEXP out = PROTECT(allocVector(tx, col));
    if(narm) {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int j = 0, i = 0; j != col; ++j) {
          while(ISNAN(px[i]) && i != end) ++i;
          pout[j] = px[i]; px += l; i = 0;
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int j = 0, i = 0; j != col; ++j) {
          while(px[i] == NA_STRING && i != end) ++i;
          pout[j] = px[i]; px += l; i = 0;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int j = 0, i = 0; j != col; ++j) {
          while(px[i] == NA_INTEGER && i != end) ++i;
          pout[j] = px[i]; px += l; i = 0;
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int j = 0, i = 0; j != col; ++j) {
          while(length(px[i]) == 0 && i != end) ++i;
          pout[j] = px[i]; px += l; i = 0;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    } else {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l];
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l];
        break;
      }
      case STRSXP:
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l];
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    matCopyAttr(out, x, Rdrop, ng);
    UNPROTECT(1);
    return out;
  } else { // with groups
    int nprotect = 1;
    if(length(g) != l) error("length(g) must match nrow(X)");
    SEXP out = PROTECT(allocVector(tx, ng * col));
    int *pg = INTEGER(g);
    if(narm) {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int i = ng * col; i--; ) pout[i] = NA_REAL;
        --pout;
        for(int j = 0; j != col; ++j) {
          for(int i = 0; i != l; ++i) if(NISNAN(px[i]) && ISNAN(pout[pg[i]])) pout[pg[i]] = px[i];
          px += l; pout += ng;
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int i = ng * col; i--; ) pout[i] = NA_STRING;
        --pout;
        for(int j = 0; j != col; ++j) {
          for(int i = 0; i != l; ++i) if(px[i] != NA_STRING && pout[pg[i]] == NA_STRING) pout[pg[i]] = px[i];
          px += l; pout += ng;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int i = ng * col; i--; ) pout[i] = NA_INTEGER;
        --pout;
        for(int j = 0; j != col; ++j) {
          for(int i = 0; i != l; ++i) if(px[i] != NA_INTEGER && pout[pg[i]] == NA_INTEGER) pout[pg[i]] = px[i];
          px += l; pout += ng;
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int i = ng * col; i--; ) pout[i] = R_NilValue;
        --pout;
        for(int j = 0; j != col; ++j) {
          for(int i = 0; i != l; ++i) if(length(px[i]) && pout[pg[i]] == R_NilValue) pout[pg[i]] = px[i];
          px += l; pout += ng;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    } else {
      int *pgl;
      if(length(gst) != ng) {
      SEXP gl = PROTECT(allocVector(INTSXP, ng)); ++nprotect;
      // int gl[ng], *pgl; pgl = &gl[0];
      pgl = INTEGER(gl);
      for(int i = ng; i--; ) pgl[i] = NA_INTEGER;
      --pgl; // gcc11 issue with plain array
      for(int i = 0; i != l; ++i) if(pgl[pg[i]] == NA_INTEGER) pgl[pg[i]] = i+1;
      ++pgl;
      } else pgl = INTEGER(gst);
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x)-1, *pout = REAL(out);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pout[i] = pgl[i] == NA_INTEGER ? NA_REAL : px[pgl[i]];
          px += l; pout += ng;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x)-1, *pout = INTEGER(out);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pout[i] = pgl[i] == NA_INTEGER ? NA_INTEGER : px[pgl[i]];
          px += l; pout += ng;
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = SEXPPTR_RO(x)-1;
        SEXP *pout = SEXPPTR(out);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pout[i] = pgl[i] == NA_INTEGER ? NA_STRING : px[pgl[i]];
          px += l; pout += ng;
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x)-1;
        SEXP *pout = SEXPPTR(out);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pout[i] = pgl[i] == NA_INTEGER ? R_NilValue : px[pgl[i]];
          px += l; pout += ng;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    matCopyAttr(out, x, Rdrop, ng);
    UNPROTECT(nprotect);
    return out;
  }
}

















