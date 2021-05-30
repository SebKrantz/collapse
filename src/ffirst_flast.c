
#include <R.h>
#include <Rinternals.h>
// #include <stdint.h>
#include <stdbool.h>

// TODO: Implemented smarter copy names ?!

#define SEXPPTR(x) ((SEXP *)DATAPTR(x))  // to avoid overhead of looped STRING_ELT and VECTOR_ELT

// About Pointers
// https://www.tutorialspoint.com/cprogramming/c_pointers.htm
// https://www.tutorialspoint.com/cprogramming/c_pointer_arithmetic.htm


// Use const ?
SEXP ffirst_impl(SEXP x, int ng, SEXP g, bool narm, int *gl) {

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
        SEXP *px = STRING_PTR(x);
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
        SEXP *px = SEXPPTR(x);
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
    copyMostAttrib(x, out); //  DUPLICATE_ATTRIB(out, x);
    if(getAttrib(x, R_NamesSymbol) != R_NilValue) {
      if(narm) namesgets(out, STRING_ELT(getAttrib(x, R_NamesSymbol), j)); // ScalarString()??
      else namesgets(out, STRING_ELT(getAttrib(x, R_NamesSymbol), 0)); // ScalarString()??
    }
    UNPROTECT(1);
    return out;
  } else { // with groups
    if(length(g) != l) error("length(g) must match nrow(X)");
    SEXP first = PROTECT(allocVector(tx, ng));
    if(narm) {
      int ngs = 0, *pg = INTEGER(g);
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pfirst = REAL(first);
        for(int i = l; i--; ) pfirst[i] = NA_REAL;
        --pfirst;
        for(int i = 0; i != l; ++i) {
          if(!ISNAN(px[i])) {
            if(ISNAN(pfirst[pg[i]])) {
              pfirst[pg[i]] = px[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *pfirst = STRING_PTR(first);
        for(int i = l; i--; ) pfirst[i] = NA_STRING;
        --pfirst;
        for(int i = 0; i != l; ++i) {
          if(px[i] != NA_STRING) {
            if(pfirst[pg[i]] == NA_STRING) {
              pfirst[pg[i]] = px[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pfirst = INTEGER(first);
        for(int i = l; i--; ) pfirst[i] = NA_INTEGER;
        --pfirst;
        for(int i = 0; i != l; ++i) {
          if(px[i] != NA_INTEGER) {
            if(pfirst[pg[i]] == NA_INTEGER) {
              pfirst[pg[i]] = px[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
        break;
      }
      case VECSXP: {
        SEXP *px = SEXPPTR(x), *pfirst = SEXPPTR(first);
        for(int i = l; i--; ) pfirst[i] = R_NilValue; // R_NilValue or just leave empty ??
        --pfirst;
        for(int i = 0; i != l; ++i) {
          if(length(px[i])) {
            if(pfirst[pg[i]] == R_NilValue) {
              pfirst[pg[i]] = px[i];
              ++ngs;
              if(ngs == ng) break;
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
      //            pfirst[pg[i]] = px[i];
      //            ++ngs;
      //            if(ngs == ng) break;
      //          }
      //        }
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pfirst = REAL(first);
        for(int i = ng; i--; ) pfirst[i] = px[gl[i]];
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *pfirst = STRING_PTR(first);
        for(int i = ng; i--; ) pfirst[i] = px[gl[i]];
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pfirst = INTEGER(first);
        for(int i = ng; i--; ) pfirst[i] = px[gl[i]];
        break;
      }
      case VECSXP: {
        SEXP *px = SEXPPTR(x), *pfirst = SEXPPTR(first);
        for(int i = ng; i--; ) pfirst[i] = px[gl[i]];
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    copyMostAttrib(x, first); // DUPLICATE_ATTRIB(first, x);
    UNPROTECT(1);
    return first;
  }
}

SEXP ffirstC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm) {
  int *pgl, ng = asInteger(Rng), narm = asLogical(Rnarm);
  if(ng == 0 || narm) {
    pgl = &ng; // TO avoid Wmaybe uninitialized
    return ffirst_impl(x, ng, g, narm, pgl);
  }

  // Using C-Array -> Error, possibly because unprotected ???
  int gl[ng], *pg = INTEGER(g), lg = length(g);
  memset(gl, 0, sizeof(int) * ng); //
  pgl = &gl[0] - 1; // Or gl-1; // Pointer to -1 array element (since g starts from 1): https://beginnersbook.com/2014/01/c-pointer-to-array-example/
  for(int i = 0; i != lg; ++i) if(!pgl[pg[i]]) pgl[pg[i]] = i; // Correct? even for first value ?

  //  SEXP gl = PROTECT(allocVector(INTSXP, ng));
  //  memset(gl, 0, sizeof(int)*ng); //
  //  int *pg = INTEGER(g);
  //  pgl = INTEGER(gl)-1; // Pointer to -1 array element (since g starts from 1): https://beginnersbook.com/2014/01/c-pointer-to-array-example/
  //  for(int i = length(g); i--; ) if(!pgl[pg[i]]) pgl[pg[i]] = i; // Correct? even for first value ?


  //  SEXP out = PROTECT(allocVector(INTSXP, ng));
  //  int *pout = INTEGER(out);
  //  for(int i = ng; i--; ) pout[i] = pgl[i+1];
  //  UNPROTECT(1);
  //  return out; // Checking pointer: appears to be correct...
  // UNPROTECT(1);
  // return gl;
  return ffirst_impl(x, ng, g, narm, ++pgl);
}

SEXP ffirstlC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm) {
  int l = length(x), *pgl, ng = asInteger(Rng), narm = asLogical(Rnarm), groups = ng != 0;
  if(groups && !narm) {
    int gl[ng], *pg = INTEGER(g), lg = length(g);
    memset(gl, 0, sizeof(int) * ng);
    pgl = &gl[0] - 1; // Or gl - 1; // Pointer to -1 array element (since g starts from 1): https://beginnersbook.com/2014/01/c-pointer-to-array-example/
    for(int i = 0; i != lg; ++i) if(!pgl[pg[i]]) pgl[pg[i]] = i; // Correct? even for first value ?
    ++pgl;
  } else pgl = &l; // To avoid Wmaybe uninitialized..
  SEXP out = PROTECT(allocVector(VECSXP, l));
  SEXP *px = SEXPPTR(x), *pout = SEXPPTR(out);
  // SEXP *px = SEXPPTR(x);
  //  return ffirst_impl(VECTOR_ELT(x, 0), ng, g, narm, pgl);
  for(int j = 0; j != l; ++j) pout[j] = ffirst_impl(px[j], ng, g, narm, pgl);
  //  // for(int j = l; j--; ) SET_VECTOR_ELT(out, j, ffirst_impl(VECTOR_ELT(x, j), ng, g, narm, pgl));
  DUPLICATE_ATTRIB(out, x);
  if(!groups) {
    setAttrib(out, R_RowNamesSymbol, ScalarInteger(1));
  } else {
    SEXP rn;
    setAttrib(out, R_RowNamesSymbol, rn = allocVector(INTSXP, 2)); // PROTECT ?? -> protected by out
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -ng;
  }
  UNPROTECT(1);
  return out;
}

// For matrix writing a separate function to increase efficiency.
SEXP ffirstmC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), ng = asInteger(Rng),
    narm = asLogical(Rnarm), drop = asInteger(Rdrop),
    l = INTEGER(dim)[0], col = INTEGER(dim)[1], end = l-1, nprotect = 1;
  if (l < 2) return x;
  if (ng == 0) {
    SEXP out = PROTECT(allocVector(tx, col));
    if(narm) {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int j = 0, i = 0; j != col; ++j) {
          while(ISNAN(px[i]) && i != end) ++i;
          pout[j] = px[i];
          px += l;
          i = 0;
        }
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *pout = STRING_PTR(out);
        for(int j = 0, i = 0; j != col; ++j) {
          while(px[i] == NA_STRING && i != end) ++i;
          pout[j] = px[i];
          px += l;
          i = 0;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int j = 0, i = 0; j != col; ++j) {
          while(px[i] == NA_INTEGER && i != end) ++i;
          pout[j] = px[i];
          px += l;
          i = 0;
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
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *pout = STRING_PTR(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l];
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l];
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    if(getAttrib(x, R_DimNamesSymbol) != R_NilValue) {
      SEXP dimnames = PROTECT(getAttrib(x, R_DimNamesSymbol)); ++nprotect; // PROTECT ??
      if(!isNull(VECTOR_ELT(dimnames, 1))) {
        if(drop) setAttrib(out, R_NamesSymbol, VECTOR_ELT(dimnames, 1));
        else {
          dim = PROTECT(allocVector(INTSXP, 2)); ++nprotect;
          INTEGER(dim)[0] = 1; INTEGER(dim)[1] = col;
          dimgets(out, dim);
          SET_VECTOR_ELT(dimnames, 0, R_NilValue);
          dimnamesgets(out, dimnames);
          if(!isObject(x)) copyMostAttrib(x, out);
        }
      } else if(!drop && !isObject(x)) copyMostAttrib(x, out);
    }
    UNPROTECT(nprotect);
    return out;
  } else { // with groups
    if(length(g) != l) error("length(g) must match nrow(X)");
    SEXP first = PROTECT(allocMatrix(tx, ng, col));
    int *pg = INTEGER(g);
    if(narm) {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pfirst = REAL(first);
        for(int i = l * col; i--; ) pfirst[i] = NA_REAL;
        --pfirst;
        for(int j = 0; j != col; ++j) {
          for(int i = 0; i != l; ++i) if(!ISNAN(px[i]) && ISNAN(pfirst[pg[i]])) pfirst[pg[i]] = px[i];
          px += l; pfirst += ng;
        }
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *pfirst = STRING_PTR(first);
        for(int i = l * col; i--; ) pfirst[i] = NA_STRING;
        --pfirst;
        for(int j = 0; j != col; ++j) {
          for(int i = 0; i != l; ++i) if(px[i] != NA_STRING && pfirst[pg[i]] == NA_STRING) pfirst[pg[i]] = px[i];
          px += l; pfirst += ng;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pfirst = INTEGER(first);
        for(int i = l * col; i--; ) pfirst[i] = NA_INTEGER;
        --pfirst;
        for(int j = 0; j != col; ++j) {
          for(int i = 0; i != l; ++i) if(px[i] != NA_INTEGER && pfirst[pg[i]] == NA_INTEGER) pfirst[pg[i]] = px[i];
          px += l; pfirst += ng;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    } else {
      int gl[ng], *pgl, lg = length(g);
      memset(gl, 0, sizeof(int) * ng);
      pgl = &gl[0] - 1; // Or gl - 1; // Pointer to -1 array element (since g starts from 1): https://beginnersbook.com/2014/01/c-pointer-to-array-example/
      for(int i = 0; i != lg; ++i) if(!pgl[pg[i]]) pgl[pg[i]] = i; // Correct? even for first value ?
      pgl++;
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pfirst = REAL(first);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pfirst[i] = px[pgl[i]];
          px += l; pfirst += ng;
        }
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *pfirst = STRING_PTR(first);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pfirst[i] = px[pgl[i]];
          px += l; pfirst += ng;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pfirst = INTEGER(first);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pfirst[i] = px[pgl[i]];
          px += l; pfirst += ng;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    if(getAttrib(x, R_DimNamesSymbol) != R_NilValue) {
      SEXP dimnames = getAttrib(x, R_DimNamesSymbol); // PROTECT ??
      if(!isNull(VECTOR_ELT(dimnames, 0))) SET_VECTOR_ELT(dimnames, 0, R_NilValue);
      dimnamesgets(first, dimnames);
    }
    if(!isObject(x)) copyMostAttrib(x, first);
    UNPROTECT(nprotect);
    return first;
  }
}


// Use const ?
SEXP flast_impl(SEXP x, int ng, SEXP g, bool narm, int *gl) {

  int l = length(x), tx = TYPEOF(x);
  if (l < 2) return x; // Prevents seqfault for numeric(0) #101
  if (ng == 0) {
    SEXP out = PROTECT(allocVector(tx, 1));
    int j = l-1;
    if(narm) {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x);
        while(ISNAN(px[j]) && j != 0) --j;
        REAL(out)[0] = px[j];
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x);
        while(px[j] == NA_STRING && j != 0) --j;
        SET_STRING_ELT(out, 0, px[j]);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        while(px[j] == NA_INTEGER && j != 0) --j;
        INTEGER(out)[0] = px[j];
        break;
      }
      case VECSXP: {
        SEXP *px = SEXPPTR(x);
        while(length(px[j]) == 0 && j != 0) --j;
        SET_VECTOR_ELT(out, 0, px[j]);
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    } else {
      switch(tx) {
      case REALSXP: REAL(out)[0] = REAL(x)[l-1];
        break;
      case STRSXP: SET_STRING_ELT(out, 0, STRING_ELT(x, l-1));
        break;
      case INTSXP:
      case LGLSXP: INTEGER(out)[0] = INTEGER(x)[l-1];
        break;
      case VECSXP: SET_VECTOR_ELT(out, 0, VECTOR_ELT(x, l-1));
        break;
      default: error("Unsupported SEXP type!");
      }
    }
    copyMostAttrib(x, out); //  DUPLICATE_ATTRIB(out, x);
    if(getAttrib(x, R_NamesSymbol) != R_NilValue) {
      if(narm) namesgets(out, STRING_ELT(getAttrib(x, R_NamesSymbol), j)); // ScalarString()??
      else namesgets(out, STRING_ELT(getAttrib(x, R_NamesSymbol), l-1)); // ScalarString()??
    }
    UNPROTECT(1);
    return out;
  } else { // with groups
    if(length(g) != l) error("length(g) must match nrow(X)");
    SEXP last = PROTECT(allocVector(tx, ng));
    if(narm) {
      int ngs = 0, *pg = INTEGER(g);
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *plast = REAL(last);
        for(int i = l; i--; ) plast[i] = NA_REAL;
        --plast;
        for(int i = l; i--; ) {
          if(!ISNAN(px[i])) {
            if(ISNAN(plast[pg[i]])) {
              plast[pg[i]] = px[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *plast = STRING_PTR(last);
        for(int i = l; i--; ) plast[i] = NA_STRING;
        --plast;
        for(int i = l; i--; ) {
          if(px[i] != NA_STRING) {
            if(plast[pg[i]] == NA_STRING) {
              plast[pg[i]] = px[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *plast = INTEGER(last);
        for(int i = l; i--; ) plast[i] = NA_INTEGER;
        --plast;
        for(int i = l; i--; ) {
          if(px[i] != NA_INTEGER) {
            if(plast[pg[i]] == NA_INTEGER) {
              plast[pg[i]] = px[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
        break;
      }
      case VECSXP: {
        SEXP *px = SEXPPTR(x), *plast = SEXPPTR(last);
        for(int i = l; i--; ) plast[i] = R_NilValue; // R_NilValue or just leave empty ??
        --plast;
        for(int i = l; i--; ) {
          if(length(px[i])) {
            if(plast[pg[i]] == R_NilValue) {
              plast[pg[i]] = px[i];
              ++ngs;
              if(ngs == ng) break;
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
      //        for(int i = l; i--; ) {
      //          if(gl[pg[i]]) {
      //            gl[pg[i]] = false;
      //            plast[pg[i]] = px[i];
      //            ++ngs;
      //            if(ngs == ng) break;
      //          }
      //        }
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *plast = REAL(last);
        for(int i = ng; i--; ) plast[i] = px[gl[i]];
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *plast = STRING_PTR(last);
        for(int i = ng; i--; ) plast[i] = px[gl[i]];
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *plast = INTEGER(last);
        for(int i = ng; i--; ) plast[i] = px[gl[i]];
        break;
      }
      case VECSXP: {
        SEXP *px = SEXPPTR(x), *plast = SEXPPTR(last);
        for(int i = ng; i--; ) plast[i] = px[gl[i]];
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    copyMostAttrib(x, last); // DUPLICATE_ATTRIB(last, x);
    UNPROTECT(1);
    return last;
  }
}

SEXP flastC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm) {
  int *pgl, ng = asInteger(Rng), narm = asLogical(Rnarm);
  if(ng == 0 || narm) {
    pgl = &ng; // TO avoid Wmaybe uninitialized
    return flast_impl(x, ng, g, narm, pgl);
  }

  // Using C-Array -> Error, possibly because unprotected ???
   int gl[ng], *pg = INTEGER(g);
   memset(gl, 0, sizeof(int) * ng); //
   pgl = &gl[0] - 1; // Or gl-1; // Pointer to -1 array element (since g starts from 1): https://beginnersbook.com/2014/01/c-pointer-to-array-example/
   for(int i = length(g); i--; ) if(!pgl[pg[i]]) pgl[pg[i]] = i; // Correct? even for last value ?

//  SEXP gl = PROTECT(allocVector(INTSXP, ng));
//  memset(gl, 0, sizeof(int)*ng); //
//  int *pg = INTEGER(g);
//  pgl = INTEGER(gl)-1; // Pointer to -1 array element (since g starts from 1): https://beginnersbook.com/2014/01/c-pointer-to-array-example/
//  for(int i = length(g); i--; ) if(!pgl[pg[i]]) pgl[pg[i]] = i; // Correct? even for last value ?


  //  SEXP out = PROTECT(allocVector(INTSXP, ng));
  //  int *pout = INTEGER(out);
  //  for(int i = ng; i--; ) pout[i] = pgl[i+1];
  //  UNPROTECT(1);
  //  return out; // Checking pointer: appears to be correct...
  // UNPROTECT(1);
  // return gl;
  return flast_impl(x, ng, g, narm, ++pgl);
}

SEXP flastlC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm) {
  int l = length(x), *pgl, ng = asInteger(Rng), narm = asLogical(Rnarm), groups = ng != 0;
  if(groups && !narm) {
    int gl[ng], *pg = INTEGER(g);
    memset(gl, 0, sizeof(int) * ng);
    pgl = &gl[0] - 1; // Or gl - 1; // Pointer to -1 array element (since g starts from 1): https://beginnersbook.com/2014/01/c-pointer-to-array-example/
    for(int i = length(g); i--; ) if(!pgl[pg[i]]) pgl[pg[i]] = i; // Correct? even for last value ?
    ++pgl;
  } else pgl = &l; // To avoid Wmaybe uninitialized..
  SEXP out = PROTECT(allocVector(VECSXP, l));
  SEXP *px = SEXPPTR(x), *pout = SEXPPTR(out);
  // SEXP *px = SEXPPTR(x);
  //  return flast_impl(VECTOR_ELT(x, 0), ng, g, narm, pgl);
  for(int j = 0; j != l; ++j) pout[j] = flast_impl(px[j], ng, g, narm, pgl);
  //  // for(int j = l; j--; ) SET_VECTOR_ELT(out, j, flast_impl(VECTOR_ELT(x, j), ng, g, narm, pgl));
  DUPLICATE_ATTRIB(out, x);
  if(!groups) {
    setAttrib(out, R_RowNamesSymbol, ScalarInteger(1));
  } else {
    SEXP rn;
    setAttrib(out, R_RowNamesSymbol, rn = allocVector(INTSXP, 2)); // PROTECT ?? -> protexted by out
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -ng;
  }
  UNPROTECT(1);
  return out;
}

// For matrix writing a separate function to increase efficiency.
SEXP flastmC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), ng = asInteger(Rng),
    narm = asLogical(Rnarm), drop = asInteger(Rdrop),
    l = INTEGER(dim)[0], col = INTEGER(dim)[1], nprotect = 1;
  if (l < 2) return x;
  if (ng == 0) {
    SEXP out = PROTECT(allocVector(tx, col));
    if(narm) {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int j = 0, i = l-1; j != col; ++j) {
          while(ISNAN(px[i]) && i != 0) --i;
          pout[j] = px[i];
          px += l;
          i = l-1;
        }
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *pout = STRING_PTR(out);
        for(int j = 0, i = l-1; j != col; ++j) {
          while(px[i] == NA_STRING && i != 0) --i;
          pout[j] = px[i];
          px += l;
          i = l-1;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int j = 0, i = l-1; j != col; ++j) {
          while(px[i] == NA_INTEGER && i != 0) --i;
          pout[j] = px[i];
          px += l;
          i = l-1;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    } else {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l + l-1];
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *pout = STRING_PTR(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l + l-1];
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l + l-1];
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    if(getAttrib(x, R_DimNamesSymbol) != R_NilValue) {
      SEXP dimnames = PROTECT(getAttrib(x, R_DimNamesSymbol)); ++nprotect; // PROTECT ??
      if(!isNull(VECTOR_ELT(dimnames, 1))) {
        if(drop) setAttrib(out, R_NamesSymbol, VECTOR_ELT(dimnames, 1));
        else {
          dim = PROTECT(allocVector(INTSXP, 2)); ++nprotect;
          INTEGER(dim)[0] = 1; INTEGER(dim)[1] = col;
          dimgets(out, dim);
          SET_VECTOR_ELT(dimnames, 0, R_NilValue);
          dimnamesgets(out, dimnames);
          if(!isObject(x)) copyMostAttrib(x, out);
        }
      } else if(!drop && !isObject(x)) copyMostAttrib(x, out);
    }
    UNPROTECT(nprotect);
    return out;
  } else { // with groups
    if(length(g) != l) error("length(g) must match nrow(X)");
    SEXP last = PROTECT(allocMatrix(tx, ng, col));
    int *pg = INTEGER(g);
    if(narm) {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *plast = REAL(last);
        for(int i = l * col; i--; ) plast[i] = NA_REAL;
        --plast;
        for(int j = 0; j != col; ++j) {
          for(int i = l; i--; ) if(!ISNAN(px[i]) && ISNAN(plast[pg[i]])) plast[pg[i]] = px[i];
          px += l; plast += ng;
        }
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *plast = STRING_PTR(last);
        for(int i = l * col; i--; ) plast[i] = NA_STRING;
        --plast;
        for(int j = 0; j != col; ++j) {
          for(int i = l; i--; ) if(px[i] != NA_STRING && plast[pg[i]] == NA_STRING) plast[pg[i]] = px[i];
          px += l; plast += ng;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *plast = INTEGER(last);
        for(int i = l * col; i--; ) plast[i] = NA_INTEGER;
        --plast;
        for(int j = 0; j != col; ++j) {
          for(int i = l; i--; ) if(px[i] != NA_INTEGER && plast[pg[i]] == NA_INTEGER) plast[pg[i]] = px[i];
          px += l; plast += ng;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    } else {
      int gl[ng], *pgl;
      memset(gl, 0, sizeof(int) * ng);
      pgl = &gl[0] - 1; // Or gl - 1; // Pointer to -1 array element (since g starts from 1): https://beginnersbook.com/2014/01/c-pointer-to-array-example/
      for(int i = length(g); i--; ) if(!pgl[pg[i]]) pgl[pg[i]] = i; // Correct? even for last value ?
      pgl++;
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *plast = REAL(last);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) plast[i] = px[pgl[i]];
          px += l; plast += ng;
        }
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x), *plast = STRING_PTR(last);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) plast[i] = px[pgl[i]];
          px += l; plast += ng;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *plast = INTEGER(last);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) plast[i] = px[pgl[i]];
          px += l; plast += ng;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    if(getAttrib(x, R_DimNamesSymbol) != R_NilValue) {
      SEXP dimnames = getAttrib(x, R_DimNamesSymbol); // PROTECT ??
      if(!isNull(VECTOR_ELT(dimnames, 0))) SET_VECTOR_ELT(dimnames, 0, R_NilValue);
      dimnamesgets(last, dimnames);
    }
    if(!isObject(x)) copyMostAttrib(x, last);
    UNPROTECT(nprotect);
    return last;
  }
}















