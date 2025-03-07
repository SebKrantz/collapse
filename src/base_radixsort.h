// #include <Defn.h> // Not available in C API !!
// #include <Internal.h> // Not available in C API !!
// #define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include "internal/R_defn.h"
// typedef uint64_t ZPOS64_T; // already defined in stdint.h

void checkEncodings(SEXP x);
SEXP Cradixsort(SEXP NA_last, SEXP decreasing, SEXP RETstrt, SEXP RETgs, SEXP SORTStr, SEXP args);
void num1radixsort(int *o, Rboolean NA_last, Rboolean decreasing, SEXP x);
void iradixsort(int *o, Rboolean NA_last, Rboolean decreasing, int n, int *x);
void dradixsort(int *o, Rboolean NA_last, Rboolean decreasing, int n, double *x);
