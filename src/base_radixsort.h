// #include <Defn.h> // Not available in C API !!
// #include <Internal.h> // Not available in C API !!
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
// typedef uint64_t ZPOS64_T; // already defined in stdint.h
// #define ASCII_MASK (1<<6) // evaluates to 64 !!
// #define IS_ASCII(x) ((x)->sxpinfo.gp & ASCII_MASK)
// #define IS_ASCII(x) (LEVELS(x) & ASCII_MASK)

// NOTE: All of this is copied from Defn.h: https://github.com/wch/r-source/blob/28de75af0541f93832c5899139b969d290bf422e/src/include/Defn.h
// to avoid checking for ALTREP in TRUELENGTH, which slows down the code unnecessarily...

#define STDVEC_TRUELENGTH(x) (((VECSEXP) (x))->vecsxp.truelength)
#define SET_STDVEC_TRUELENGTH(x, v) (STDVEC_TRUELENGTH(x)=(v))

/* It would be better to find a way to avoid abusing TRUELENGTH, but
 in the meantime replace TRUELENGTH/SET_TRUELENGTH with
 TRLEN/SET_TRLEN that cast to int to avoid warnings. */
#define TRLEN(x) ((int) STDVEC_TRUELENGTH(x)) // ((int) TRUELENGTH(x))
#define SET_TRLEN(x, v) SET_STDVEC_TRUELENGTH(x, ((int) (v)))

void checkEncodings(SEXP x);


SEXP Cradixsort(SEXP NA_last, SEXP decreasing, SEXP RETstrt, SEXP RETgs, SEXP SORTStr, SEXP args);
void num1radixsort(int *o, Rboolean NA_last, Rboolean decreasing, SEXP x);
void iradixsort(int *o, Rboolean NA_last, Rboolean decreasing, int n, int *x);
void dradixsort(int *o, Rboolean NA_last, Rboolean decreasing, int n, double *x);
