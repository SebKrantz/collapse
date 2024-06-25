// #include <Defn.h> // Not available in C API !!
// #include <Internal.h> // Not available in C API !!
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
// typedef uint64_t ZPOS64_T; // already defined in stdint.h
// #define ASCII_MASK (1<<6) // evaluates to 64 !!
// #define IS_ASCII(x) ((x)->sxpinfo.gp & ASCII_MASK)
// #define IS_ASCII(x) (LEVELS(x) & ASCII_MASK)

#define SEXPPTR(x) ((SEXP *)DATAPTR(x))

// NOTE: All of this is copied from Defn.h: https://github.com/wch/r-source/blob/28de75af0541f93832c5899139b969d290bf422e/src/include/Defn.h
// to avoid checking for ALTREP in TRUELENGTH, which slows down the code unnecessarily...

#define STDVEC_TRUELENGTH(x) (((VECSEXP) (x))->vecsxp.truelength)
#define SET_STDVEC_TRUELENGTH(x, v) (STDVEC_TRUELENGTH(x)=(v))

/* It would be better to find a way to avoid abusing TRUELENGTH, but
 in the meantime replace TRUELENGTH/SET_TRUELENGTH with
 TRLEN/SET_TRLEN that cast to int to avoid warnings. */
#define TRLEN(x) ((int) STDVEC_TRUELENGTH(x)) // ((int) TRUELENGTH(x))
#define SET_TRLEN(x, v) SET_STDVEC_TRUELENGTH(x, ((int) (v)))

#define MYLEV(x)	(((SEXPREC_partial *)(x))->sxpinfo.gp)
#define IS_ASCII(x) (MYLEV(x) & 64) // from data.table.h

#define SETTOF(x,v)	((((SEXPREC_partial *)(x))->sxpinfo.type)=(v))

// NOTE: All of this is copied from Defn.h: https://github.com/wch/r-source/blob/28de75af0541f93832c5899139b969d290bf422e/src/include/Defn.h
// to avoid checking for ALTREP in TRUELENGTH, which slows down the code unnecessarily...

#define NAMED_BITS 16

struct sxpinfo_struct {
  SEXPTYPE type      :  TYPE_BITS;
  /* ==> (FUNSXP == 99) %% 2^5 == 3 == CLOSXP
   * -> warning: `type' is narrower than values
   *              of its type
   * when SEXPTYPE was an enum */
  unsigned int scalar:  1;
  unsigned int obj   :  1;
  unsigned int alt   :  1;
  unsigned int gp    : 16;
  unsigned int mark  :  1;
  unsigned int debug :  1;
  unsigned int trace :  1;  /* functions and memory tracing */
  unsigned int spare :  1;  /* used on closures and when REFCNT is defined */
  unsigned int gcgen :  1;  /* old generation number */
  unsigned int gccls :  3;  /* node class */
  unsigned int named : NAMED_BITS;
  unsigned int extra : 32 - NAMED_BITS; /* used for immediate bindings */
}; /*		    Tot: 64 */

#define SEXPREC_HEADER           \
  struct sxpinfo_struct sxpinfo; \
  struct SEXPREC *attrib;        \
  struct SEXPREC *gengc_next_node, *gengc_prev_node

struct vecsxp_struct {
  R_xlen_t	length;
  R_xlen_t	truelength;
};

typedef struct VECTOR_SEXPREC {
  SEXPREC_HEADER;
  struct vecsxp_struct vecsxp;
} VECTOR_SEXPREC, *VECSEXP;

typedef struct {
  SEXPREC_HEADER;
} SEXPREC_partial;


void checkEncodings(SEXP x);
SEXP Cradixsort(SEXP NA_last, SEXP decreasing, SEXP RETstrt, SEXP RETgs, SEXP SORTStr, SEXP args);
void num1radixsort(int *o, Rboolean NA_last, Rboolean decreasing, SEXP x);
void iradixsort(int *o, Rboolean NA_last, Rboolean decreasing, int n, int *x);
void dradixsort(int *o, Rboolean NA_last, Rboolean decreasing, int n, double *x);
