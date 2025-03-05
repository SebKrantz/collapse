#ifndef CONNECTION_H  // Check if CONNECTION_H is not defined
#define CONNECTION_H  // Define CONNECTION_H

#include <R.h>
#include <Rinternals.h>

#define DPTR DATAPTR
#define SEXPPTR(x) ((SEXP *)DATAPTR(x))  // to avoid overhead of looped VECTOR_ELT
#define SEXPPTR_RO(x) ((const SEXP *)DATAPTR_RO(x))  // to avoid overhead of looped VECTOR_ELT

#define MYLEV LEVELS
#define IS_UTF8(x)  (MYLEV(x) & 8)
#define IS_ASCII(x) (MYLEV(x) & 64) // from data.table.h
// #define ASCII_MASK (1<<6) // evaluates to 64 !!
// #define IS_ASCII(x) ((x)->sxpinfo.gp & ASCII_MASK)
// #define IS_ASCII(x) (LEVELS(x) & ASCII_MASK)

#define SETTOF SET_TYPEOF
#define SET_LEN SETLENGTH
#define TRULEN TRUELENGTH
#define SET_TRULEN SET_TRUELENGTH

// NOTE: All of this is copied from Defn.h: https://github.com/wch/r-source/blob/28de75af0541f93832c5899139b969d290bf422e/src/include/Defn.h

#ifndef SEXPREC_HEADER

#ifndef NAMED_BITS
#define NAMED_BITS 16
#endif

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

#endif

typedef struct {
  SEXPREC_HEADER;
} SEXPREC_partial;

// to avoid checking for ALTREP in TRUELENGTH, which slows down the code unnecessarily...
#ifndef STDVEC_TRUELENGTH
#define STDVEC_TRUELENGTH(x) (((VECSEXP) (x))->vecsxp.truelength)
#define SET_STDVEC_TRUELENGTH(x, v) (STDVEC_TRUELENGTH(x)=(v))
#endif
  /* It would be better to find a way to avoid abusing TRUELENGTH, but
   in the meantime replace TRUELENGTH/SET_TRUELENGTH with
   TRLEN/SET_TRLEN that cast to int to avoid warnings. */
#define TRLEN(x) ((int) STDVEC_TRUELENGTH(x)) // ((int) TRUELENGTH(x))
#define SET_TRLEN(x, v) SET_STDVEC_TRUELENGTH(x, ((int) (v)))

#define MYEFL(x) (((SEXPREC_partial *)(x))->sxpinfo.gp)
#define MYSEFL(x,v)	((((SEXPREC_partial *)(x))->sxpinfo.gp)=(v))


#endif // End of CONNECTION_H guard
