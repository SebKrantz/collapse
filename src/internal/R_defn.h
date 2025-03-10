#ifndef R_DEFINITIONS_H  // Check if R_DEFINITIONS_H is not defined
#define R_DEFINITIONS_H  // Define R_DEFINITIONS_H

// #define USE_RINTERNALS

#include <R.h>
#include <Rinternals.h>

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

  struct vecsxp_struct {
    R_xlen_t	length;
    R_xlen_t	truelength;
  };

  struct primsxp_struct {
    int offset;
  };

  struct symsxp_struct {
    struct SEXPREC *pname;
    struct SEXPREC *value;
    struct SEXPREC *internal;
  };

  struct listsxp_struct {
    struct SEXPREC *carval;
    struct SEXPREC *cdrval;
    struct SEXPREC *tagval;
  };

  struct envsxp_struct {
    struct SEXPREC *frame;
    struct SEXPREC *enclos;
    struct SEXPREC *hashtab;
  };

  struct closxp_struct {
    struct SEXPREC *formals;
    struct SEXPREC *body;
    struct SEXPREC *env;
  };

  struct promsxp_struct {
    struct SEXPREC *value;
    struct SEXPREC *expr;
    struct SEXPREC *env;
  };


#define SEXPREC_HEADER           \
  struct sxpinfo_struct sxpinfo; \
  struct SEXPREC *attrib;        \
  struct SEXPREC *gengc_next_node, *gengc_prev_node

typedef struct SEXPREC {
  SEXPREC_HEADER;
  union {
    struct primsxp_struct primsxp;
    struct symsxp_struct symsxp;
    struct listsxp_struct listsxp;
    struct envsxp_struct envsxp;
    struct closxp_struct closxp;
    struct promsxp_struct promsxp;
  } u;
} SEXPREC;

// typedef struct {
//   SEXPREC_HEADER;
// } SEXPREC_partial;

typedef struct VECTOR_SEXPREC {
  SEXPREC_HEADER;
  struct vecsxp_struct vecsxp;
} VECTOR_SEXPREC, *VECSEXP;

typedef union { VECTOR_SEXPREC s; double align; } SEXPREC_ALIGN;

#endif

#undef OOBJ
#define OOBJ(x)	((x)->sxpinfo.obj)
#define SET_OOBJ(x,v) (OOBJ(x)=(v))
#undef ATTTR
#define ATTTR(x)	((x)->attrib)
#define SET_ATTTR(x,v) (ATTR(x)=(v))

#undef MYLEV
#define MYLEV(x)	((x)->sxpinfo.gp)
#undef IS_UTF8
#define IS_UTF8(x)  (MYLEV(x) & 8)
#undef IS_ASCII
#define IS_ASCII(x) (MYLEV(x) & 64) // from data.table.h
// #define ASCII_MASK (1<<6) // evaluates to 64 !!
// #define IS_ASCII(x) ((x)->sxpinfo.gp & ASCII_MASK)
// #define IS_ASCII(x) (LEVELS(x) & ASCII_MASK)
#undef SETTOF
#define SETTOF(x,v)	(((x)->sxpinfo.type)=(v))

// to avoid checking for ALTREP in TRUELENGTH, which slows down the code unnecessarily...
#ifndef STDVEC_TRUELENGTH
#define STDVEC_TRUELENGTH(x) (((VECSEXP) (x))->vecsxp.truelength)
#define SET_STDVEC_TRUELENGTH(x, v) (STDVEC_TRUELENGTH(x)=(v))
#endif
  /* It would be better to find a way to avoid abusing TRUELENGTH, but
   in the meantime replace TRUELENGTH/SET_TRUELENGTH with
   TRLEN/SET_TRLEN that cast to int to avoid warnings. */

#undef TRULEN
#define TRULEN(x) (ALTREP(x) ? 0 : STDVEC_TRUELENGTH(x))
#undef SET_TRULEN
#define SET_TRULEN(x, v) (STDVEC_TRUELENGTH(x)=(v))

#undef TRLEN
#define TRLEN(x) ((int) STDVEC_TRUELENGTH(x)) // ((int) TRUELENGTH(x))
#undef SET_TRLEN
#define SET_TRLEN(x, v) SET_STDVEC_TRUELENGTH(x, ((int) (v)))

#ifndef STDVEC_LENGTH
#define STDVEC_LENGTH(x) (((VECSEXP) (x))->vecsxp.length)
#endif
  // Needed for SETLENGTH
#ifndef SETSCAL
#define SETSCAL(x, v) (((x)->sxpinfo.scalar) = (v))
#endif

#ifndef SET_STDVEC_LENGTH
#define SET_STDVEC_LENGTH(x,v) do {		      \
  SEXP __x__ = (x);			                    \
  R_xlen_t __v__ = (v);			                \
  STDVEC_LENGTH(__x__) = __v__;		          \
  SETSCAL(__x__, __v__ == 1 ? 1 : 0);	      \
} while (0)
#endif

#undef SET_LEN
#define SET_LEN(x, v) SET_STDVEC_LENGTH((x), (v))

#undef MYEFL
#define MYEFL(x) ((x)->sxpinfo.gp)
#undef MYSEFL
#define MYSEFL(x,v)	(((x)->sxpinfo.gp)=(v))

// For super efficient access, e.g. in gsplit()
#undef SEXP_DATAPTR
#define SEXP_DATAPTR(x) ((SEXP *) (((SEXPREC_ALIGN *) (x)) + 1))

#undef DPTR
#define DPTR(x) ((void *)DATAPTR_RO(x))
#undef SEXPPTR
#define SEXPPTR(x) ((SEXP *)DATAPTR_RO(x))  // to avoid overhead of looped VECTOR_ELT
#undef SEXPPTR_RO
#define SEXPPTR_RO(x) ((const SEXP *)DATAPTR_RO(x))  // to avoid overhead of looped VECTOR_ELT

// #define STDVEC_DATAPTR(x) ((void *) (((SEXPREC_ALIGN *) (x)) + 1))
//
// static R_INLINE void *DPTR(SEXP x) {
//   if (ALTREP(x)) error("Cannot get writable DATAPTR from ALTREP string or list");
//   else if (LENGTH(x) == 0 && TYPEOF(x) != CHARSXP) return (void *) 1;
//   else return STDVEC_DATAPTR(x);
// }
// External symbols not in DLL?
// extern inline void *DPTR(SEXP x) {
//   return DATAPTR(x);
// }

#endif // End of R_DEFINITIONS_H guard
