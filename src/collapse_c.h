#include <R.h>
#include <Rinternals.h>

#define SEXPPTR(x) ((SEXP *)DATAPTR(x))  // to avoid overhead of looped VECTOR_ELT
#define NISNAN(x) ((x) == (x))  // opposite of ISNAN for doubles

void matCopyAttr(SEXP out, SEXP x, SEXP Rdrop, int ng);
void DFcopyAttr(SEXP out, SEXP x, int ng);
