#include "collapse_c.h"

static void eptrFinalizer(SEXP eptr) {
  if(!R_ExternalPtrAddr(eptr)) return;
  // R_SetExternalPtrProtected(eptr, R_NilValue);
  R_ClearExternalPtr(eptr);
}

SEXP createeptr(SEXP x) {
  SEXP eptr = PROTECT(R_MakeExternalPtr(x, R_NilValue, R_NilValue)); // x // Using the 'prot' or 'tag' fields includes the object in the pointer, which obscures the purpose of this which is memory efficiency.
  R_RegisterCFinalizerEx(eptr, eptrFinalizer, TRUE);
  UNPROTECT(1);
  return eptr;
}

SEXP geteptr(SEXP x) {
  if(TYPEOF(x) != EXTPTRSXP) return x;
  void * res = R_ExternalPtrAddr(x);
  if(!res) error("Invalid pointer to 'index': external pointers are only valid within the current R session. Please reindex() your data: data = reindex(data)");
  return (SEXP)res;
  // return R_ExternalPtrProtected(x);
}
