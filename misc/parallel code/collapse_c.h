// Copied from kit package...
#ifdef _OPENMP
#include <omp.h>
#define omp_enabled true
#define max_thread omp_get_num_procs()
#define min_thread 1
#define OMP_PARALLEL_FOR(nth) _Pragma("omp parallel for num_threads(nth)")
#else
#define omp_enabled false
#define max_thread 1
#define min_thread 1
#define omp_get_thread_num() 0
#define OMP_PARALLEL_FOR(n)
#endif

#include <R.h>
#include <Rinternals.h>

#define SEXPPTR(x) ((SEXP *)DATAPTR(x))  // to avoid overhead of looped VECTOR_ELT
#define NISNAN(x) ((x) == (x))  // opposite of ISNAN for doubles
// Faster than Rinternals version (which uses math library version)
#undef ISNAN
#define ISNAN(x) ((x) != (x))

void matCopyAttr(SEXP out, SEXP x, SEXP Rdrop, int ng);
void DFcopyAttr(SEXP out, SEXP x, int ng);
