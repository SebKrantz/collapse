#ifdef _OPENMP
  #include <omp.h>
  #define OMP_NUM_PROCS omp_get_num_procs()
  #define OMP_THREAD_LIMIT omp_get_thread_limit()
  #define OMP_MAX_THREADS omp_get_max_threads()
#else
  #define OMP_NUM_PROCS 1
  #define OMP_THREAD_LIMIT 1
  #define OMP_MAX_THREADS 1
#endif

#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>

#define SEXPPTR(x) ((SEXP *)DATAPTR(x))  // to avoid overhead of looped VECTOR_ELT
#define SEXPPTR_RO(x) ((const SEXP *)DATAPTR_RO(x))  // to avoid overhead of looped VECTOR_ELT

#define NISNAN(x) ((x) == (x))  // opposite of ISNAN for doubles
// Faster than Rinternals version (which uses math library version)
#undef ISNAN
#define ISNAN(x) ((x) != (x))

extern int max_threads;

// from base_radixsort.h (with significant modifications)
SEXP Cradixsort(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
void num1radixsort(int *, Rboolean, Rboolean, SEXP);
void iradixsort(int *, Rboolean, Rboolean, int, int *);
void dradixsort(int *, Rboolean, Rboolean, int, double *);

// from stats_mAR.c
void multi_yw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
SEXP pacf1(SEXP, SEXP);

// from data.table.h (with major modifications)
SEXP collapse_init(SEXP);
SEXP dt_na(SEXP, SEXP, SEXP, SEXP);
SEXP allNAv(SEXP, SEXP);
SEXP frankds(SEXP, SEXP, SEXP, SEXP);
SEXP rbindlist(SEXP, SEXP, SEXP, SEXP);
SEXP setcolorder(SEXP, SEXP);
SEXP subsetDT(SEXP, SEXP, SEXP, SEXP);
SEXP subsetCols(SEXP, SEXP, SEXP);
SEXP subsetVector(SEXP, SEXP, SEXP);
void subsetVectorRaw(SEXP, SEXP, SEXP, const bool);
SEXP Calloccol(SEXP);
void writeValue(SEXP, SEXP, const int, const int);
void writeNA(SEXP, const int, const int);

// Native collapse functions
void matCopyAttr(SEXP out, SEXP x, SEXP Rdrop, int ng);
void DFcopyAttr(SEXP out, SEXP x, int ng);
SEXP falloc(SEXP, SEXP, SEXP);
SEXP frange(SEXP x, SEXP Rnarm);
SEXP fdist(SEXP x, SEXP vec, SEXP Rret, SEXP Rnthreads);
SEXP fnrowC(SEXP x);
// SEXP CasChar(SEXP x);
SEXP setAttributes(SEXP x, SEXP a);
SEXP setattributes(SEXP x, SEXP a);
// SEXP CsetAttr(SEXP object, SEXP a, SEXP v); -> mot more efficeint than attr i.e. for row.names...
// void setattr(SEXP x, SEXP a, SEXP v);
SEXP duplAttributes(SEXP x, SEXP y);
// void duplattributes(SEXP x, SEXP y);
// SEXP cond_duplAttributes(SEXP x, SEXP y);
SEXP CsetAttrib(SEXP object, SEXP a);
SEXP CcopyAttrib(SEXP to, SEXP from);
SEXP CcopyMostAttrib(SEXP to, SEXP from);
SEXP copyMostAttributes(SEXP to, SEXP from);
SEXP lassign(SEXP x, SEXP s, SEXP rows, SEXP fill);
SEXP groups2GRP(SEXP x, SEXP lx, SEXP gs);
SEXP gsplit(SEXP x, SEXP gobj, SEXP toint);
SEXP greorder(SEXP x, SEXP gobj);
SEXP Cna_rm(SEXP x);
SEXP whichv(SEXP x, SEXP val, SEXP Rinvert);
SEXP anyallv(SEXP x, SEXP val, SEXP Rall);
SEXP setcopyv(SEXP x, SEXP val, SEXP rep, SEXP Rinvert, SEXP Rset, SEXP Rind1);
SEXP setop(SEXP x, SEXP val, SEXP op, SEXP roww);
SEXP vtypes(SEXP x, SEXP isnum);
SEXP vlengths(SEXP x, SEXP usenam);
SEXP multiassign(SEXP lhs, SEXP rhs, SEXP envir);
SEXP vlabels(SEXP x, SEXP attrn, SEXP usenam);
SEXP setvlabels(SEXP x, SEXP attrn, SEXP value, SEXP ind);
SEXP setnames(SEXP x, SEXP nam);
SEXP Cissorted(SEXP x, SEXP strictly);
SEXP groupVec(SEXP X, SEXP starts, SEXP sizes);
SEXP groupAtVec(SEXP X, SEXP starts, SEXP naincl);
SEXP funiqueC(SEXP x);
SEXP fmatchC(SEXP x, SEXP table, SEXP nomatch, SEXP count, SEXP overid);
SEXP coerce_to_equal_types(SEXP x, SEXP table);
void count_match(SEXP res, int nt, int nmv);
SEXP createeptr(SEXP x);
SEXP geteptr(SEXP x);
SEXP fcrosscolon(SEXP x, SEXP ngp, SEXP y, SEXP ckna);
SEXP fwtabulate(SEXP x, SEXP w, SEXP ngp, SEXP ckna);
SEXP vecgcd(SEXP x);
SEXP all_funs(SEXP x);
SEXP unlock_collapse_namespace(SEXP env);
void writeValueByIndex(SEXP target, SEXP source, const int from, SEXP index);
SEXP pivot_long(SEXP data, SEXP ind, SEXP idcol);
SEXP pivot_wide(SEXP index, SEXP id, SEXP column, SEXP fill, SEXP Rnthreads);
SEXP sort_merge_join(SEXP x, SEXP table, SEXP ot, SEXP count);
SEXP replace_outliers(SEXP x, SEXP limits, SEXP value, SEXP single_limit, SEXP set);
SEXP multi_match(SEXP m, SEXP g);
// fnobs rewritten in C:
SEXP fnobsC(SEXP x, SEXP Rng, SEXP g);
SEXP fnobsmC(SEXP x, SEXP Rng, SEXP g, SEXP Rdrop);
SEXP fnobslC(SEXP x, SEXP Rng, SEXP g, SEXP Rdrop);
// ffirst and flast rewritten in C:
SEXP ffirstC(SEXP x, SEXP Rng, SEXP g, SEXP gst, SEXP Rnarm);
SEXP ffirstmC(SEXP x, SEXP Rng, SEXP g, SEXP gst, SEXP Rnarm, SEXP Rdrop);
SEXP ffirstlC(SEXP x, SEXP Rng, SEXP g, SEXP gst, SEXP Rnarm);
SEXP flastC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm);
SEXP flastmC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop);
SEXP flastlC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm);
// fsum rewritten in C:
SEXP fsumC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP fill, SEXP Rnthreads);
SEXP fsummC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP fill, SEXP Rdrop, SEXP Rnthreads);
SEXP fsumlC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP fill, SEXP Rdrop, SEXP Rnthreads);
// fprod rewritten in C:
SEXP fprodC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm);
SEXP fprodmC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop);
SEXP fprodlC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop);
// fmean rewritten in C:
SEXP fmeanC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rnthreads);
SEXP fmeanmC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads);
SEXP fmeanlC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads);
// fmin and fmax rewritten in C:
SEXP fminC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm);
SEXP fminmC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop);
SEXP fminlC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop);
SEXP fmaxC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm);
SEXP fmaxmC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop);
SEXP fmaxlC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop);
// Added fcumsum, written in C:
SEXP fcumsumC(SEXP x, SEXP Rng, SEXP g, SEXP o, SEXP Rnarm, SEXP Rfill);
SEXP fcumsummC(SEXP x, SEXP Rng, SEXP g, SEXP o, SEXP Rnarm, SEXP Rfill);
SEXP fcumsumlC(SEXP x, SEXP Rng, SEXP g, SEXP o, SEXP Rnarm, SEXP Rfill);
// TRA, rewritten in C and extended:
SEXP TRAC(SEXP x, SEXP xAG, SEXP g, SEXP Rret, SEXP Rset);
SEXP TRAmC(SEXP x, SEXP xAG, SEXP g, SEXP Rret, SEXP Rset);
SEXP TRAlC(SEXP x, SEXP xAG, SEXP g, SEXP Rret, SEXP Rset);
// fndistinct, rewritten in C:
SEXP fndistinctC(SEXP x, SEXP g, SEXP Rnarm, SEXP Rnthreads);
SEXP fndistinctlC(SEXP x, SEXP g, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads);
SEXP fndistinctmC(SEXP x, SEXP g, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads);
// fmode, rewritten in C:
SEXP fmodeC(SEXP x, SEXP g, SEXP w, SEXP Rnarm, SEXP Rret, SEXP Rnthreads);
SEXP fmodelC(SEXP x, SEXP g, SEXP w, SEXP Rnarm, SEXP Rret, SEXP Rnthreads);
SEXP fmodemC(SEXP x, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rret, SEXP Rnthreads);
// fnth, rewritten in C:
SEXP fnthC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rret, SEXP Rnthreads, SEXP o, SEXP checko);
SEXP fnthlC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rret, SEXP Rnthreads);
SEXP fnthmC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rret, SEXP Rnthreads);
// New: fquantile:
SEXP fquantileC(SEXP x, SEXP Rprobs, SEXP w, SEXP o, SEXP Rnarm, SEXP Rtype, SEXP Rnames, SEXP checko);
// Helper functions for C API
double dquickselect_elem(double *x, const int n, const unsigned int elem, double h);
double iquickselect_elem(int *x, const int n, const unsigned int elem, double h);
double dquickselect(double *x, const int n, const int ret, const double Q);
double iquickselect(int *x, const int n, const int ret, const double Q);
double nth_int(const int *restrict px, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q);
double nth_double(const double *restrict px, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q);
double nth_int_ord(const int *restrict px, const int *restrict po, int l, const int narm, const int ret, const double Q);
double nth_double_ord(const double *restrict px, const int *restrict po, int l, const int narm, const int ret, const double Q);
double w_nth_int_ord(const int *restrict px, const double *restrict pw, const int *restrict po, double h, int l, const int narm, const int ret, const double Q);
double w_nth_double_ord(const double *restrict px, const double *restrict pw, const int *restrict po, double h, int l, const int narm, const int ret, const double Q);
double w_nth_int_qsort(const int *restrict px, const double *restrict pw, const int *restrict po, double h,
                       const int l, const int sorted, const int narm, const int ret, const double Q);
double w_nth_double_qsort(const double *restrict px, const double *restrict pw, const int *restrict po, double h,
                          const int l, const int sorted, const int narm, const int ret, const double Q);
SEXP nth_impl(SEXP x, int narm, int ret, double Q);
SEXP nth_ord_impl(SEXP x, int *pxo, int narm, int ret, double Q);
SEXP w_nth_ord_impl(SEXP x, int *pxo, double *pw, int narm, int ret, double Q, double h);

