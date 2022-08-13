#include <R.h>
#include <Rinternals.h>

#define SEXPPTR(x) ((SEXP *)DATAPTR(x))  // to avoid overhead of looped VECTOR_ELT
#define NISNAN(x) ((x) == (x))  // opposite of ISNAN for doubles
// Faster than Rinternals version (which uses math library version)
#undef ISNAN
#define ISNAN(x) ((x) != (x))

void matCopyAttr(SEXP out, SEXP x, SEXP Rdrop, int ng);
void DFcopyAttr(SEXP out, SEXP x, int ng);

void multi_yw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
SEXP collapse_init(SEXP);
SEXP dt_na(SEXP, SEXP);
SEXP allNAv(SEXP, SEXP);
SEXP Cradixsort(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP frankds(SEXP, SEXP, SEXP, SEXP);
SEXP pacf1(SEXP, SEXP);
SEXP rbindlist(SEXP, SEXP, SEXP, SEXP);
SEXP setcolorder(SEXP, SEXP);
SEXP subsetDT(SEXP, SEXP, SEXP, SEXP);
SEXP subsetCols(SEXP, SEXP, SEXP);
SEXP subsetVector(SEXP, SEXP, SEXP);
SEXP Calloccol(SEXP dt); // , SEXP Rn
SEXP falloc(SEXP, SEXP);
SEXP frange(SEXP x, SEXP Rnarm);
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
SEXP createeptr(SEXP x);
SEXP geteptr(SEXP x);
SEXP fcrosscolon(SEXP x, SEXP ngp, SEXP y, SEXP ckna);
SEXP fwtabulate(SEXP x, SEXP w, SEXP ngp, SEXP ckna);
SEXP vecgcd(SEXP x);
SEXP all_funs(SEXP x);
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
SEXP fsumC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP fill, SEXP Rnth);
SEXP fsummC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP fill, SEXP Rdrop, SEXP Rnth);
SEXP fsumlC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP fill, SEXP Rdrop, SEXP Rnth);
// fprod rewritten in C:
SEXP fprodC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm);
SEXP fprodmC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop);
SEXP fprodlC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop);
// fmean rewritten in C:
SEXP fmeanC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rnth);
SEXP fmeanmC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnth);
SEXP fmeanlC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnth);
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

