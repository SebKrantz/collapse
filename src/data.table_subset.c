/*
 This code is adapted from the data.table package: http://r-datatable.com
 and licensed under a Mozilla Public License 2.0 (MPL-2.0) license.
*/

#include "data.table.h"

// selfref stuff is taken from data.tables assign.c

static void finalizer(SEXP p)
{
  SEXP x;
  R_len_t n, l, tl;
  if(!R_ExternalPtrAddr(p)) error("Internal error: finalizer hasn't received an ExternalPtr"); // # nocov
  p = R_ExternalPtrTag(p);
  if (!isString(p)) error("Internal error: finalizer's ExternalPtr doesn't see names in tag"); // # nocov
  l = LENGTH(p);
  tl = TRUELENGTH(p);
  if (l<0 || tl<l) error("Internal error: finalizer sees l=%d, tl=%d",l,tl); // # nocov
  n = tl-l;
  if (n==0) {
    // gc's ReleaseLargeFreeVectors() will have reduced R_LargeVallocSize by the correct amount
    // already, so nothing to do (but almost never the case).
    return;
  }
  x = PROTECT(allocVector(INTSXP, 50));  // 50 so it's big enough to be on LargeVector heap. See NodeClassSize in memory.c:allocVector
  // INTSXP rather than VECSXP so that GC doesn't inspect contents after LENGTH (thanks to Karl Miller, Jul 2015)
  SETLENGTH(x,50+n*2*sizeof(void *)/4);  // 1*n for the names, 1*n for the VECSXP itself (both are over allocated).
  UNPROTECT(1);
  return;
}

void setselfref(SEXP x) {
  SEXP p;
  // Store pointer to itself so we can detect if the object has been copied. See
  // ?copy for why copies are not just inefficient but cause a problem for over-allocated data.tables.
  // Called from C only, not R level, so returns void.
  // setAttrib(x, SelfRefSymbol, R_NilValue); // Probably not needed but I like it here.
  setAttrib(x, SelfRefSymbol, p=R_MakeExternalPtr(
    R_NilValue,                  // for identical() to return TRUE. identical() doesn't look at tag and prot
    PROTECT(getAttrib(x, R_NamesSymbol)),  // to detect if names has been replaced and its tl lost, e.g. setattr(DT,"names",...)
    PROTECT(R_MakeExternalPtr(   // to avoid an infinite loop in object.size(), if prot=x here
        x,                         // to know if this data.table has been copied by key<-, attr<-, names<-, etc.
        R_NilValue,                // this tag and prot currently unused
        R_NilValue
    ))
  ));
  R_RegisterCFinalizerEx(p, finalizer, FALSE);
  UNPROTECT(2);

  /*
   *  base::identical doesn't check prot and tag of EXTPTR, just that the ptr itself is the
   same in both objects. R_NilValue is always equal to R_NilValue.  R_NilValue is a memory
   location constant within an R session, but can vary from session to session. So, it
   looks like a pointer to a user looking at attributes(DT), but they might wonder how it
   works if they realise the selfref of all data.tables all point to the same address (rather
   than to the table itself which would be reasonable to expect given the attribute's name).
   *  p=NULL rather than R_NilValue works too, other than we need something other than NULL
   so we can detect tables loaded from disk (which set p to NULL, see 5.13 of R-exts).
   *  x is wrapped in another EXTPTR because of object.size (called by tables(), and by users).
   If the prot (or tag) was x directly it sends object.size into an infinite loop and then
   "segfault from C stack overflow" (object.size does count tag and prot, unlike identical,
   but doesn't count what's pointed to).
   *  Could use weak reference possibly, but the fact that they can get be set to R_NilValue
   by gc (?) didn't seem appropriate.
   *  If the .internal.selfref attribute is removed (e.g. by user code), nothing will break, but
   an extra copy will just be taken on next :=, with warning, with a new selfref.
   *  object.size will count size of names twice, but that's ok as only small.
   *  Thanks to Steve L for suggesting ExtPtr for this, rather than the previous REALSXP
   vector which required data.table to do a show/hide dance in a masked identical.
   */
}

// also need this stuff from assign.c
static SEXP shallow(SEXP dt, SEXP cols, R_len_t n)
{
  // NEW: cols argument to specify the columns to shallow copy on. If NULL, all columns.
  // called from alloccol where n is checked carefully, or from shallow() at R level
  // where n is set to truelength (i.e. a shallow copy only with no size change)
  int protecti=0;
  SEXP newdt = PROTECT(allocVector(VECSXP, n)); protecti++;   // to do, use growVector here?
  SET_ATTRIB(newdt, shallow_duplicate(ATTRIB(dt)));
  SET_OBJECT(newdt, OBJECT(dt));
  IS_S4_OBJECT(dt) ? SET_S4_OBJECT(newdt) : UNSET_S4_OBJECT(newdt);  // To support S4 objects that incude data.table
  //SHALLOW_DUPLICATE_ATTRIB(newdt, dt);  // SHALLOW_DUPLICATE_ATTRIB would be a bit neater but is only available from R 3.3.0

  // TO DO: keepattr() would be faster, but can't because shallow isn't merely a shallow copy. It
  //        also increases truelength. Perhaps make that distinction, then, and split out, but marked
  //        so that the next change knows to duplicate.
  //        keepattr() also merely points to the entire attrbutes list and thus doesn't allow replacing
  //        some of its elements.

  // We copy all attributes that refer to column names so that calling setnames on either
  // the original or the shallow copy doesn't break anything.
  SEXP index = PROTECT(getAttrib(dt, sym_index)); protecti++;
  setAttrib(newdt, sym_index, shallow_duplicate(index));

  SEXP sorted = PROTECT(getAttrib(dt, sym_sorted)); protecti++;
  setAttrib(newdt, sym_sorted, duplicate(sorted));

  SEXP names = PROTECT(getAttrib(dt, R_NamesSymbol)); protecti++;
  SEXP newnames = PROTECT(allocVector(STRSXP, n)); protecti++;

  SEXP *pdt = SEXPPTR(dt), *pnewdt = SEXPPTR(newdt),
    *pnam = STRING_PTR(names), *pnnam = STRING_PTR(newnames);

  const int l = isNull(cols) ? LENGTH(dt) : length(cols);
  if (isNull(cols)) {
    for (int i=0; i != l; ++i) pnewdt[i] = pdt[i];
    if (length(names)) {
      if (length(names) < l) error("Internal error: length(names)>0 but <length(dt)"); // # nocov
      for (int i=0; i != l; ++i) pnnam[i] = pnam[i];
    }
    // else an unnamed data.table is valid e.g. unname(DT) done by ggplot2, and .SD may have its names cleared in dogroups, but shallow will always create names for data.table(NULL) which has 100 slots all empty so you can add to an empty data.table by reference ok.
  } else {
    int *pcols = INTEGER(cols);
    for (int i=0; i != l; ++i) pnewdt[i] = pdt[pcols[i]-1];
    if (length(names)) {
      // no need to check length(names) < l here. R-level checks if all value
      // in 'cols' are valid - in the range of 1:length(names(x))
      for (int i=0; i != l; ++i) pnnam[i] = pnam[pcols[i]-1];
    }
  }
  setAttrib(newdt, R_NamesSymbol, newnames);
  // setAttrib appears to change length and truelength, so need to do that first _then_ SET next,
  // otherwise (if the SET were were first) the 100 tl is assigned to length.
  SETLENGTH(newnames,l);
  SET_TRUELENGTH(newnames,n);
  SETLENGTH(newdt,l);
  SET_TRUELENGTH(newdt,n);
  setselfref(newdt);
  UNPROTECT(protecti);
  return(newdt);
}

// Can allocate conditionally on size, for export... use in collap, qDT etc.
SEXP Calloccol(SEXP dt) // , SEXP Rn
{
  R_len_t tl, n, l;
  l = LENGTH(dt);
  n = l + 100; // asInteger(GetOption1(sym_collapse_DT_alloccol));  // asInteger(Rn);
  tl = TRUELENGTH(dt);
  // R <= 2.13.2 and we didn't catch uninitialized tl somehow
  if (tl < 0) error("Internal error, tl of class is marked but tl<0."); // # nocov
  // better disable these...
  // if (tl > 0 && tl < l) error("Internal error, please report (including result of sessionInfo()) to collapse issue tracker: tl (%d) < l (%d) but tl of class is marked.", tl, l); // # nocov
  // if (tl > l+10000) warning("tl (%d) is greater than 10,000 items over-allocated (l = %d). If you didn't set the collapse_DT_alloccol option to be very large, please report to collapse issue tracker including the result of sessionInfo().",tl,l);

// TODO:  MAKE THIS WORK WITHOUT SHALLOW COPYING EVERY TIME !!!

// if (n > tl)
    return shallow(dt, R_NilValue, n); // usual case (increasing alloc)

//  SEXP nam = PROTECT(getAttrib(dt, R_NamesSymbol));
//  if(LENGTH(nam) != l) SETLENGTH(nam, l);
//  SET_TRUELENGTH(nam, n);
//  setselfref(dt); // better, otherwise may be invalid !!
//  UNPROTECT(1);
//  return(dt);
}
// #pragma GCC diagnostic ignored "-Wunknown-pragmas" // don't display this warning!! // https://stackoverflow.com/questions/1867065/how-to-suppress-gcc-warnings-from-library-headers?noredirect=1&lq=1

static void subsetVectorRaw(SEXP ans, SEXP source, SEXP idx, const bool anyNA)
// Only for use by subsetDT() or subsetVector() below, hence static
{
  const int n = length(idx);
  if (length(ans)!=n) error("Internal error: subsetVectorRaw length(ans)==%d n=%d", length(ans), n);

  const int *restrict idxp = INTEGER(idx);
  // anyNA refers to NA _in idx_; if there's NA in the data (source) that's just regular data to be copied
  // negatives, zeros and out-of-bounds have already been dealt with in convertNegAndZero so we can rely
  // here on idx in range [1,length(ans)].
  //  _Pragma("omp parallel for num_threads(getDTthreads())") (in PARLOOP below)
  //  _Pragma("omp parallel for num_threads(getDTthreads())")

  #define PARLOOP(_NAVAL_)                                        \
  if (anyNA) {                                                    \
    for (int i = 0; i != n; ++i) {                                \
      int elem = idxp[i];                                         \
      ap[i] = elem==NA_INTEGER ? _NAVAL_ : sp[elem];              \
    }                                                             \
  } else {                                                        \
    for (int i = 0; i != n; ++i) {                                \
      ap[i] = sp[idxp[i]];                                        \
    }                                                             \
  }
  // For small n such as 2,3,4 etc we hope OpenMP will be sensible inside it and not create a team with each thread doing just one item. Otherwise,
  // call overhead would be too high for highly iterated calls on very small subests. Timings were tested in #3175
  // Futher, we desire (currently at least) to stress-test the threaded code (especially in latest R-devel) on small data to reduce chance that bugs
  // arise only over a threshold of n.

  switch(TYPEOF(source)) {
    case INTSXP:
    case LGLSXP: {
      int *restrict sp = INTEGER(source)-1, *restrict ap = INTEGER(ans);
      PARLOOP(NA_INTEGER);
    } break;
    case REALSXP : {
      if (INHERITS(source, char_integer64)) {
        int64_t *restrict sp = (int64_t *)REAL(source)-1, *restrict ap = (int64_t *)REAL(ans);
        PARLOOP(INT64_MIN);
      } else {
        double *restrict sp = REAL(source)-1, *restrict ap = REAL(ans);
        PARLOOP(NA_REAL);
      }
    } break;
    case STRSXP : {
      // write barrier (assigning strings/lists) is not thread safe. Hence single threaded.
      // To go parallel here would need access to NODE_IS_OLDER, at least. Given gcgen, mark and named
      // are upper bounded and max 3, REFCNT==REFCNTMAX could be checked first and then critical SET_ if not.
      // Inside that critical just before SET_ it could check REFCNT<REFCNTMAX still held. Similarly for gcgen.
      // TODO - discuss with Luke Tierney. Produce benchmarks on integer/double to see if it's worth making a safe
    //        API interface for package use for STRSXP.
    // Aside: setkey() is a separate special case (a permutation) and does do this in parallel without using SET_*.
    SEXP *restrict sp = STRING_PTR(source)-1, *restrict ap = STRING_PTR(ans);
    PARLOOP(NA_STRING);
  } break;
  case VECSXP : {
    // VECTOR_PTR does exist but returns 'not safe to return vector pointer' when USE_RINTERNALS is not defined.
    // VECTOR_DATA and LIST_POINTER exist too but call VECTOR_PTR. All are clearly not intended to be used by packages.
    // The concern is overhead inside VECTOR_ELT() biting when called repetitively in a loop like we do here. That's why
      // we take the R API (INTEGER()[i], REAL()[i], etc) outside loops for the simple types even when not parallel. For this
      // type list case (VECSXP) it might be that some items are ALTREP for example, so we really should use the heavier
      // _ELT accessor (VECTOR_ELT) inside the loop in this case.
      SEXP *restrict sp = SEXPPTR(source)-1, *restrict ap = SEXPPTR(ans);
      PARLOOP(R_NilValue);
    } break;
    case CPLXSXP : {
      Rcomplex *restrict sp = COMPLEX(source)-1, *restrict ap = COMPLEX(ans);
      PARLOOP(NA_CPLX);
    } break;
    case RAWSXP : {
      Rbyte *restrict sp = RAW(source)-1, *restrict ap = RAW(ans);
      PARLOOP(0);
    } break;
    default :
      error("Internal error: column type '%s' not supported by data.table subset. All known types are supported so please report as bug.", type2char(TYPEOF(source)));  // # nocov
  }
}

static const char *check_idx(SEXP idx, int max, bool *anyNA_out) // , bool *orderedSubset_out) Not needed
// set anyNA for branchless subsetVectorRaw
// error if any negatives, zeros or >max since they should have been dealt with by convertNegAndZeroIdx() called ealier at R level.
// single cache efficient sweep with prefetch, so very low priority to go parallel
{
  if (!isInteger(idx)) error("Internal error. 'idx' is type '%s' not 'integer'", type2char(TYPEOF(idx))); // # nocov
    bool anyNA = false; // anyLess=false,
  // int last = INT32_MIN;
  int *idxp = INTEGER(idx), n = LENGTH(idx);
  for (int i = 0; i != n; ++i) {
    int elem = idxp[i];
    if (elem<=0 && elem!=NA_INTEGER) return "Internal inefficiency: idx contains negatives or zeros. Should have been dealt with earlier.";  // e.g. test 762  (TODO-fix)
    if (elem>max) return "Internal inefficiency: idx contains an item out-of-range. Should have been dealt with earlier.";                   // e.g. test 1639.64
    anyNA |= elem==NA_INTEGER;
    // anyLess |= elem<last;
    // last = elem;
  }
  *anyNA_out = anyNA;
  // *orderedSubset_out = !anyLess; // for the purpose of ordered keys elem==last is allowed
  return NULL;
}

SEXP convertNegAndZeroIdx(SEXP idx, SEXP maxArg, SEXP allowOverMax)
{
  // called from [.data.table to massage user input, creating a new strictly positive idx if there are any negatives or zeros
                  // + more precise and helpful error messages telling user exactly where the problem is (saving user debugging time)
                  // + a little more efficient than negativeSubscript in src/main/subscript.c (it's private to R so we can't call it anyway)
                  // allowOverMaxArg is false when := (test 1024), otherwise true for selecting

                  if (!isInteger(idx)) error("Internal error. 'idx' is type '%s' not 'integer'", type2char(TYPEOF(idx))); // # nocov
                    if (!isInteger(maxArg) || length(maxArg)!=1) error("Internal error. 'maxArg' is type '%s' and length %d, should be an integer singleton", type2char(TYPEOF(maxArg)), length(maxArg)); // # nocov
                    if (!isLogical(allowOverMax) || LENGTH(allowOverMax)!=1 || LOGICAL(allowOverMax)[0]==NA_LOGICAL) error("Internal error: allowOverMax must be TRUE/FALSE");  // # nocov
                    int max = INTEGER(maxArg)[0], n=LENGTH(idx);
                  if (max<0) error("Internal error. max is %d, must be >= 0.", max); // # nocov    includes NA which will print as INT_MIN
                    int *idxp = INTEGER(idx);

                  bool stop = false;
                  // #pragma omp parallel for num_threads(getDTthreads())
                  for (int i = 0; i != n; ++i) {
                    if (stop) continue;
                    int elem = idxp[i];
                    if ((elem<1 && elem!=NA_INTEGER) || elem>max) stop=true;
                  }
                  if (!stop) return(idx); // most common case to return early: no 0, no negative; all idx either NA or in range [1-max]

                  // ---------
                    // else massage the input to a standard idx where all items are either NA or in range [1,max] ...

                  int countNeg=0, countZero=0, countNA=0, firstOverMax=0;
                  for (int i = 0; i != n; ++i) {
                    int elem = idxp[i];
                    if (elem==NA_INTEGER) countNA++;
                    else if (elem<0) countNeg++;
                    else if (elem==0) countZero++;
                    else if (elem>max && firstOverMax==0) firstOverMax=i+1;
                  }
                  if (firstOverMax && LOGICAL(allowOverMax)[0]==FALSE) {
                    error("i[%d] is %d which is out of range [1,nrow=%d]", firstOverMax, idxp[firstOverMax-1], max);
                  }

                  int countPos = n-countNeg-countZero-countNA;
                  if (countPos && countNeg) {
                    int i = 0, firstNeg=0, firstPos=0;
                    while (i != n && (firstNeg==0 || firstPos==0)) {
                      int elem = idxp[i];
                      if (firstPos==0 && elem>0) firstPos=i+1;
                      if (firstNeg==0 && elem<0 && elem!=NA_INTEGER) firstNeg=i+1;
                      i++;
                    }
                    error("Item %d of i is %d and item %d is %d. Cannot mix positives and negatives.", firstNeg, idxp[firstNeg-1], firstPos, idxp[firstPos-1]);
                  }
                  if (countNeg && countNA) {
                    int i = 0, firstNeg=0, firstNA=0;
                    while (i != n && (firstNeg==0 || firstNA==0)) {
                      int elem = idxp[i];
                      if (firstNeg==0 && elem<0 && elem!=NA_INTEGER) firstNeg=i+1;
                      if (firstNA==0 && elem==NA_INTEGER) firstNA=i+1;
                      i++;
                    }
                    error("Item %d of i is %d and item %d is NA. Cannot mix negatives and NA.", firstNeg, idxp[firstNeg-1], firstNA);
                  }

                  SEXP ans;
                  if (countNeg==0) {
                    // just zeros to remove, or >max to convert to NA
                    ans = PROTECT(allocVector(INTSXP, n - countZero));
                    int *ansp = INTEGER(ans);
                    for (int i = 0, ansi = 0; i != n; ++i) {
                      int elem = idxp[i];
                      if (elem==0) continue;
                      ansp[ansi++] = elem>max ? NA_INTEGER : elem;
                    }
                  } else {
                    // idx is all negative without any NA but perhaps some zeros
                    bool *keep = (bool *)R_alloc(max, sizeof(bool));    // 4 times less memory that INTSXP in src/main/subscript.c
                    for (int i = 0; i != max; ++i) keep[i] = true;
                    int countRemoved=0, countDup=0, countBeyond=0;   // idx=c(-10,-5,-10) removing row 10 twice
                    int firstBeyond=0, firstDup=0;
                    for (int i = 0; i != n; ++i) {
                      int elem = -idxp[i];
                      if (elem==0) continue;
                      if (elem>max) {
                        countBeyond++;
                        if (firstBeyond==0) firstBeyond=i+1;
                        continue;
                      }
                      if (!keep[elem-1]) {
                        countDup++;
                        if (firstDup==0) firstDup=i+1;
                      } else {
                        keep[elem-1] = false;
                        countRemoved++;
                      }
                    }
                    if (countBeyond)
                      warning("Item %d of i is %d but there are only %d rows. Ignoring this and %d more like it out of %d.", firstBeyond, idxp[firstBeyond-1], max, countBeyond-1, n);
                    if (countDup)
                      warning("Item %d of i is %d which removes that item but that has occurred before. Ignoring this dup and %d other dups.", firstDup, idxp[firstDup-1], countDup-1);
                    int ansn = max-countRemoved;
                    ans = PROTECT(allocVector(INTSXP, ansn));
                    int *ansp = INTEGER(ans);
                    for (int i = 0, ansi = 0; i != max; ++i) {
                      if (keep[i]) ansp[ansi++] = i+1;
                    }
                  }
                  UNPROTECT(1);
                  return ans;
}

static void checkCol(SEXP col, int colNum, int nrow, SEXP x)
{
  if (isNull(col)) error("Column %d is NULL; malformed data.table.", colNum);
  if (isNewList(col) && INHERITS(col, char_dataframe)) {
    SEXP names = getAttrib(x, R_NamesSymbol);
    error("Column %d ['%s'] is a data.frame or data.table; malformed data.table.",
          colNum, isNull(names)?"":CHAR(STRING_ELT(names,colNum-1)));
  }
  if (length(col)!=nrow) {
    SEXP names = getAttrib(x, R_NamesSymbol);
    error("Column %d ['%s'] is length %d but column 1 is length %d; malformed data.table.",
          colNum, isNull(names)?"":CHAR(STRING_ELT(names,colNum-1)), length(col), nrow);
  }
}


/* helper */
SEXP extendIntVec(SEXP x, int len, int val) {
  SEXP out = PROTECT(allocVector(INTSXP, len + 1));
  int *pout = INTEGER(out), *px = INTEGER(x);
  for(int i = len; i--; ) pout[i] = px[i];
  pout[len] = val;
  UNPROTECT(1);
  return out;
}


/* subset columns of a list efficiently */

SEXP subsetCols(SEXP x, SEXP cols, SEXP checksf) { // SEXP fretall
  if(TYPEOF(x) != VECSXP) error("x is not a list.");
  int l = LENGTH(x), nprotect = 3, oxl = OBJECT(x) != 0;
  if(l == 0) return x; //  ncol == 0 -> Nope, need emty selections such as cat_vars(mtcars) !!
  PROTECT_INDEX ipx;
  PROTECT_WITH_INDEX(cols = convertNegAndZeroIdx(cols, ScalarInteger(l), ScalarLogical(FALSE)), &ipx);
  int ncol = LENGTH(cols);
  int *pcols = INTEGER(cols);
  // if(ncol == 0 || (asLogical(fretall) && l == ncol)) return(x);
  // names
  SEXP nam = PROTECT(getAttrib(x, R_NamesSymbol));
  // sf data frames: Need to add sf_column
  if(oxl && asLogical(checksf) && INHERITS(x, char_sf)) {
    int sfcoln = NA_INTEGER, sf_col_sel = 0;
    SEXP *pnam = STRING_PTR(nam), sfcol = asChar(getAttrib(x, sym_sf_column));
    for(int i = l; i--; ) {
      if(pnam[i] == sfcol) {
        sfcoln = i + 1;
        break;
      }
    }
    if(sfcoln == NA_INTEGER) error("sf data frame has no attribute 'sf_column'");
    for(int i = ncol; i--; ) {
      if(pcols[i] == sfcoln) {
        sf_col_sel = 1;
        break;
      }
    }
    if(sf_col_sel == 0) {
      REPROTECT(cols = extendIntVec(cols, ncol, sfcoln), ipx);
      ++ncol;
      pcols = INTEGER(cols);
    }
  }
  SEXP ans = PROTECT(allocVector(VECSXP, ncol));
  SEXP *px = SEXPPTR(x), *pans = SEXPPTR(ans);
  for(int i = 0; i != ncol; ++i) {
    pans[i] = px[pcols[i]-1]; // SET_VECTOR_ELT(ans, i, VECTOR_ELT(x, pcols[i]-1));
  }
  if(!isNull(nam)) {
    SEXP tmp = PROTECT(allocVector(STRSXP, ncol));
    setAttrib(ans, R_NamesSymbol, tmp);
    subsetVectorRaw(tmp, nam, cols, /*anyNA=*/false);
    ++nprotect;
  }
  copyMostAttrib(x, ans); // includes row.names and class...
  // clear any index that was copied over by copyMostAttrib(), e.g. #1760 and #1734 (test 1678)
  // setAttrib(ans, sym_index, R_NilValue); -> deletes "index" attribute of pdata.frame -> don't use!!

  if(oxl && INHERITS(x, char_datatable)) {
    setAttrib(ans, sym_datatable_locked, R_NilValue);
    // int n = asInteger(GetOption1(sym_collapse_DT_alloccol));
    // UNPROTECT(nprotect); // This needs to be here !! (asInteger and GetOption1 are allocating functions)
    SEXP res = shallow(ans, R_NilValue, ncol + 100); // n // 1024 is data.table default..
    UNPROTECT(nprotect);
    return res;
    // setselfref(ans); // done by shallow
  }
  UNPROTECT(nprotect);
  return ans;
}

/*
  * subsetDT - Subsets a data.table
* NOTE:
*   1) 'rows' and 'cols' are 1-based, passed from R level
*   2) Originally for subsetting vectors in fcast and now the beginnings of [.data.table ported to C
*   3) Immediate need is for R 3.1 as lglVec[1] now returns R's global TRUE and we don't want := to change that global [think 1 row data.tables]
*   4) Could do it other ways but may as well go to C now as we were going to do that anyway
*/

SEXP subsetDT(SEXP x, SEXP rows, SEXP cols, SEXP checkrows) { // , SEXP fastret
    int nprotect=0, oxl = OBJECT(x) != 0;
    if (!isNewList(x)) error("Internal error. Argument 'x' to CsubsetDT is type '%s' not 'list'", type2char(TYPEOF(rows))); // # nocov
      if (!length(x)) return x;  // return empty list

    if (!isInteger(cols)) error("Internal error. Argument 'cols' to Csubset is type '%s' not 'integer'", type2char(TYPEOF(cols))); // # nocov
    int ncol = LENGTH(cols), l = LENGTH(x), *pcols = INTEGER(cols);
    for (int i = 0; i != ncol; ++i) {
      if (pcols[i] < 1 || pcols[i] > l) error("Item %d of 'cols' is %d which is outside 1-based range [1,ncol(x)=%d]", i+1, pcols[i], l);
    }

    const int nrow = length(VECTOR_ELT(x, pcols[0]-1)); // Allows checking just subsetted columns for right length
    // if fast return, return data.table if all rows selected through positive indices...
    // if(asLogical(fastret) && nrow == LENGTH(rows) && INTEGER(rows)[0] > 0) {
    //  if(LENGTH(cols) == length(x)) return x;
    //  return subsetCols(x, cols);
    // }
    // check index once up front for 0 or NA, for branchless subsetVectorRaw which is repeated for each column
    bool anyNA=false; // , orderedSubset=true;   // true for when rows==null (meaning all rows)
    if (asLogical(checkrows) && !isNull(rows) && check_idx(rows, nrow, &anyNA)!=NULL) { // , &orderedSubset
      SEXP max = PROTECT(ScalarInteger(nrow)); nprotect++;
      rows = PROTECT(convertNegAndZeroIdx(rows, max, ScalarLogical(TRUE))); nprotect++;
      const char *err = check_idx(rows, nrow, &anyNA); // , &orderedSubset
      if (err!=NULL) error(err);
    }

      // Adding sf geometry column if not already selected...
      if(oxl && INHERITS(x, char_sf)) {
        int sfcoln = NA_INTEGER, sf_col_sel = 0;
        SEXP nam = PROTECT(getAttrib(x, R_NamesSymbol));
        SEXP *pnam = STRING_PTR(nam), sfcol = asChar(getAttrib(x, sym_sf_column));
        for(int i = l; i--; ) {
          if(pnam[i] == sfcol) {
            sfcoln = i + 1;
            break;
          }
        }
        UNPROTECT(1);
        if(sfcoln == NA_INTEGER) error("sf data frame has no attribute 'sf_column'");
        for(int i = ncol; i--; ) {
          if(pcols[i] == sfcoln) {
            sf_col_sel = 1;
            break;
          }
        }
        if(sf_col_sel == 0) {
          cols = PROTECT(extendIntVec(cols, LENGTH(cols), sfcoln));
          ++ncol; ++nprotect;
          pcols = INTEGER(cols);
        }
      }


    // int overAlloc = 1024; // checkOverAlloc(GetOption(install("datatable.alloccol"), R_NilValue));
    SEXP ans = PROTECT(allocVector(VECSXP, ncol)); nprotect++; // +overAlloc  // doing alloc.col directly here; eventually alloc.col can be deprecated.

    // user-defined and superclass attributes get copied as from v1.12.0
    copyMostAttrib(x, ans);
    // most means all except R_NamesSymbol, R_DimSymbol and R_DimNamesSymbol
    // includes row.names (oddly, given other dims aren't) and "sorted" dealt with below
  // class is also copied here which retains superclass name in class vector as has been the case for many years; e.g. tests 1228.* for #5296

  // This is because overalloc.. creating columns by reference stuff..
  // SET_TRUELENGTH(ans, LENGTH(ans));
  // SETLENGTH(ans, LENGTH(cols));
  int ansn;
  SEXP *px = SEXPPTR(x), *pans = SEXPPTR(ans);
  if (isNull(rows)) {
    ansn = nrow;
    for (int i = 0; i != ncol; ++i) {
      SEXP thisCol = px[pcols[i]-1];
      checkCol(thisCol, pcols[i], nrow, x);
      pans[i] = thisCol; // copyAsPlain(thisCol) -> No deep copy
      // materialize the column subset as we have always done for now, until REFCNT is on by default in R (TODO)
    }
  } else {
    ansn = LENGTH(rows);  // has been checked not to contain zeros or negatives, so this length is the length of result
    for (int i = 0; i != ncol; ++i) {
      SEXP source = px[pcols[i]-1];
      checkCol(source, pcols[i], nrow, x);
      SEXP target;
      SET_VECTOR_ELT(ans, i, target = allocVector(TYPEOF(source), ansn));
      copyMostAttrib(source, target);
      subsetVectorRaw(target, source, rows, anyNA);  // parallel within column
    }
  }

  SEXP colnam = getAttrib(x, R_NamesSymbol);
  if(TYPEOF(colnam) == STRSXP) {
    PROTECT(colnam);
    SEXP tmp = PROTECT(allocVector(STRSXP, ncol)); nprotect++;
    // SET_TRUELENGTH(tmp, LENGTH(tmp));
    // SETLENGTH(tmp, LENGTH(cols));
    setAttrib(ans, R_NamesSymbol, tmp);
    subsetVectorRaw(tmp, colnam, cols, /*anyNA=*/false);
    UNPROTECT(1);
  }

  if(oxl) {
    SEXP tmp = PROTECT(allocVector(INTSXP, 2)); nprotect++;
    INTEGER(tmp)[0] = NA_INTEGER;
    INTEGER(tmp)[1] = -ansn;
    setAttrib(ans, R_RowNamesSymbol, tmp);  // The contents of tmp must be set before being passed to setAttrib(). setAttrib looks at tmp value and copies it in the case of R_RowNamesSymbol. Caused hard to track bug around 28 Sep 2014.
    // clear any index that was copied over by copyMostAttrib() above, e.g. #1760 and #1734 (test 1678)
    setAttrib(ans, sym_index, R_NilValue); // also ok for pdata.frame (can't use on subsetted or ordered data frame)
    setAttrib(ans, sym_index_df, R_NilValue);
  }

  if(oxl && INHERITS(x, char_datatable)) {
    setAttrib(ans, sym_sorted, R_NilValue);
    setAttrib(ans, sym_datatable_locked, R_NilValue);
    // int n = asInteger(GetOption1(sym_collapse_DT_alloccol));
    SEXP res = shallow(ans, R_NilValue, ncol + 100); // n // 1024 is data.table default..
    UNPROTECT(nprotect); // This needs to be here !! (asInteger and GetOption1 are allocating functions)
    return res;
    // setselfref(ans); // done by shallow
  }
  UNPROTECT(nprotect);
  return ans;
}

SEXP subsetVector(SEXP x, SEXP idx, SEXP checkidx) { // idx is 1-based passed from R level
  bool anyNA = false; //, orderedSubset=false;
  int nprotect=0;
  if (isNull(x))
    error("Internal error: NULL can not be subset. It is invalid for a data.table to contain a NULL column.");      // # nocov
  if (asLogical(checkidx) && check_idx(idx, length(x), &anyNA) != NULL) // , &orderedSubset
    error("Internal error: CsubsetVector is internal-use-only but has received negatives, zeros or out-of-range");  // # nocov
  SEXP ans = PROTECT(allocVector(TYPEOF(x), length(idx))); nprotect++;
  copyMostAttrib(x, ans);
  subsetVectorRaw(ans, x, idx, anyNA);
  UNPROTECT(nprotect);
  return ans;
}
