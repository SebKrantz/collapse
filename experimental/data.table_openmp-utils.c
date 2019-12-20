#include "data.table.h"
#ifdef _OPENMP
#include <pthread.h>
#endif
#include <errno.h>     // errno
#include <ctype.h>     // isspace

static int  DTthreads = -1;           // Never read directly hence static; use getDTthreads(). -1 so we know for sure initDTthreads() ran and set it >= 1.
static bool RestoreAfterFork = true;  // see #2885 in v1.12.0

static int getIntEnv(const char *name, int def) // do I need this ?? -> yes, initDTthreads requires it !!
{
  const char *val = getenv(name);
  if (val==NULL) return def;
  size_t nchar = strlen(val);
  if (nchar==0) return def;
  char *end;
  errno = 0;
  long int ans = strtol(val, &end, 10);  // ignores leading whitespace. If it fully consumed the string, *end=='\0' and isspace('\0')==false
  while (isspace(*end)) end++;  // ignore trailing whitespace
  if (errno || (size_t)(end-val)!=nchar || ans<1 || ans>INT_MAX) {
    warning("Ignoring invalid %s==\"%s\". Not an integer >= 1. Please remove any characters that are not a digit [0-9]. See ?data.table::setDTthreads.", name, val);
    return def;
  }
  return (int)ans;
}

static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }

void initDTthreads() { // void
  // called at package startup from init.c
  // also called by setDTthreads(threads=NULL) (default) to reread environment variables; see setDTthreads below
  // No verbosity here in this setter. Verbosity is in getDTthreads(verbose=TRUE)
  int ans = omp_get_num_procs();  // starting point is all logical CPUs. This is a hard limit; user cannot achieve more than this.
                                  // ifndef _OPENMP then myomp.h defines this to be 1
  int perc = getIntEnv("R_DATATABLE_NUM_PROCS_PERCENT", 50); // use "NUM_PROCS" to use the same name as the OpenMP function this uses
  // 50% of logical CPUs by default; half of 8 is 4 on laptop with 4 cores. Leaves plenty of room for other processes: #3395 & #3298
  if (perc<=1 || perc>100) {
    warning("Ignoring invalid R_DATATABLE_NUM_PROCS_PERCENT==%d. If used it must be an integer between 2 and 100. Default is 50. See ?setDTtheads.", perc);
    // not allowing 1 is to catch attempts to use 1 or 1.0 to represent 100%.
    perc = 50;
  }
  ans = imax(ans*perc/100, 1);
  ans = imin(ans, omp_get_thread_limit());  // honors OMP_THREAD_LIMIT when OpenMP started; e.g. CRAN sets this to 2. Often INT_MAX meaning unlimited/unset
  ans = imin(ans, omp_get_max_threads());   // honors OMP_NUM_THREADS when OpenMP started, plus reflects any omp_set_* calls made since
  ans = imax(ans, 1);                       // just in case omp_get_* returned <= 0 for any reason
  // max_threads() -vs- num_procs(): https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/302866
  ans = imin(ans, getIntEnv("R_DATATABLE_NUM_THREADS", INT_MAX));
  ans = imin(ans, getIntEnv("OMP_THREAD_LIMIT", INT_MAX));  // user might expect `Sys.setenv(OMP_THREAD_LIMIT=2);setDTthreads()` to work. Satisfy this
  ans = imin(ans, getIntEnv("OMP_NUM_THREADS", INT_MAX));   //   expectation by reading them again now. OpenMP just reads them on startup (quite reasonably)
  DTthreads = ans;
  // return ScalarInteger(ans);
}

int getDTthreads() {
  // this is the main getter used by all parallel regions; they specify num_threads(getDTthreads())
  // Therefore keep it light, simple and robust. Local static variable. initDTthreads() ensures 1 <= DTthreads <= omp_get_num_proc()
  return DTthreads;
}


/*
  Automatically drop down to 1 thread when called from parallel package (e.g. mclapply) to avoid
  deadlock when data.table is used from within parallel::mclapply; #1745 and #1727.
  GNU OpenMP seems ok with just setting DTthreads to 1 which limits the next parallel region
  if data.table is used within the fork'd proceess. This is tested by test 1705.

  From v1.12.0 we're trying again to RestoreAferFork (#2285) with optional-off due to success
  reported by Ken Run and Mark Klik in fst#110 and fst#112. We had tried that before but had
  experienced problems likely on Intel's OpenMP only (Mac).

  DO NOT call omp_set_num_threads(1) inside when_fork()!! That causes a different crash/hang on MacOS
  upon mclapply's fork even if data.table is merely loaded and neither used yet in the session nor by
  what mclapply is calling. Even when passing on CRAN's MacOS all-OK. As discovered by several MacOS
  and RStudio users here when 1.10.4-2 went to CRAN :
     https://github.com/Rdatatable/data.table/issues/2418
  v1.10.4-3 removed the call to omp_set_num_threads().
  We do not know why calling (mere) omp_set_num_threads(1) in the when_fork() would cause a problem
  on MacOS.
*/

static int pre_fork_DTthreads = 0;

void when_fork() {
  pre_fork_DTthreads = DTthreads;
  DTthreads = 1;
}

void after_fork() {
  if (RestoreAfterFork) DTthreads = pre_fork_DTthreads;
}

void avoid_openmp_hang_within_fork() {
  // Called once on loading data.table from init.c
#ifdef _OPENMP
  pthread_atfork(&when_fork, &after_fork, NULL);
#endif
}
