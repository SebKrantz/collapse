# .Call(setSizes)
# .Call(initDTthreads)
# init_collapse <- function() cat(.Call(collapse_init, "Welcome to collapse! See ?collapse"))
#

# See https://github.com/tidyverse/dplyr/blob/bbcfe99e29fe737d456b0d7adc33d3c445a32d9d/R/zzz.r
.onLoad <- function(libname, pkgname) {
  .Call(collapse_init)
  invisible()
}

.onAttach <- function(libname, pkgname) {
  cat("Welcome to collapse! See ?collapse")
}

