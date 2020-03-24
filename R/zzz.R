# .Call(setSizes)
# .Call(initDTthreads)
# init_collapse <- function() cat(.Call(collapse_init, "Welcome to collapse! See ?collapse"))
#

# See https://github.com/tidyverse/dplyr/blob/bbcfe99e29fe737d456b0d7adc33d3c445a32d9d/R/zzz.r
.onLoad <- function(libname, pkgname) {
  .Call(C_collapse_init)
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("collapse 1.1.0, see ?`collapse-package` or ?`collapse-documentation`")
}

.onUnload <- function (libpath) {
  library.dynam.unload("collapse", libpath)
}


release_questions <- function() {
  c(
    "Have you updated the version number in DESCRIPTION, NEWS.md, NEWS.Rd, cran.comments and .onAttach?",
    "Spell check ?",
    "built vignettes properly with Sys.setenv(NCRAN = TRUE)?"
  )
}
