
.onLoad <- function(libname, pkgname) {

  res <- .Call(C_collapse_init, "init.success")
  if(!is.character(res) || res != "init.success") stop("collapse not succesfully loaded!")

  # https://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r
  # https://stackoverflow.com/questions/49056642/r-how-to-make-variable-available-to-namespace-at-loading-time?noredirect=1&lq=1
  .collapse_env <- new.env()

  suppressMessages({

  # .collapse_env$lfe_demeanlist <-
  #        if(requireNamespace("lfe", quietly = TRUE)) # lfe::demeanlist else NULL
  #            get0("demeanlist", envir = getNamespace("lfe"))

  .collapse_env$fixest_demean <-
        if(requireNamespace("fixest", quietly = TRUE)) # fixest::demean else NULL
           get0("demean", envir = getNamespace("fixest")) else NULL

  .collapse_env$weights_wtd.cors <-
        if(requireNamespace("weights", quietly = TRUE)) # weights::wtd.cors else NULL
           get0("wtd.cors", envir = getNamespace("weights")) else NULL

  .collapse_env$RcppArmadillo_fastLm <-
         if(requireNamespace("RcppArmadillo", quietly = TRUE)) # RcppArmadillo::fastLmPure else NULL
           get0("fastLmPure", envir = getNamespace("RcppArmadillo")) else NULL # _RcppArmadillo_fastLm_impl

  .collapse_env$RcppEigen_fastLm <-
         if(requireNamespace("RcppEigen", quietly = TRUE)) # RcppEigen::fastLmPure else NULL
           get0("fastLmPure", envir = getNamespace("RcppEigen")) else NULL # RcppEigen_fastLm_Impl

  })

  assign(".collapse_env", .collapse_env, envir = parent.env(environment()))

  # Old solution: does not dynamically update, would have to re-install collapse after installing these packages
  # assign(".RcppArmadillo_fastLm",
  #        if(requireNamespace("RcppArmadillo", quietly = TRUE))
  #        get0("_RcppArmadillo_fastLm_impl", envir = getNamespace("RcppArmadillo")) else NULL, envir = parent.env(environment()))

  options(collapse_unused_arg_action = "warning", # error, warning, message or none
          collapse_DT_alloccol = 100L)

  invisible(res)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("collapse ",packageVersion("collapse"),", see ?`collapse-package` or ?`collapse-documentation`\nNote: stats::D  ->  D.expression, D.call, D.name"))
}

.onUnload <- function (libpath) {
  library.dynam.unload("collapse", libpath)
}

# Note: To create local dev version of package change package name in DESCRIPTION, NAMESPACE, this file (including C_collapse_init),
# replace all instances of `_collapse_` in source files, and also rename `R_init_collapse` in ExportSymbols.cpp.
# and in vignetter / Rd files replace library(collapse)
release_questions <- function() {
  c(
    "Have you updated the version number in DESCRIPTION, NEWS.md, NEWS.Rd, cran.comments and .onAttach?",
    "Updated Readme?",
    "Spell check ?",
    "built vignettes properly with Sys.setenv(RUNBENCH = TRUE)?",
    "Have you updated all help files with code changes, even if it's only documenting arguments or links?",
    "updated collapse-package.Rd and collapse-documentation.Rd?",
    "All function in global_macros.R?",
    "checked all depreciated functions and arguments?",
    "any changess to arguments or order of arguments in key functions (GRP etc.). Does everything work?"
  )
}
