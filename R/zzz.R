
.onLoad <- function(libname, pkgname) {

  res <- .Call(C_collapsedev2_init, "init.success")
  if(!is.character(res) || res != "init.success") stop("collapse not succesfully loaded!")

  # https://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r
  # https://stackoverflow.com/questions/49056642/r-how-to-make-variable-available-to-namespace-at-loading-time?noredirect=1&lq=1
  .collapse_env <- new.env()

  suppressMessages({

  .collapse_env$lfe_demeanlist <-
         if(requireNamespace("lfe", quietly = TRUE)) # lfe::demeanlist else NULL
             get0("demeanlist", envir = getNamespace("lfe")) else NULL

  .collapse_env$RcppArmadillo_fastLm <-
         if(requireNamespace("RcppArmadillo", quietly = TRUE)) # RcppArmadillo::fastLmPure else NULL
           get0("fastLmPure", envir = getNamespace("RcppArmadillo")) else NULL # _RcppArmadillo_fastLm_impl

  .collapse_env$RcppEigen_fastLm <-
         if(requireNamespace("RcppEigen", quietly = TRUE)) # RcppEigen::fastLmPure else NULL
           get0("fastLmPure", envir = getNamespace("RcppEigen")) else NULL # RcppEigen_fastLm_Impl

  })

  assign(".collapse_env", .collapse_env, envir = parent.env(environment()))

  # Old solution: does not dynamically update, would have to re-install collapse after installing these packages
  # assign(".lfe_demeanlist",
  #        if(requireNamespace("lfe", quietly = TRUE))
  #        get0("demeanlist", envir = getNamespace("lfe")) else NULL, envir = parent.env(environment()))
  # assign(".RcppArmadillo_fastLm",
  #        if(requireNamespace("RcppArmadillo", quietly = TRUE))
  #        get0("_RcppArmadillo_fastLm_impl", envir = getNamespace("RcppArmadillo")) else NULL, envir = parent.env(environment()))
  # assign(".RcppEigen_fastLm",
  #        if(requireNamespace("RcppEigen", quietly = TRUE))
  #        get0("RcppEigen_fastLm_Impl", envir = getNamespace("RcppEigen")) else NULL, envir = parent.env(environment()))

  options("collapse_unused_arg_action" = "warning") # error, warning, message or none

  invisible(res)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("collapsedev2 ",packageVersion("collapsedev2"),", see ?`collapse-package` or ?`collapse-documentation`\nNote: stats::D  ->  D.expression, D.call, D.name"))
}

.onUnload <- function (libpath) {
  library.dynam.unload("collapsedev2", libpath)
}

# Note: To create local dev version of package change package name in DESCRIPTION, NAMESPACE, this file,
# replace all instances of `_collapse_` in source files, and also rename `R_init_collapse` in ExportSymbols.cpp.

release_questions <- function() {
  c(
    "Have you updated the version number in DESCRIPTION, NEWS.md, NEWS.Rd, cran.comments and .onAttach?",
    "Updated Readme?",
    "Spell check ?",
    "built vignettes properly with Sys.setenv(NCRAN = TRUE)?",
    "Have you updated all help files with code changes, even if it's only documenting arguments or links?",
    "updated collapse-package.Rd and collapse-documentation.Rd?",
    "All function in global_macros.R?",
    "checked all depreciated functions and arguments?",
    "any changess to arguments or order of arguments in key functions (GRP etc.). Does everything work?"
  )
}
