
.onLoad <- function(libname, pkgname) {

  res <- .Call(C_collapse_init, "init.success")
  if(!is.character(res) || res != "init.success") stop("collapse not succesfully loaded!")

  # https://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r
  # https://stackoverflow.com/questions/49056642/r-how-to-make-variable-available-to-namespace-at-loading-time?noredirect=1&lq=1
  .collapse_env <- new.env()

  # This slows down th eloading of collapse too much. Therefore we load those when needed.
  # suppressMessages({
  #
  # .collapse_env$fixest_demean <-
  #       if(requireNamespace("fixest", quietly = TRUE)) # fixest::demean else NULL
  #          get0("demean", envir = getNamespace("fixest")) else NULL
  #
  # .collapse_env$weights_wtd.cors <-
  #       if(requireNamespace("weights", quietly = TRUE)) # weights::wtd.cors else NULL
  #          get0("wtd.cors", envir = getNamespace("weights")) else NULL
  #
  # .collapse_env$RcppArmadillo_fastLm <-
  #        if(requireNamespace("RcppArmadillo", quietly = TRUE)) # RcppArmadillo::fastLmPure else NULL
  #          get0("fastLmPure", envir = getNamespace("RcppArmadillo")) else NULL # _RcppArmadillo_fastLm_impl
  #
  # .collapse_env$RcppEigen_fastLm <-
  #        if(requireNamespace("RcppEigen", quietly = TRUE)) # RcppEigen::fastLmPure else NULL
  #          get0("fastLmPure", envir = getNamespace("RcppEigen")) else NULL # RcppEigen_fastLm_Impl
  #
  # })

  clpns <- parent.env(environment())
  assign(".collapse_env", .collapse_env, envir = clpns)

  # Old solution: does not dynamically update, would have to re-install collapse after installing these packages
  # assign(".RcppArmadillo_fastLm",
  #        if(requireNamespace("RcppArmadillo", quietly = TRUE))
  #        get0("_RcppArmadillo_fastLm_impl", envir = getNamespace("RcppArmadillo")) else NULL, envir = parent.env(environment()))

  if(length(mask <- getOption("collapse_mask")) && is.character(mask)) {
    # if(!is.character(mask)) stop("Option collapse_mask needs to be character typed")
    if(any(mask == "all")) mask <- c("helper", "manip", "fast-fun", if(length(mask) > 1L) mask[mask != "all"] else NULL)
    manipfun <- c("fsubset", "ftransform", "ftransform<-", "ftransformv", "fcompute", "fcomputev", "fselect", "fselect<-", "fgroup_by", "fgroup_vars", "fungroup", "fsummarise", "fmutate", "frename")
    helperfun <- c("fdroplevels", "finteraction", "fnlevels", "funique", "fnrow", "fncol") # , "fdim": Problem of infinite recursion...
    if(any(mask == "helper")) mask <- unique.default(c(helperfun, mask[mask != "helper"]))
    if(any(mask == "manip")) mask <- unique.default(c(manipfun, mask[mask != "manip"]))
    if(any(mask == "fast-fun")) {
      mask <- unique.default(c(.FAST_FUN, mask[mask != "fast-fun"]))
      fsfnonold <- .FAST_STAT_FUN_EXT[!startsWith(.FAST_STAT_FUN_EXT, "fN")]
      assign(".FAST_STAT_FUN_EXT", c(.FAST_STAT_FUN_EXT, substr(fsfnonold, 2L, 100L)), envir = clpns)
      assign(".FAST_STAT_FUN_POLD", c(.FAST_STAT_FUN_POLD, substr(.FAST_STAT_FUN, 2L, 100L)), envir = clpns)
      ffnops <-  setdiff(.FAST_FUN_MOPS, .OPERATOR_FUN)
      assign(".FAST_FUN_MOPS", c(.FAST_FUN_MOPS, substr(ffnops, 2L, 100L)), envir = clpns)
    } else {
      if(any(mask == "fast-stat-fun")) {
        mask <- unique.default(c(.FAST_STAT_FUN, mask[mask != "fast-stat-fun"]))
        fsfnonold <- .FAST_STAT_FUN_EXT[!startsWith(.FAST_STAT_FUN_EXT, "fN")]
        assign(".FAST_STAT_FUN_EXT", c(.FAST_STAT_FUN_EXT, substr(fsfnonold, 2L, 100L)), envir = clpns)
        assign(".FAST_STAT_FUN_POLD", c(.FAST_STAT_FUN_POLD, substr(.FAST_STAT_FUN, 2L, 100L)), envir = clpns)
        assign(".FAST_FUN_MOPS", c(.FAST_FUN_MOPS, substr(.FAST_STAT_FUN, 2L, 100L)), envir = clpns)
      }
      if(any(mask == "fast-trfm-fun")) {
        ftf <- fsetdiff(.FAST_FUN, .FAST_STAT_FUN)
        mask <- unique.default(c(ftf, mask[mask != "fast-trfm-fun"]))
        assign(".FAST_FUN_MOPS", c(.FAST_FUN_MOPS, substr(fsetdiff(ftf, c("fhdbetween", "fhdwithin")), 2L, 100L)), envir = clpns)
      }
    }
    if(!all(m <- mask %in% names(clpns))) stop("Unknown collapse functions supplied to option 'collapse_mask': ", paste(mask[!m], collapse = ", "))
    if(!all(m <- startsWith(mask, "f"))) stop("All functions to me masked must start with 'f'. You supplied: ", paste(mask[!m], collapse = ", "))
    unmask <- substr(mask, 2L, 100L)
    for(i in seq_along(mask)) assign(unmask[i], clpns[[mask[i]]], envir = clpns)
    namespaceExport(clpns, unmask)
  }

  if(isTRUE(getOption("collapse_F_to_FALSE"))) {
    assign("F", FALSE, envir = clpns)
  }

  # Experimental collapse_remove option: doesn't work because namespace exports not defined yet.
  # if(length(crem <- getOption("collapse_remove")) && is.character(crem)) {
  #   # clpns <- getNamespace("collapse")
  #   exports <- getNamespaceInfo(clpns, "exports") # clpns[[".__NAMESPACE__."]][["exports"]] # .getNamespaceInfo(clpns, "exports")
  #   stop("length:", length(exports))
  #   remove(list = crem, envir = exports)
  #   # setNamespaceInfo(clpns, "exports", exports)
  #   # detach("package:collapse")
  #   # attachNamespace(clpns)
  #   # clpns[[".__NAMESPACE__."]][["exports"]] <- exports
  # }



  options(collapse_unused_arg_action = "warning", # error, warning, message or none
          collapse_DT_alloccol = 100L)

  invisible(res)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("collapse ",packageVersion("collapse"),", see ?`collapse-package` or ?`collapse-documentation`")) # \nNote: stats::D  ->  D.expression, D.call, D.name
}

.onUnload <- function (libpath) {
  library.dynam.unload("collapse", libpath)
}

# Note: To create local dev version of package change package name in DESCRIPTION, NAMESPACE, this file (including C_collapse_init),
# replace all instances of `_collapse_` in source files (except for _collapse_DT_alloccol`), and also rename `R_init_collapse` in ExportSymbols.cpp.
# and in vignetter / Rd files replace library(collapse)
release_questions <- function() {
  c(
    "Have you updated the version number in DESCRIPTION, NEWS.md, NEWS.Rd, cran.comments and .onAttach?",
    "Updated Readme?",
    "Spell check ?",
    "built vignettes properly with Sys.setenv(RUNBENCH = TRUE)?",
    "Have you updated all help files with code changes, even if it's only documenting arguments or links?",
    "updated collapse-package.Rd and collapse-documentation.Rd?",
    "All functions in global_macros.R?",
    "checked all depreciated functions and arguments?",
    "any changess to arguments or order of arguments in key functions (GRP etc.). Does everything work?"
  )
}
