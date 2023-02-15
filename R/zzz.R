
.onLoad <- function(libname, pkgname) {

  res <- .Call(C_collapse_init, "init.success")
  if(!is.character(res) || res != "init.success") stop("collapse not succesfully loaded!")

  # https://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r
  # https://stackoverflow.com/questions/49056642/r-how-to-make-variable-available-to-namespace-at-loading-time?noredirect=1&lq=1
  clpns <- parent.env(environment())
  .collapse_env <- new.env()
  assign(".collapse_env", .collapse_env, envir = clpns)
  .op <- new.env()
  .op$nthreads <- if(is.null(getOption("collapse_nthreads"))) 1L else as.integer(getOption("collapse_nthreads"))
  .op$na.rm <- if(is.null(getOption("collapse_na.rm"))) TRUE else as.logical(getOption("collapse_na.rm"))
  .op$sort <- if(is.null(getOption("collapse_sort"))) TRUE else as.logical(getOption("collapse_sort"))
  assign(".op", .op, envir = clpns)

  mask <- getOption("collapse_mask")

  # This checks if a .fastverse config file is there: to make sure collapse cannot be loaded without masking in project
  if(!(length(mask) && is.character(mask))) {
    if(file.exists(".fastverse")) {
      fileConn <- file(".fastverse")
      contents <- readLines(fileConn, warn = FALSE, skipNul = TRUE)
      close(fileConn)
      contents <- trimws(contents[nzchar(contents)])
      mask <- which(startsWith(contents, "_opt_collapse_mask")) # Also works with if-clause below
      if(length(mask)) {
        if(length(mask) > 1L) stop("Multiple collapse_mask options set in .fastverse config file")
        mask <- paste0("options(", substr(contents[mask], 6L, 100000L), ")")
        eval(str2lang(mask))
        mask <- getOption("collapse_mask")
      }
    }
  }

  if(length(mask) && is.character(mask)) {
    # if(!is.character(mask)) stop("Option collapse_mask needs to be character typed")
    mask_all <- any(mask == "all")
    if(mask_all) mask <- c("helper", "manip", "fast-fun", if(length(mask) > 1L) mask[mask != "all"] else NULL)
    manipfun <- c("fsubset", "ftransform", "ftransform<-", "ftransformv", "fcompute", "fcomputev", "fselect", "fselect<-", "fgroup_by", "fgroup_vars", "fungroup", "fsummarise", "fsummarize", "fmutate", "frename", "findex_by", "findex")
    helperfun <- c("fdroplevels", "finteraction", "fnlevels", "funique", "fnunique", "fduplicated", "fcount", "fcountv", "fquantile", "frange", "fdist", "fnrow", "fncol") # , "fdim": Problem of infinite recursion...
    if(any(mask == "helper")) mask <- unique.default(c(helperfun, mask[mask != "helper"]))
    if(any(mask == "manip")) mask <- unique.default(c(manipfun, mask[mask != "manip"]))
    if(any(mask == "fast-fun")) {
      mask <- unique.default(c(.FAST_FUN, mask[mask != "fast-fun"]))
      FSF_mask <- substr(.FAST_STAT_FUN, 2L, 100L)
      assign(".FAST_STAT_FUN_EXT", c(.FAST_STAT_FUN_EXT, FSF_mask, paste0(FSF_mask, "_uw")), envir = clpns)
      assign(".FAST_STAT_FUN_POLD", c(.FAST_STAT_FUN_POLD, FSF_mask), envir = clpns)
      ffnops <-  fsetdiff(.FAST_FUN_MOPS, c(.OPERATOR_FUN, "fNobs", "fNdistinct", "GRPN", "GRPid", "n"))
      assign(".FAST_FUN_MOPS", c(.FAST_FUN_MOPS, substr(ffnops, 2L, 100L)), envir = clpns)
    } else {
      if(any(mask == "fast-stat-fun")) {
        mask <- unique.default(c(.FAST_STAT_FUN, mask[mask != "fast-stat-fun"]))
        FSF_mask <- substr(.FAST_STAT_FUN, 2L, 100L)
        assign(".FAST_STAT_FUN_EXT", c(.FAST_STAT_FUN_EXT, FSF_mask, paste0(FSF_mask, "_uw")), envir = clpns)
        assign(".FAST_STAT_FUN_POLD", c(.FAST_STAT_FUN_POLD, FSF_mask), envir = clpns)
        assign(".FAST_FUN_MOPS", c(.FAST_FUN_MOPS, FSF_mask), envir = clpns)
      }
      if(any(mask == "fast-trfm-fun")) {
        ftf <- fsetdiff(.FAST_FUN, .FAST_STAT_FUN)
        mask <- unique.default(c(ftf, mask[mask != "fast-trfm-fun"]))
        assign(".FAST_FUN_MOPS", c(.FAST_FUN_MOPS, substr(fsetdiff(ftf, c("fhdbetween", "fhdwithin")), 2L, 100L)), envir = clpns)
      }
    }
    if(!all(m <- mask %in% names(clpns))) stop("Unknown collapse functions supplied to option 'collapse_mask': ", paste(mask[!m], collapse = ", "))
    unmask_special <- NULL
    # Special Cases
    if(mask_all || any(mask == "n")) {
      unmask_special <- "n"
      mask <- mask[mask != "n"]
      assign(".FAST_STAT_FUN_EXT", c(.FAST_STAT_FUN_EXT, "n"), envir = clpns)
      assign(".FAST_STAT_FUN_POLD", c(.FAST_STAT_FUN_POLD, "n"), envir = clpns)
      assign(".FAST_FUN_MOPS", c(.FAST_FUN_MOPS, "n"), envir = clpns)
    }
    if(mask_all || any(mask %in% c("qtab", "qtable"))) {
      assign("table", clpns[["qtab"]], envir = clpns)
      unmask_special <- c(unmask_special, "table")
      mask <- mask[!mask %in% c("qtab", "qtable")]
    }
    if(!all(m <- startsWith(mask, "f"))) stop("All functions to me masked must start with 'f', except for 'n' and 'qtab'. You supplied: ", paste(mask[!m], collapse = ", "))
    unmask <- substr(mask, 2L, 100L)
    for(i in seq_along(mask)) assign(unmask[i], clpns[[mask[i]]], envir = clpns)
    namespaceExport(clpns, c(unmask, unmask_special))
  }

  if(isTRUE(getOption("collapse_export_F"))) namespaceExport(clpns, "F")

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

  if(is.null(getOption("collapse_unused_arg_action"))) options(collapse_unused_arg_action = "warning") # error, warning, message or none
  # if(is.null(getOption("collapse_DT_alloccol"))) options(collapse_DT_alloccol = 100L)

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
# and in vignettes / Rd files replace library(collapse)
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
    "any changes to arguments or order of arguments in key functions (GRP etc.). Does everything work?"
  )
}
