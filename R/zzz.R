
do_collapse_mask <- function(clpns, mask) {
  if(!is.character(mask)) stop("Option collapse_mask needs to be character typed")
  # This ensures that you can pass functions with or without f- prefix to the option
  mask_ffunl <- mask %!in% c("all", "helper", "manip", "special", "fast-fun", "fast-stat-fun", "fast-trfm-fun", "n", "qtab", "qtable", "table", "%in%")
  if(any(mask_ffunl)) {
    mask_ffun <- mask[mask_ffunl]
    has_f_prefix <- startsWith(mask_ffun, "f")
    if(!all(has_f_prefix)) {
      mask_ffun[!has_f_prefix] <- paste0("f", mask_ffun[!has_f_prefix])
      mask[mask_ffunl] <- mask_ffun
    }
  }
  # This now does the preprocessing (interpreting keywords and changing internal optimization flags as required)
  if(any(mask == "all")) mask <- c("helper", "manip", "special", "fast-fun", if(length(mask) > 1L) mask[mask != "all"] else NULL)
  manipfun <- c("fsubset", "fslice", "fslicev", "ftransform", "ftransform<-", "ftransformv", "fcompute", "fcomputev", "fselect", "fselect<-", "fgroup_by", "fgroup_vars", "fungroup", "fsummarise", "fsummarize", "fmutate", "frename", "findex_by", "findex")
  helperfun <- c("fdroplevels", "finteraction", "fnlevels", "fmatch", "funique", "fnunique", "fduplicated", "fcount", "fcountv", "fquantile", "frange", "fdist", "fnrow", "fncol") # , "fdim": Problem of infinite recursion...
  specialfun <- c("n", "table", "%in%")
  if(any(mask == "helper")) mask <- unique.default(c(helperfun, mask[mask != "helper"]))
  if(any(mask == "manip")) mask <- unique.default(c(manipfun, mask[mask != "manip"]))
  if(any(mask == "special")) mask <- unique.default(c(specialfun, mask[mask != "special"]))
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
  unmask_special <- NULL
  # Special Cases / Functions
  if(any(mask == "n")) {
    unmask_special <- "n"
    mask <- mask[mask != "n"]
    if(is.null(clpns[["n"]])) assign("n", clpns[["n_internal"]], envir = clpns)
    assign(".FAST_STAT_FUN_EXT", c(.FAST_STAT_FUN_EXT, "n"), envir = clpns)
    assign(".FAST_STAT_FUN_POLD", c(.FAST_STAT_FUN_POLD, "n"), envir = clpns)
    assign(".FAST_FUN_MOPS", c(.FAST_FUN_MOPS, "n"), envir = clpns)
  }
  if(any(mask %in% c("qtab", "qtable", "table"))) {
    if(is.null(clpns[["table"]])) assign("table", clpns[["qtab"]], envir = clpns)
    unmask_special <- c(unmask_special, "table")
    mask <- mask[!mask %in% c("qtab", "qtable", "table")]
  }
  if(any(mask == "%in%")) {
    if(is.null(clpns[["%in%"]])) assign("%in%", clpns[["%fin%"]], envir = clpns)
    unmask_special <- c(unmask_special, "%in%")
    mask <- mask[mask != "%in%"]
  }
  if(!all(m <- mask %in% names(clpns))) stop("Unsupported functions supplied to option 'collapse_mask': ", paste(mask[!m], collapse = ", "))
  if(!all(m <- startsWith(mask, "f"))) stop("All functions to me masked must start with 'f', except for 'n' and 'qtab'/'table'. You supplied: ", paste(mask[!m], collapse = ", "))
  # This now creates the additional functions (does the masking)
  unmask <- substr(mask, 2L, 100L)
  unmask_ind <- unmask %!iin% names(clpns) # Important: cannot change locked bindings in namespace!
  for(i in unmask_ind) assign(unmask[i], clpns[[mask[i]]], envir = clpns)
  # Internals of namespaceExport(clpns, c(unmask, unmask_special)):
  export_names <- c(unmask, unmask_special)
  names(export_names) <- export_names
  list2env(as.list(export_names), .getNamespaceInfo(clpns, "exports"))
}

do_collapse_remove_core <- function(clpns, rmfun, exports = TRUE, namespace = TRUE) { # exports = FALSE in .onLoad, because exports not defined yet
  if(exports) {
    clpns_exports <- .getNamespaceInfo(clpns, "exports")
    rmfun <- rmfun[rmfun %in% names(clpns_exports)] # ckmatch(rmfun, names(clpns_exports), e = "Unknown functions to be removed:")
  }
  if(any(tmp <- .FAST_STAT_FUN_EXT %in% rmfun)) assign(".FAST_STAT_FUN_EXT", .FAST_STAT_FUN_EXT[!tmp], envir = clpns)
  if(any(tmp <- .FAST_STAT_FUN_POLD %in% rmfun)) assign(".FAST_STAT_FUN_POLD", .FAST_STAT_FUN_POLD[!tmp], envir = clpns)
  if(any(tmp <- .FAST_FUN_MOPS %in% rmfun)) assign(".FAST_FUN_MOPS", .FAST_FUN_MOPS[!tmp], envir = clpns)
  if(exports) remove(list = rmfun, envir = clpns_exports)
  if(namespace) {
    assign(".COLLAPSE_ALL_EXPORTS", .COLLAPSE_ALL_EXPORTS[match(.COLLAPSE_ALL_EXPORTS, rmfun, 0L) == 0L], envir = clpns)
    remove(list = rmfun, envir = clpns)
  }
}

do_collapse_remove <- function(clpns, rmfun, ...) {
  kwd <- c("shorthand", "operator", "infix", "old") %in% rmfun
  if(kwd[1L]) rmfun <- c(rmfun[rmfun != "shorthand"], .SHORTHANDS)
  if(kwd[2L]) rmfun <- c(rmfun[rmfun != "operator"], .OPERATOR_FUN)
  if(kwd[3L]) rmfun <- c(rmfun[rmfun != "infix"], c(.COLLAPSE_ALL[startsWith(.COLLAPSE_ALL, "%")], if(any(c("%in%", "special") %in% .op[["mask"]])) "%in%"))
  if(kwd[4L]) rmfun <- c(rmfun[rmfun != "old"], .COLLAPSE_OLD)
  do_collapse_remove_core(clpns, unique.default(rmfun), ...)
}

# Used in set_collapse(), defined in global_macros.R
do_collapse_unmask <- function(clpns) {
  nam <- getNamespaceExports(clpns)
  ffuns <- nam[startsWith(nam, "f")]
  rmfun <- nam[nam %in% substr(ffuns, 2L, 100L)]
  if(any(ntab <- c("n", "table", "%in%") %in% nam)) rmfun <- c(rmfun, c("n", "table", "%in%")[ntab])
  do_collapse_remove_core(clpns, rmfun)
}

do_collapse_restore_exports <- function(clpns) {
  clpns_exports <- .getNamespaceInfo(clpns, "exports")
  missing <- fsetdiff(.COLLAPSE_ALL_EXPORTS, names(clpns_exports))
  if(length(missing)) {
  names(missing) <- missing
  list2env(as.list(missing), clpns_exports) # = namespaceExport(clpns, missing)
  }
}


.onLoad <- function(libname, pkgname) {

  res <- .Call(C_collapse_init, "init.success")
  if(!is.character(res) || res != "init.success") stop("collapse not successfully loaded!")

  # https://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r
  # https://stackoverflow.com/questions/49056642/r-how-to-make-variable-available-to-namespace-at-loading-time?noredirect=1&lq=1
  clpns <- parent.env(environment())
  assign(".collapse_env", new.env(), envir = clpns)
  .op <- new.env()
  .op$nthreads <- if(is.null(getOption("collapse_nthreads"))) 1L else as.integer(getOption("collapse_nthreads"))
  .op$na.rm <- if(is.null(getOption("collapse_na_rm")) && is.null(getOption("collapse_na.rm"))) TRUE else if(length(getOption("collapse_na_rm")))
                  as.logical(getOption("collapse_na_rm")) else as.logical(getOption("collapse_na.rm"))
  .op$sort <- if(is.null(getOption("collapse_sort"))) TRUE else as.logical(getOption("collapse_sort"))
  .op$stable.algo <- if(is.null(getOption("collapse_stable_algo"))) TRUE else as.logical(getOption("collapse_stable_algo"))
  .op$mask <- if(is.null(getOption("collapse_mask"))) NULL else getOption("collapse_mask")
  .op$remove <- if(is.null(getOption("collapse_remove"))) NULL else getOption("collapse_remove")
  .op$stub <- if(is.null(getOption("collapse_stub"))) TRUE else as.logical(getOption("collapse_stub"))
  .op$verbose <- if(is.null(getOption("collapse_verbose"))) 1L else as.integer(getOption("collapse_verbose"))
  .op$digits <- if(is.null(getOption("collapse_digits"))) 2L else as.integer(getOption("collapse_digits"))
  assign(".op", .op, envir = clpns)

  # TODO: option to save .collapse config file in install directory?? -> Nah, .RProfile is better...
  mask <- .op$mask

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
        .op$mask <- mask <- getOption("collapse_mask")
      }
    }
  }

  if(length(mask) && is.character(mask)) do_collapse_mask(clpns, mask)
  if(length(.op$remove) && is.character(.op$remove)) do_collapse_remove(clpns, .op$remove, exports = FALSE)

  if(isTRUE(getOption("collapse_export_F"))) namespaceExport(clpns, "F")
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
