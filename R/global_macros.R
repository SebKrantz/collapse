
# Global Options
set_collapse <- function(...) {
  opts <- if(...length() == 1L && is.list(..1)) ..1 else list(...)
  op_old <- as.list(.op)
  nam <- names(opts)
  ckmatch(nam, c("nthreads", "na.rm", "sort", "stable.algo", "mask", "remove", "stub", "verbose", "digits"), e = "Unknown option:")
  if(length(opts$nthreads)) {
    nthreads <- as.integer(opts$nthreads)
    if(is.na(nthreads) || nthreads <= 0L) stop("nthreads needs to be a positive integer")
    .op$nthreads <- nthreads
  }
  if(length(opts$na.rm)) {
    na.rm <- as.logical(opts$na.rm)
    if(is.na(na.rm)) stop("na.rm needs to be TRUE or FALSE")
    .op$na.rm <- na.rm
  }
  if(length(opts$sort)) {
    sort <- as.logical(opts$sort)
    if(is.na(sort)) stop("sort needs to be TRUE or FALSE")
    .op$sort <- sort
  }
  if(length(opts$stable.algo)) {
    stable.algo <- as.logical(opts$stable.algo)
    if(is.na(stable.algo)) stop("stable.algo needs to be TRUE or FALSE")
    .op$stable.algo <- stable.algo
  }
  if(length(opts$stub)) {
    stub <- as.logical(opts$stub)
    if(is.na(stub)) stop("stub needs to be TRUE or FALSE")
    .op$stub <- stub
  }
  if(length(opts$verbose)) {
    verbose <- as.integer(opts$verbose)
    if(is.na(verbose) || verbose < 0L) stop("verbose needs to be a non-negative integer")
    .op$verbose <- verbose
  }
  if(length(opts$digits)) {
    digits <- as.integer(opts$digits)
    if(is.na(digits) || digits < 0L) stop("digits needs to be a non-negative integer")
    .op$digits <- digits
  }
  if(any(mrl <- c("mask", "remove") %in% nam)) { # either can be NULL
    maskl <- mrl[1L] && !identical(op_old$mask, opts$mask)
    removel <- mrl[2L] && !identical(op_old$remove, opts$remove)
    if(maskl || removel) {
      clpns <- getNamespace("collapse")
      .Call(C_unlock_collapse_namespace, clpns)
      if(!maskl) opts$mask <- op_old$mask
      # problem: option remove does not restore masked exports, e.g. when moving from remove = "between" to remove = NULL when mask = "all" (and not changing)
      if(maskl && length(op_old$mask)) do_collapse_unmask(clpns) # Fixed in do_collapse_mask(): not overriding already masked function in namespace anymore
      if(length(opts$mask)) do_collapse_mask(clpns, opts$mask)
      .op$mask <- opts$mask
      if(removel || (maskl && length(op_old$remove))) { # When changing mask setting also need to change remove again if specified
        if(!removel) opts$remove <- op_old$remove
        if(removel && length(op_old$remove)) do_collapse_restore_exports(clpns) # Also adjusted do_collapse_remove() to only remove existing funs
        if(length(opts$remove)) do_collapse_remove(clpns, opts$remove, namespace = FALSE)
        .op$remove <- opts$remove
      }
      lockEnvironment(clpns, bindings = TRUE)
      if(anyv(search(), "package:collapse")) {
        detach("package:collapse")
        suppressPackageStartupMessages(attachNamespace(clpns))
      }
    }
  }
  invisible(op_old)
}



get_collapse <- function(opts = NULL) if(is.null(opts)) as.list(.op) else if(length(opts) == 1L) .op[[opts]] else `names<-`(lapply(opts, function(x) .op[[x]]), opts)

# Global Macros


.COLLAPSE_TOPICS <- c("collapse-documentation","fast-statistical-functions","fast-grouping-ordering",
                      "fast-data-manipulation","quick-conversion","advanced-aggregation",
                      "data-transformations","time-series-panel-series","list-processing",
                      "summary-statistics","recode-replace","efficient-programming","small-helpers","collapse-options")

# .COLLAPSE_TOPICS <- c("collapse-documentation","A1-fast-statistical-functions","A2-fast-grouping-ordering",
#                       "A3-fast-data-manipulation","A4-quick-conversion","A5-advanced-aggregation",
#                       "A6-data-transformations","A7-time-series-panel-series","A8-list-processing",
#                       "A9-summary-statistics","AA1-recode-replace","AA2-efficient-programming","AA3-small-helpers")

# rd <- tools::Rd_db("collapse")
# .COLLAPSE_HELP <- unlist(lapply(rd, tools:::.Rd_get_metadata, "name"), use.names = FALSE)
# grep("^A|depreciated", unlist(lapply(rd, tools:::.Rd_get_metadata, "name"), use.names = FALSE), invert = TRUE, value = TRUE)

# # Get updated .COLLAPSE_ALL:
# # ".default$|.matrix$|.data.frame$"
# v <- grep("\\.|N|HD", objects("package:collapse"), invert = TRUE, value = TRUE) # getNamespaceExports("collapse")
# # grep("N", objects("package:collapse"), value = TRUE)
# v <- c(v, "GRPN", "GRPid", "HDB", "HDW", "allNA", "whichNA", "replace_NA")
# # TODO: also remove Date_vars...
# cat(unique(sort(v)), sep = '", "')

# all package objects..
# allobj <- ls(getNamespace("collapse"), all.names=TRUE)

# dput(setdiff(objects("package:collapse"), .COLLAPSE_DATA))
.COLLAPSE_ALL_EXPORTS <-  c("%-=%", "%!=%", "%!iin%", "%!in%", "%*=%", "%/=%", "%+=%",
                            "%=%", "%==%", "%c-%", "%c*%", "%c/%", "%c+%", "%cr%", "%iin%",
                            "%r-%", "%r*%", "%r/%", "%r+%", "%rr%", "add_stub", "add_vars",
                            "add_vars<-", "all_funs", "all_identical", "all_obj_equal", "allNA",
                            "alloc", "allv", "any_duplicated", "anyv", "as_character_factor",
                            "as_factor_GRP", "as_factor_qG", "as_numeric_factor", "as.character_factor",
                            "as.factor_GRP", "as.factor_qG", "as.numeric_factor", "atomic_elem",
                            "atomic_elem<-", "av", "av<-", "B", "BY", "BY.data.frame", "BY.default",
                            "BY.matrix", "cat_vars", "cat_vars<-", "char_vars", "char_vars<-",
                            "cinv", "ckmatch", "collap", "collapg", "collapv", "colorder",
                            "colorderv", "copyAttrib", "copyMostAttrib", "copyv", "D", "dapply",
                            "date_vars", "Date_vars", "date_vars<-", "Date_vars<-", "descr",
                            "descr.default", "Dlog", "fact_vars", "fact_vars<-", "fbetween",
                            "fbetween.data.frame", "fbetween.default", "fbetween.matrix",
                            "fcompute", "fcomputev", "fcount", "fcountv", "fcumsum", "fcumsum.data.frame",
                            "fcumsum.default", "fcumsum.matrix", "fdiff", "fdiff.data.frame",
                            "fdiff.default", "fdiff.matrix", "fdim", "fdist", "fdroplevels",
                            "fdroplevels.data.frame", "fdroplevels.factor", "fduplicated",
                            "ffirst", "ffirst.data.frame", "ffirst.default", "ffirst.matrix",
                            "fFtest", "fFtest.default", "fgroup_by", "fgroup_vars", "fgrowth",
                            "fgrowth.data.frame", "fgrowth.default", "fgrowth.matrix", "fhdbetween",
                            "fHDbetween", "fhdbetween.data.frame", "fhdbetween.default",
                            "fhdbetween.matrix", "fhdwithin", "fHDwithin", "fhdwithin.data.frame",
                            "fhdwithin.default", "fhdwithin.matrix", "findex", "findex_by",
                            "finteraction", "flag", "flag.data.frame", "flag.default", "flag.matrix",
                            "flast", "flast.data.frame", "flast.default", "flast.matrix",
                            "flm", "flm.default", "fmatch", "fmax", "fmax.data.frame", "fmax.default",
                            "fmax.matrix", "fmean", "fmean.data.frame", "fmean.default",
                            "fmean.matrix", "fmedian", "fmedian.data.frame", "fmedian.default",
                            "fmedian.matrix", "fmin", "fmin.data.frame", "fmin.default",
                            "fmin.matrix", "fmode", "fmode.data.frame", "fmode.default",
                            "fmode.matrix", "fmutate", "fncol", "fndistinct", "fNdistinct",
                            "fndistinct.data.frame", "fndistinct.default", "fndistinct.matrix",
                            "fnlevels", "fnobs", "fNobs", "fnobs.data.frame", "fnobs.default",
                            "fnobs.matrix", "fnrow", "fnth", "fnth.data.frame", "fnth.default",
                            "fnth.matrix", "fnunique", "fprod", "fprod.data.frame", "fprod.default",
                            "fprod.matrix", "fquantile", "frange", "frename", "fscale", "fscale.data.frame",
                            "fscale.default", "fscale.matrix", "fsd", "fsd.data.frame", "fsd.default",
                            "fsd.matrix", "fselect", "fselect<-", "fsubset", "fsubset.data.frame",
                            "fsubset.default", "fsubset.matrix", "fsum", "fsum.data.frame",
                            "fsum.default", "fsum.matrix", "fsummarise", "fsummarize", "ftransform",
                            "ftransform<-", "ftransformv", "fungroup", "funique", "funique.data.frame",
                            "funique.default", "fvar", "fvar.data.frame", "fvar.default",
                            "fvar.matrix", "fwithin", "fwithin.data.frame", "fwithin.default",
                            "fwithin.matrix", "G", "gby", "get_collapse", "get_elem", "get_vars",
                            "get_vars<-", "greorder", "group", "groupid", "GRP", "GRP.default",
                            "GRPid", "GRPN", "GRPnames", "gsplit", "gv", "gv<-", "gvr", "gvr<-",
                            "has_elem", "HDB", "HDW", "iby", "irreg_elem", "is_categorical",
                            "is_date", "is_GRP", "is_irregular", "is_qG", "is_unlistable",
                            "is.categorical", "is.Date", "is.GRP", "is.qG", "is.unlistable",
                            "itn", "ix", "join", "L", "ldepth", "list_elem", "list_elem<-",
                            "logi_vars", "logi_vars<-", "massign", "mctl", "missing_cases",
                            "mrtl", "mtt", "na_insert", "na_omit", "na_rm", "namlab", "num_vars",
                            "num_vars<-", "nv", "nv<-", "pad", "pivot", "plot.psmat", "print.pwcor",
                            "print.pwcov", "print.qsu", "psacf", "psacf.data.frame", "psacf.default",
                            "psccf", "psccf.default", "psmat", "psmat.data.frame", "psmat.default",
                            "pspacf", "pspacf.data.frame", "pspacf.default", "pwcor", "pwcov",
                            "pwnobs", "pwNobs", "qDF", "qDT", "qF", "qG", "qM", "qsu", "qsu.data.frame",
                            "qsu.default", "qsu.matrix", "qtab", "qtable", "qTBL", "radixorder",
                            "radixorderv", "rapply2d", "recode_char", "recode_num", "reg_elem",
                            "reindex", "relabel", "replace_inf", "replace_Inf", "replace_na",
                            "replace_NA", "replace_outliers", "rm_stub", "rnm", "rowbind",
                            "roworder", "roworderv", "rsplit", "rsplit.data.frame", "rsplit.default",
                            "rsplit.matrix", "sbt", "seq_col", "seq_row", "seqid", "set_collapse",
                            "setattrib", "setAttrib", "setColnames", "setDimnames", "setLabels",
                            "setop", "setrelabel", "setrename", "setRownames", "settfm",
                            "settfmv", "setTRA", "settransform", "settransformv", "setv",
                            "slt", "slt<-", "smr", "ss", "STD", "t_list", "tfm", "tfm<-",
                            "tfmv", "timeid", "to_plm", "TRA", "TRA.data.frame", "TRA.default",
                            "TRA.matrix", "unattrib", "unindex", "unlist2d", "varying", "varying.data.frame",
                            "varying.default", "varying.matrix", "vclasses", "vec", "vgcd",
                            "vlabels", "vlabels<-", "vlengths", "vtypes", "W", "whichNA",
                            "whichv")

.COLLAPSE_ALL <- sort(unique(c("%-=%", "%!=%", "%!iin%", "%!in%", "%*=%", "%/=%", "%+=%", "%=%", "%==%", "%c-%", "%c*%", "%c/%", "%c+%",
                               "%cr%", "%iin%", "%r-%", "%r*%", "%r/%", "%r+%", "%rr%", "add_stub", "add_vars", "add_vars<-", "all_funs",
                               "all_identical", "all_obj_equal", "allNA", "alloc", "allv", "any_duplicated", "anyv", "as_character_factor",
                               "as_factor_GRP", "as_factor_qG", "as_numeric_factor", "atomic_elem", "atomic_elem<-", "av", "av<-", "B", "BY",
                               "cat_vars", "cat_vars<-", "char_vars", "char_vars<-", "cinv", "ckmatch", "collap", "collapg", "collapv", "colorder",
                               "colorderv", "copyAttrib", "copyMostAttrib", "copyv", "D", "dapply", "date_vars", "Date_vars", "date_vars<-",
                               "Date_vars<-", "descr", "Dlog", "fact_vars", "fact_vars<-", "fbetween", "fcompute", "fcomputev", "fcount",
                               "fcountv", "fcumsum", "fdiff", "fdim", "fdist", "fdroplevels", "fduplicated", "ffirst", "fFtest", "fgroup_by",
                               "fgroup_vars", "fgrowth", "fhdbetween", "fhdwithin", "findex", "findex_by", "finteraction", "flag", "flast", "flm",
                               "fmatch", "fmax", "fmean", "fmedian", "fmin", "fmode", "fmutate", "fncol", "fndistinct", "fnlevels", "fnobs", "fnrow",
                               "fnth", "fnunique", "fprod", "fquantile", "frange", "frename", "fscale", "fsd", "fselect", "fselect<-", "fsubset", "fsum",
                               "fsummarise", "fsummarize", "ftransform", "ftransform<-", "ftransformv", "fungroup", "funique", "fvar", "fwithin", "G",
                               "gby", "get_collapse", "get_elem", "get_vars", "get_vars<-", "GGDC10S", "greorder", "group", "groupid", "GRP", "GRPid",
                               "GRPN", "GRPnames", "gsplit", "gv", "gv<-", "gvr", "gvr<-", "has_elem", "HDB", "HDW", "iby", "irreg_elem", "is_categorical",
                               "is_date", "is_GRP", "is_irregular", "is_qG", "is_unlistable", "itn", "ix", "join", "L", "ldepth", "list_elem", "list_elem<-",
                               "logi_vars", "logi_vars<-", "massign", "mctl", "missing_cases", "mrtl", "mtt", "na_insert", "na_omit", "na_rm", "namlab",
                               "num_vars", "num_vars<-", "nv", "nv<-", "pad", "pivot", "psacf", "psccf", "psmat", "pspacf", "pwcor", "pwcov", "pwnobs",
                               "qDF", "qDT", "qF", "qG", "qM", "qsu", "qtab", "qtable", "qTBL", "radixorder", "radixorderv", "rapply2d", "recode_char",
                               "recode_num", "reg_elem", "reindex", "relabel", "replace_inf", "replace_Inf", "replace_na", "replace_NA", "replace_outliers",
                               "rm_stub", "rnm", "rowbind", "roworder", "roworderv", "rsplit", "sbt", "seq_col", "seq_row", "seqid", "set_collapse",
                               "setattrib", "setAttrib", "setColnames", "setDimnames", "setLabels", "setop", "setrelabel", "setrename", "setRownames",
                               "settfm", "settfmv", "setTRA", "settransform", "settransformv", "setv", "slt", "slt<-", "smr", "ss", "STD", "t_list",
                               "tfm", "tfm<-", "tfmv", "timeid", "to_plm", "TRA", "unattrib", "unindex", "unlist2d", "varying", "vclasses", "vec", "vgcd",
                               "vlabels", "vlabels<-", "vlengths", "vtypes", "W", "whichNA", "whichv", "wlddev")))

.COLLAPSE_GENERIC   <-   sort(unique(c("B","BY","D","Dlog","F","fsubset","fbetween","fdiff","ffirst","fgrowth","fhdbetween",
                           "fhdwithin","flag","flast","fmax","fmean","fmedian","fnth","fmin","fmode","varying",
                           "fndistinct","fnobs","fprod","fscale","fsd","fsum","fcumsum","fvar","fwithin","funique",
                           "G","GRP","HDB","HDW","L","psacf","psccf","psmat","pspacf","qsu", "rsplit","fdroplevels",
                           "STD","TRA","W", "descr")))

.COLLAPSE_DATA <- c("GGDC10S", "wlddev")

.FAST_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
               "fmin","fmax","fnth","ffirst","flast","fnobs","fndistinct",
               "fcumsum","fscale","fbetween","fwithin","fhdbetween","fhdwithin",
               "flag","fdiff","fgrowth")

.FAST_STAT_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
                    "fmin","fmax","fnth","ffirst","flast","fnobs","fndistinct")

.OPERATOR_FUN <- c("STD","B","W","HDB","HDW","L","F","D","Dlog","G")

.SHORTHANDS <- c("gv", "gv<-", "av", "av<-", "nv", "nv<-", "gvr", "gvr<-", "itn", "ix",
                 "slt", "slt<-", "sbt", "gby", "iby", "mtt", "smr",
                 "tfm", "tfmv", "tfm<-", "settfm", "settfmv", "rnm")

.COLLAPSE_OLD <- c("fNobs", "fNdistinct", "pwNobs", "fHDwithin", "fHDbetween", "as.factor_GRP", "as.factor_qG",
             "is.GRP", "is.qG", "is.unlistable", "is.categorical", "is.Date", "as.numeric_factor",
              "as.character_factor", "Date_vars", "Date_vars<-", "replace_NA", "replace_Inf")

.FAST_STAT_FUN_POLD <- c(.FAST_STAT_FUN, "fNobs","fNdistinct", "GRPN", "GRPid") # "n"

.FAST_FUN_MOPS <- c(.FAST_STAT_FUN_POLD, "fcumsum","fscale","fbetween","fwithin",
                    "flag","fdiff","fgrowth","STD","B","W","L","F","D","Dlog","G")

.FAST_STAT_FUN_EXT <- c(.FAST_STAT_FUN_POLD, paste0(setdiff(.FAST_STAT_FUN_POLD, c("GRPN", "GRPid")), "_uw")) # "n"



