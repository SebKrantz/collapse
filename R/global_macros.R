
# Global Options
set_collapse <- function(...) {
  opts <- if(...length() == 1L && is.list(..1)) ..1 else list(...)
  op_old <- as.list(.op)
  nam <- names(opts)
  if(any(nam %!in% c("nthreads", "na.rm"))) stop("Currently only supports options 'nthreads' and 'na.rm'")
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

# Get updated .COLLAPSE_ALL:
# ".default$|.matrix$|.data.frame$"
# v <- grep("\\.|N|HD", objects("package:collapse"), invert = TRUE, value = TRUE) # getNamespaceExports("collapse")
# # grep("N", objects("package:collapse"), value = TRUE)
# v <- c(v, "GRPN", "GRPid", "HDB", "HDW", "allNA", "whichNA", "replace_NA")
# TODO: also remove Date_vars...
# cat(v, sep = '", "')

# all package objects..
# allobj <- ls(getNamespace("collapse"), all.names=TRUE)

.COLLAPSE_ALL <- sort(unique(c("%-=%", "%!=%", "%!in%", "%*=%", "%/=%", "%+=%", "%=%", "%==%", "%c-%", "%c*%",
                               "%c/%", "%c+%", "%cr%", "%r-%", "%r*%", "%r/%", "%r+%", "%rr%", "add_stub",
                               "add_vars", "add_vars<-", "all_funs", "all_identical", "all_obj_equal", "alloc",
                               "allv", "anyv", "as_character_factor", "as_factor_GRP", "as_factor_qG", "as_numeric_factor",
                               "atomic_elem", "atomic_elem<-", "av", "av<-", "B", "BY", "cat_vars", "cat_vars<-", "char_vars",
                               "char_vars<-", "cinv", "ckmatch", "collap", "collapg", "collapv", "colorder", "colorderv",
                               "copyAttrib", "copyMostAttrib", "copyv", "D", "dapply", "date_vars", "date_vars<-",
                               "descr", "Dlog", "fact_vars", "fact_vars<-", "fbetween", "fcompute", "fcomputev",
                               "fcount", "fcountv", "fcumsum", "fdiff", "fdim", "fdist", "fdroplevels", "ffirst", "fFtest",
                               "fgroup_by", "fgroup_vars", "fgrowth", "fhdbetween", "fhdwithin", "findex", "findex_by",
                               "finteraction", "flag", "flast", "flm", "fmax", "fmean", "fmedian", "fmin", "fmode", "fmutate",
                               "fncol", "fndistinct", "fnlevels", "fnobs", "fnrow", "fnth", "fnunique", "fprod", "fquantile",
                               "frange", "frename", "fscale", "fsd", "fselect", "fselect<-", "fsubset", "fsum", "fsummarise",
                               "fsummarize", "ftransform", "ftransform<-", "ftransformv", "fungroup", "funique", "fvar",
                               "fwithin", "G", "gby", "get_elem", "get_vars", "get_vars<-", "GGDC10S", "greorder", "group",
                               "groupid", "GRP", "GRPid", "GRPnames", "gsplit", "gv", "gv<-", "gvr", "gvr<-", "has_elem", "iby",
                               "irreg_elem", "is_categorical", "is_date", "is_GRP", "is_irregular", "is_qG", "is_unlistable", "itn",
                               "ix", "L", "ldepth", "list_elem", "list_elem<-", "logi_vars", "logi_vars<-", "massign", "mctl", "missing_cases",
                               "mrtl", "mtt", "na_insert", "na_omit", "na_rm", "namlab", "num_vars", "num_vars<-", "nv", "nv<-", "pad", "psacf",
                               "psccf", "psmat", "pspacf", "pwcor", "pwcov", "pwnobs", "qDF", "qDT", "qF", "qG", "qM", "qsu", "qtab", "qtable",
                               "qTBL", "radixorder", "radixorderv", "rapply2d", "recode_char", "recode_num", "reg_elem", "reindex", "relabel",
                               "replace_Inf", "replace_outliers", "rm_stub", "rnm", "roworder", "roworderv", "rsplit", "sbt", "seq_col", "seq_row",
                               "seqid", "setattrib", "setAttrib", "setColnames", "setDimnames", "setLabels", "setop", "setrelabel", "setrename",
                               "setRownames", "settfm", "settfmv", "setTRA", "settransform", "settransformv", "setv", "slt", "slt<-", "smr", "ss",
                               "STD", "t_list", "tfm", "tfm<-", "tfmv", "timeid", "to_plm", "TRA", "unattrib", "unindex", "unlist2d", "varying",
                               "vclasses", "vgcd", "vlabels", "vlabels<-", "vlengths", "vtypes", "W", "whichv", "wlddev", "GRPN", "GRPid", "HDB",
                               "HDW", "allNA", "whichNA", "replace_NA")))

.COLLAPSE_GENERIC   <-   sort(unique(c("B","BY","D","Dlog","F","fsubset","fbetween","fdiff","ffirst","fgrowth","fhdbetween",
                           "fhdwithin","flag","flast","fmax","fmean","fmedian","fnth","fmin","fmode","varying",
                           "fndistinct","fnobs","fprod","fscale","fsd","fsum","fcumsum","fvar","fwithin","funique",
                           "G","GRP","HDB","HDW","L","psacf","psccf","psmat","pspacf","qsu", "rsplit","fdroplevels",
                           "STD","TRA","W")))

.COLLAPSE_DATA <- c("GGDC10S", "wlddev")

.FAST_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
               "fmin","fmax","fnth","ffirst","flast","fnobs","fndistinct",
               "fcumsum","fscale","fbetween","fwithin","fhdbetween","fhdwithin",
               "flag","fdiff","fgrowth")

.FAST_STAT_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
                    "fmin","fmax","fnth","ffirst","flast","fnobs","fndistinct")

.OPERATOR_FUN <- c("STD","B","W","HDB","HDW","L","F","D","Dlog","G")

.FAST_STAT_FUN_POLD <- c(.FAST_STAT_FUN, "fNobs","fNdistinct", "GRPN", "GRPid") # "n"

.FAST_FUN_MOPS <- c(.FAST_STAT_FUN_POLD, "fcumsum","fscale","fbetween","fwithin",
                    "flag","fdiff","fgrowth","STD","B","W","L","F","D","Dlog","G")

.FAST_STAT_FUN_EXT <- c(.FAST_STAT_FUN_POLD, paste0(setdiff(.FAST_STAT_FUN_POLD, c("GRPN", "GRPid")), "_uw")) # "n"



