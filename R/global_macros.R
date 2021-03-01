
# .NA_RM <- TRUE

# global macros

.COLLAPSE_TOPICS <- c("collapse-documentation","A1-fast-statistical-functions","A2-fast-grouping-ordering",
                      "A3-fast-data-manipulation","A4-quick-conversion","A5-advanced-aggregation",
                      "A6-data-transformations","A7-time-series-panel-series","A8-list-processing",
                      "A9-summary-statistics","AA1-recode-replace","AA2-small-helpers")

# rd <- tools::Rd_db("collapse")
# .COLLAPSE_HELP <- unlist(lapply(rd, tools:::.Rd_get_metadata, "name"), use.names = FALSE)
# grep("^A|depreciated", unlist(lapply(rd, tools:::.Rd_get_metadata, "name"), use.names = FALSE), invert = TRUE, value = TRUE)

# Get updated .COLLAPSE_ALL:
# v <- grep(".default$|.matrix$|.data.frame$", objects("package:collapse"), invert = TRUE, value = TRUE)
# cat(v, sep = '", "')

# all package objects..
# allobj <- ls(getNamespace("collapse"), all.names=TRUE)


.COLLAPSE_ALL <- sort(unique(c("%!in%", "%c-%", "%c*%", "%c/%", "%c+%", "%cr%", "%r-%", "%r*%", "%r/%", "%r+%", "%rr%", "add_stub", "add_vars", "add_vars<-", "all_identical", "all_obj_equal", "allNA", "as.character_factor", "as.factor_GRP", "as.factor_qG", "as.numeric_factor", "atomic_elem", "atomic_elem<-", "av", "av<-", "B", "BY", "cat_vars", "cat_vars<-", "char_vars", "char_vars<-", "cinv", "ckmatch", "collap", "collapg", "collapv", "colorder", "colorderv", "copyAttrib", "copyMostAttrib", "D", "dapply", "Date_vars", "Date_vars<-", "descr", "Dlog", "F", "fact_vars", "fact_vars<-", "fbetween", "fcompute", "fdiff", "fdim", "fdroplevels", "fdroplevels.factor", "ffirst", "fFtest", "fgroup_by", "fgroup_vars", "fgrowth", "fHDbetween", "fHDwithin", "finteraction", "flag", "flast", "flm", "fmax", "fmean", "fmedian", "fmin", "fmode", "fncol", "fNdistinct", "fnlevels", "fNobs", "fnrow", "fnth", "fprod", "frename", "fscale", "fsd", "fselect", "fselect<-", "fsubset", "fsum", "fsummarise", "ftransform", "ftransform<-", "ftransformv", "fungroup", "funique", "fvar", "fwithin", "G", "gby", "get_elem", "get_vars", "get_vars<-", "GGDC10S", "groupid", "GRP", "GRPnames", "gv", "gv<-", "gvr", "gvr<-", "has_elem", "HDB", "HDW", "irreg_elem", "is.categorical", "is.Date", "is.GRP", "is.qG", "is.regular", "is.unlistable", "L", "ldepth", "list_elem", "list_elem<-", "logi_vars", "logi_vars<-", "mctl", "missing_cases", "mrtl", "na_insert", "na_omit", "na_rm", "namlab", "num_vars", "num_vars<-", "nv", "nv<-", "print.pwcor", "print.pwcov", "print.qsu", "psacf", "psccf", "psmat", "pspacf", "pwcor", "pwcov", "pwNobs", "qDF", "qDT", "qF", "qG", "qM", "qsu", "qTBL", "radixorder", "radixorderv", "rapply2d", "Recode", "recode_char", "recode_num", "reg_elem", "replace_Inf", "replace_NA", "replace_non_finite", "replace_outliers", "rm_stub", "roworder", "roworderv", "rsplit", "sbt", "seq_col", "seq_row", "seqid", "setAttrib", "setColnames", "setDimnames", "setrename", "setRownames", "settfm", "settfmv", "settransform", "settransformv", "slt", "slt<-", "smr", "ss", "STD", "t_list", "tfm", "tfm<-", "tfmv", "TRA", "unattrib", "unlist2d", "varying", "vclasses", "vlabels", "vlabels<-", "vtypes", "W", "wlddev")))

.COLLAPSE_GENERIC   <-   sort(unique(c("B","BY","D","Dlog","F","fsubset","fbetween","fdiff","ffirst","fgrowth","fHDbetween",
                           "fHDwithin","flag","flast","fmax","fmean","fmedian","fnth","fmin","fmode","varying",
                           "fNdistinct","fNobs","fprod","fscale","fsd","fsum","fvar","fwithin","funique",
                           "G","GRP","HDB","HDW","L","psacf","psccf","psmat","pspacf","qsu", "rsplit","fdroplevels",
                           "STD","TRA","W")))

.COLLAPSE_DATA <- c("GGDC10S", "wlddev")

.FAST_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
               "fmin","fmax","fnth","ffirst","flast","fNobs","fNdistinct",
               "fscale","fbetween","fwithin","fHDbetween","fHDwithin",
               "flag","fdiff","fgrowth")

.FAST_STAT_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
                    "fmin","fmax","fnth","ffirst","flast","fNobs","fNdistinct")

.OPERATOR_FUN <- c("STD","B","W","HDB","HDW","L","F","D","Dlog","G")
