# global macros

.COLLAPSE_TOPICS <- c("collapse-documentation","A1-fast-statistical-functions","A2-fast-grouping",
                      "A3-data-frame-manipulation","A4-quick-conversion","A5-advanced-aggregation",
                      "A6-data-transformations","A7-time-series-panel-series","A8-list-processing",
                      "A9-summary-statistics","AA1-recode-replace","AA2-small-helpers")

# Get updated .COLLAPSE_ALL:
# v <- grep(".default$|.matrix$|.data.frame$", objects("package:collapse"), invert = TRUE, value = TRUE)
# cat(v, sep = '", "')

# all package objects..
# allobj <- ls(getNamespace("collapse"), all.names=TRUE)


.COLLAPSE_ALL <- sort(unique(c("%!in%", "add_stub", "add_vars", "add_vars<-", "all_identical", "all_obj_equal",
                               "as.character_factor", "as.factor.GRP", "as.numeric_factor", "atomic_elem",
                               "atomic_elem<-", "av", "av<-", "B", "BY", "cat_vars", "cat_vars<-", "char_vars",
                               "char_vars<-", "ckmatch", "collap", "collapg", "collapv", "D", "dapply", "Date_vars",
                               "Date_vars<-", "descr", "Dlog", "F", "fact_vars", "fact_vars<-", "fbetween",
                               "fcompute", "fdiff", "fdim", "ffirst", "fFtest", "fgroup_by", "fgroup_vars",
                               "fgrowth", "fHDbetween", "fHDwithin", "finteraction", "flag", "flast",
                               "fmax", "fmean", "fmedian", "fmin", "fmode", "fncol", "fNdistinct",
                               "fnlevels", "fNobs", "fnrow", "fprod", "fscale", "fsd", "fselect",
                               "fselect<-", "fsubset", "fsum", "ftransform", "funique", "fvar",
                               "fwithin", "G", "get_elem", "get_vars", "get_vars<-", "GGDC10S",
                               "group_names.GRP", "groupid", "GRP", "gv", "gv<-", "has_elem", "HDB",
                               "HDW", "irreg_elem", "is.categorical", "is.Date", "is.GRP", "is.qG",
                               "is.regular", "is.unlistable", "L", "ldepth", "list_elem", "list_elem<-",
                               "logi_vars", "logi_vars<-", "mctl", "mrtl", "na_insert", "na_omit",
                               "na_rm", "namlab", "num_vars", "num_vars<-", "nv", "nv<-", "print.pwcor",
                               "print.pwcov", "print.qsu", "psacf", "psccf", "psmat", "pspacf", "pwcor",
                               "pwcov", "pwNobs", "qDF", "qDT", "qF", "qG", "qM", "qsu", "radixorder",
                               "radixorderv", "rapply2d", "Recode", "recode_char", "recode_num",
                               "reg_elem", "replace_Inf", "replace_NA", "replace_non_finite", "replace_outliers",
                               "rm_stub", "sbt", "seq_col", "seq_row", "seqid", "setColnames", "setDimnames",
                               "setRownames", "settfm", "settransform", "slt", "slt<-", "ss", "STD", "tfm",
                               "TRA", "unattrib", "unlist2d", "varying", "vclasses", "vlabels", "vlabels<-",
                               "vtypes", "W", "wlddev")))

.COLLAPSE_GENERIC   <-   sort(unique(c("B","BY","D","Dlog","F","fsubset","fbetween","fdiff","ffirst","fgrowth","fHDbetween",
                           "fHDwithin","flag","flast","fmax","fmean","fmedian","fmin","fmode","varying",
                           "fNdistinct","fNobs","fprod","fscale","fsd","fsum","fvar","fwithin",
                           "G","GRP","HDB","HDW","L","psacf","psccf","psmat","pspacf","qsu",
                           "STD","TRA","W")))

.COLLAPSE_DATA <- c("GGDC10S","wlddev")

.FAST_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
               "fmin","fmax","ffirst","flast","fNobs","fNdistinct",
               "fscale","fbetween","fwithin","fHDbetween","fHDwithin",
               "flag","fdiff","fgrowth")

.FAST_STAT_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
                    "fmin","fmax","ffirst","flast","fNobs","fNdistinct")

.OPERATOR_FUN <- c("STD","B","W","HDB","HDW","L","F","D","Dlog","G")
