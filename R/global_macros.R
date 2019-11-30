# global macros

.COLLAPSE_TOPICS <- c("collapse-documentation","fast-statistical-functions","quick-grouping",
                      "select-replace-vars","quick-conversion","advanced-data-aggregation",
                      "data-transformations","time-series-panel-series","list-processing",
                      "quick-summary","recode-replace","small-helpers")

.COLLAPSE_ALL    <-      c("add_stub","all.identical","as.factor.GRP","atomic_elem",
                           "atomic_elem<-","B","BY","cat_vars","cat_vars<-","char_vars",
                           "char_vars<-","collap","D","dapply","Date_vars","Date_vars<-",
                           "F","fact_vars","fact_vars<-","fdiff","ffirst","fgrowth","flag",
                           "flast","fmax","fmean","fmedian","fmin","fmode","fNdistinct",
                           "fNobs","fprod","fscale","fsd","fsum","fvar","fbetween","fwithin",
                           "fHDbetween","fHDwithin","Recode","replace_non_finite","replace_outliers",
                           "G","get_elem",
                           "get_vars","get_vars<-","GGDC","group_names.GRP","GRP","has_elem",
                           "HDB","HDW","irreg_elem","is.categorical","is.Date","is.GRP",
                           "is.regular","is.unlistable","L","ldepth","list_elem",
                           "list_elem<-","logi_vars","logi_vars<-","mctl","mrtl","namlab",
                           "num_vars","num_vars<-","psacf","psccf","psmat","pspacf","qDF",
                           "qDT","qF","qG","qM","qsu","rapply2d","reg_elem","STD","TRA",
                           "unlist2d","vlabels","vlabels<-","W","wlddev")

.COLLAPSE_GENERIC   <-   c("B","BY","D","F","fdiff","ffirst","fgrowth","flag",
                           "flast","fmax","fmean","fmedian","fmin","fmode","fNdistinct",
                           "fNobs","fprod","fscale","fsd","fsum","fvar",
                           "fbetween","fwithin","fHDbetween","fHDwithin","G","GRP",
                           "HDB","HDW","L","psacf","psccf","psmat","pspacf","qsu",
                           "STD","TRA","W")

.COLLAPSE_DATA <- c("GGDC10S","wlddev")

.FAST_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
               "fmin","fmax","ffirst","flast","fNobs","fNdistinct",
               "fscale","fbetween","fwithin","fHDbetween","fHDwithin",
               "flag","fdiff","fgrowth")

.FAST_STAT_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
                    "fmin","fmax","ffirst","flast","fNobs","fNdistinct")

.OPERATOR_FUN <- c("STD","B","W","HDB","HDW","L","F","D","G")
