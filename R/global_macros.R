# global macros

.COLLAPSE_TOPICS <- c("collapse-documentation","fast-statistical-functions","quick-grouping",
                      "select-replace-vars","quick-conversion","advanced-data-aggregation",
                      "data-transformations","time-series-panel-series","list-processing",
                      "quick-summary","recode-replace","small-helpers")

.COLLAPSE_ALL    <-      c("%!in%","add_stub","all.identical","as.factor.GRP","as.numeric_factor",
                           "atomic_elem","atomic_elem<-","B","BY","cat_vars","cat_vars<-","char_vars",
                           "char_vars<-","collap","D","dapply","Date_vars","Date_vars<-","F","fact_vars",
                           "fact_vars<-","fbetween","fdiff","ffirst","fgrowth","fHDbetween","fHDwithin",
                           "flag","flast","fmax","fmean","fmedian","fmin","fmode","fNdistinct","fnlevels",
                           "fNobs","fprod","fscale","fsd","fsum","funique","fvar","fwithin","G","get_elem",
                           "get_vars","get_vars<-","GGDC10S","group_names.GRP","GRP","has_elem","HDB","HDW",
                           "irreg_elem","is.categorical","is.Date","is.GRP","is.regular","is.unlistable",
                           "L","ldepth","list_elem","list_elem<-","logi_vars","logi_vars<-","mctl","mrtl",
                           "na_insert","na_rm","namlab","num_vars","num_vars<-","psacf","psccf","psmat",
                           "pspacf","pwcor","pwcov","qDF","qDT","qF","qG","qM","qsu","rapply2d","Recode",
                           "reg_elem","replace_non_finite","replace_outliers","seq_col","seq_row",
                           "setColnames","setDimnames","setRownames","STD","TRA","unlist2d","vlabels",
                           "vlabels<-","W","wlddev")

.COLLAPSE_GENERIC   <-   c("B","BY","D","F","fbetween","fdiff","ffirst","fgrowth","fHDbetween",
                           "fHDwithin","flag","flast","fmax","fmean","fmedian","fmin","fmode",
                           "fNdistinct","fNobs","fprod","fscale","fsd","fsum","fvar","fwithin",
                           "G","GRP","HDB","HDW","L","psacf","psccf","psmat","pspacf","qsu",
                           "STD","TRA","W")

.COLLAPSE_DATA <- c("GGDC10S","wlddev")

.FAST_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
               "fmin","fmax","ffirst","flast","fNobs","fNdistinct",
               "fscale","fbetween","fwithin","fHDbetween","fHDwithin",
               "flag","fdiff","fgrowth")

.FAST_STAT_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
                    "fmin","fmax","ffirst","flast","fNobs","fNdistinct")

.OPERATOR_FUN <- c("STD","B","W","HDB","HDW","L","F","D","G")
