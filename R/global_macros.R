# global macros

.COLLAPSE_TOPICS <- c("collapse-documentation","fast-statistical-functions","quick-grouping",
                      "select-replace-vars","quick-conversion","advanced-data-aggregation",
                      "data-transformations","time-series-panel-series","list-processing",
                      "quick-summary","recode-replace","small-helpers")
#
# obj <- objects("package:collapse")
# export <- setdiff(obj[-grep("Cpp|_collapse_", obj)],
#                   c("change_row_names","cols2int","cond_duplicate_attributes","at2GRP",
#                     "G_t","give_nam","library.dynam.unload","remove_attributes","setRow.names",
#                     "system.file","TRAtoInt","View.qsu","group.names.GRP","na.rm","pwcor","pwcov",
#                     "seq_col","seq_row","setAttributes","setColnames","setDimnames","setRownames","duplAttributes",
#                     "fnlevels","anyNAerror","charorNULL","collapse_init","complete.cases","cols2log",
#                     "colsubset","cond_duplattributes","cond_duplAttributes","condsetn","cor","cov","Csplit",
#                     "duplattributes","hist","list_extract_FUN","list_extract_ind","list_extract_names",
#                     "list_extract_regex","myModFrame","NCOL2","NROW2","par","plot","rbindlist","as.formula",
#                     "forder","forderv","getfl","getPartData","init_collapse","rgrep","setAttr","setattr_clp",
#                     "setattributes","setcolorder","setNames","subsetDT","subsetfl","subsetVector","terms.formula",
#                     "uniqlengths","unique_factor","frank"))
# cat(paste0("export(", export,")\n"))
#
# methods <- obj[grep(".default|.matrix|.data.frame|.pseries|.pdata.frame|.grouped_df|GRP.|print.|plot.", obj)]
# methods <- unlist(lapply(strsplit(methods,".", fixed = TRUE), function(x) paste0("S3method(",x[1],", ",paste(x[-1], collapse = "."),")\n" )))
# cat(methods)



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
