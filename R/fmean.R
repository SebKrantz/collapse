# library(Rcpp)
# sourceCpp('src/fsum.cpp')
# sourceCpp('src/fmean.cpp')
# sourceCpp('src/fmeana.cpp')
# sourceCpp('src/fmeanl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# Note: for principal innovations of this code see fsum.R !!

fmean <- function(x, ...) UseMethod("fmean") # , x

fmean.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fmean,x,0L,0L,NULL,w,na.rm)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fmean,x,length(lev),g,NULL,w,na.rm), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fmean,x,fnlevels(g),g,NULL,w,na.rm)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fmean,x,attr(g,"N.groups"),g,NULL,w,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fmean,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_fmean,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmean,x,0L,0L,NULL,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmean,x,fnlevels(g),g,NULL,w,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fmean,x,attr(g,"N.groups"),g,NULL,w,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fmean,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmean.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fmeanm,x,0L,0L,NULL,w,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fmeanm,x,length(lev),g,NULL,w,na.rm,drop), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fmeanm,x,fnlevels(g),g,NULL,w,na.rm,drop)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fmeanm,x,attr(g,"N.groups"),g,NULL,w,na.rm,drop))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fmeanm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fmeanm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fmeanm,x,0L,0L,NULL,w,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fmeanm,x,fnlevels(g),g,NULL,w,na.rm,drop),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fmeanm,x,attr(g,"N.groups"),g,NULL,w,na.rm,drop),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fmeanm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmean.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fmeanl,x,0L,0L,NULL,w,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fmeanl,x,length(lev),g,NULL,w,na.rm,drop), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fmeanl,x,fnlevels(g),g,NULL,w,na.rm,drop)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fmeanl,x,attr(g,"N.groups"),g,NULL,w,na.rm,drop))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop), groups)) else
          return(.Call(Cpp_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fmeanl,x,0L,0L,NULL,w,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fmeanl,x,fnlevels(g),g,NULL,w,na.rm,drop),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fmeanl,x,attr(g,"N.groups"),g,NULL,w,na.rm,drop),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmean.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                             keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  g <- GRP.grouped_df(x)
  wsym <- l1orn(all.vars(substitute(w)), "NULL") # deparse
  nam <- attr(x, "names")
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  sumw <- NULL

  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- unclass(x)[[wn]] # faster using unclass??
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    gn <- c(gn, wn)
    if(keep.w) {
      if(nTRAl) sumw <- `names<-`(list(fsumCpp(w,g[[1L]],g[[2L]],NULL,na.rm)), paste0("sum.", wsym)) else if(keep.group_vars)
        gn2 <- gn else sumw <- gn2 <- wn
    }
  }

  gl <- length(gn) > 0L # necessary here, not before !!!

  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], names(sumw), nam[-gn])
          return(setAttributes(c(g[[4L]], sumw, .Call(Cpp_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE)), ax))
        } else {
          ax[["names"]] <- c(names(sumw), nam[-gn])
          return(setAttributes(c(sumw, .Call(Cpp_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE)), ax))
        }
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(Cpp_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE)), ax))
      } else return(setAttributes(.Call(Cpp_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- nam[-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
}

# Previous Version: With deparse(substitute(w)) and only keeping grouping columns if found in x.
# fmean.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
#                              keep.group_vars = TRUE, keep.w = TRUE, ...) {
#   if(!missing(...)) unused_arg_warning(match.call(), ...)
#   g <- GRP.grouped_df(x)
#   wsym <- deparse(substitute(w))
#   nam <- attr(x, "names")
#   gn2 <- gn <- which(nam %in% g[[5L]])
#   nTRAl <- is.null(TRA)
#   sumw <- NULL
#
#   if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
#     w <- unclass(x)[[wn]] # faster using unclass??
#     if(any(gn == wn)) stop("Weights coincide with grouping variables!")
#     onlyw <- !length(gn)
#     gn <- c(gn, wn)
#     if(keep.w) {
#       if(nTRAl) sumw <- `names<-`(list(fsumCpp(w,g[[1L]],g[[2L]],NULL,na.rm)), paste0("sum.", wsym)) else if(keep.group_vars)
#         gn2 <- gn else sumw <- gn2 <- wn
#     }
#   } else onlyw <- FALSE
#
#   gl <- length(gn) > 0L # necessary here, not before !!!
#
#   if(gl || nTRAl) {
#     ax <- attributes(x)
#     attributes(x) <- NULL
#     if(nTRAl) {
#       ax[["groups"]] <- NULL
#       ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
#       ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
#       if(gl) {
#         if(keep.group_vars && !onlyw) {
#           ax[["names"]] <- c(g[[5L]], names(sumw), ax[["names"]][-gn])
#           return(setAttributes(c(g[[4L]], sumw, .Call(Cpp_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE)), ax))
#         } else {
#           ax[["names"]] <- c(names(sumw), ax[["names"]][-gn])
#           return(setAttributes(c(sumw, .Call(Cpp_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE)), ax))
#         }
#       } else return(setAttributes(.Call(Cpp_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE), ax))
#     } else if(keep.group_vars || (keep.w && length(sumw))) {
#       ax[["names"]] <- c(ax[["names"]][gn2], ax[["names"]][-gn])
#       return(setAttributes(c(x[gn2],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
#     } else {
#       ax[["names"]] <- ax[["names"]][-gn]
#       return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
#     }
#   } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
# }
