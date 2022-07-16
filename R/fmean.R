
# Note: for principal innovations of this code see fsum.R

fmean <- function(x, ...) UseMethod("fmean") # , x

fmean.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, nthreads = 1L, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fmean.matrix(x, g, w, TRA, na.rm, use.g.names, nthreads = nthreads, ...))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fmean,x,0L,0L,NULL,w,na.rm,nthreads))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(C_fmean,x,length(lev),g,NULL,w,na.rm,nthreads), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fmean,x,fnlevels(g),g,NULL,w,na.rm,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fmean,x,attr(g,"N.groups"),g,NULL,w,na.rm,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(C_fmean,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,nthreads), GRPnames(g)))
    return(.Call(C_fmean,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,nthreads))
  }
  if(is.null(g)) return(TRAC(x,.Call(C_fmean,x,0L,0L,NULL,w,na.rm,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(C_fmean,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,nthreads),g[[2L]],TRA, ...)
}

fmean.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, nthreads = 1L, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fmeanm,x,0L,0L,NULL,w,na.rm,drop,nthreads))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(C_fmeanm,x,length(lev),g,NULL,w,na.rm,drop,nthreads), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(C_fmeanm,x,fnlevels(g),g,NULL,w,na.rm,drop,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fmeanm,x,attr(g,"N.groups"),g,NULL,w,na.rm,drop,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(C_fmeanm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop,nthreads), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(C_fmeanm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop,nthreads))
  }
  if(is.null(g)) return(TRAmC(x,.Call(C_fmeanm,x,0L,0L,NULL,w,na.rm,TRUE,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(C_fmeanm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop,nthreads),g[[2L]],TRA, ...)
}

fmean.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, nthreads = 1L, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fmeanl,x,0L,0L,NULL,w,na.rm,drop,nthreads))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(C_fmeanl,x,length(lev),g,NULL,w,na.rm,drop,nthreads), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fmeanl,x,fnlevels(g),g,NULL,w,na.rm,drop,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fmeanl,x,attr(g,"N.groups"),g,NULL,w,na.rm,drop,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
    return(setRnDF(.Call(C_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop,nthreads), groups))
    return(.Call(C_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop,nthreads))
  }
  if(is.null(g)) return(TRAlC(x,.Call(C_fmeanl,x,0L,0L,NULL,w,na.rm,TRUE,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(C_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop,nthreads),g[[2L]],TRA, ...)
}

fmean.list <- function(x, ...) fmean.data.frame(x, ...)

fmean.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                             keep.group_vars = TRUE, keep.w = TRUE, nthreads = 1L, ...) {
  g <- GRP.grouped_df(x, call = FALSE)
  if(is.null(g[[4L]])) keep.group_vars <- FALSE
  wsym <- substitute(w)
  nam <- attr(x, "names")
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  sumw <- NULL

  if(!is.null(wsym)) {
    w <- eval(wsym, x, parent.frame())
    if(length(wn <- which(nam %in% all.vars(wsym)))) {
      if(any(gn %in% wn)) stop("Weights coincide with grouping variables!")
      gn <- c(gn, wn)
      if(keep.w) {
        if(nTRAl) sumw <- `names<-`(list(fsumC(w,g[[1L]],g[[2L]],NULL,na.rm)), paste0("sum.", if(length(wsym) == 1L) wsym else deparse(wsym))) else if(keep.group_vars)
          gn2 <- gn else sumw <- gn2 <- wn
      }
    }
  }

  gl <- length(gn) > 0L # necessary here, not before !

  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      if(!missing(...)) unused_arg_action(match.call(), ...)
      ax[["groups"]] <- NULL
      ax[["class"]] <- fsetdiff(ax[["class"]], c("GRP_df", "grouped_df"))
      ax[["row.names"]] <- if(use.g.names) GRPnames(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], names(sumw), nam[-gn])
          return(setAttributes(c(g[[4L]], sumw, .Call(C_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,nthreads)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(C_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,nthreads)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(C_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,nthreads)), ax))
      } else return(setAttributes(.Call(C_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,nthreads), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(C_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,nthreads),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,nthreads),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,nthreads),g[[2L]],TRA, ...))
}

# Previous Version: With deparse(substitute(w)) and only keeping grouping columns if found in x.
# fmean.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
#                              keep.group_vars = TRUE, keep.w = TRUE, ...) {
#   if(!missing(...)) unused_arg_action(match.call(), ...)
#   g <- GRP.grouped_df(x, call = FALSE)
#   wsym <- deparse(substitute(w))
#   nam <- attr(x, "names")
#   gn2 <- gn <- which(nam %in% g[[5L]])
#   nTRAl <- is.null(TRA)
#   sumw <- NULL
#
#   if(length(wsym) && length(wn <- whichv(nam, wsym))) {
#     w <- .subset2(x, wn) # faster using unclass??
#     if(any(gn %in% wn)) stop("Weights coincide with grouping variables!")
#     onlyw <- !length(gn)
#     gn <- c(gn, wn)
#     if(keep.w) {
#       if(nTRAl) sumw <- `names<-`(list(fsumC(w,g[[1L]],g[[2L]],NULL,na.rm)), paste0("sum.", if(length(wsym) == 1L) wsym else deparse(wsym))) else if(keep.group_vars)
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
#       ax[["class"]] <- fsetdiff(ax[["class"]], c("GRP_df", "grouped_df"))
#       ax[["row.names"]] <- if(use.g.names) GRPnames(g) else .set_row_names(g[[1L]])
#       if(gl) {
#         if(keep.group_vars && !onlyw) {
#           ax[["names"]] <- c(g[[5L]], names(sumw), ax[["names"]][-gn])
#           return(setAttributes(c(g[[4L]], sumw, .Call(C_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE)), ax))
#         } else {
#           ax[["names"]] <- c(names(sumw), ax[["names"]][-gn])
#           return(setAttributes(c(sumw, .Call(C_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE)), ax))
#         }
#       } else return(setAttributes(.Call(C_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE), ax))
#     } else if(keep.group_vars || (keep.w && length(sumw))) {
#       ax[["names"]] <- c(ax[["names"]][gn2], ax[["names"]][-gn])
#       return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(C_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE),g[[2L]],TRA, ...)), ax))
#     } else {
#       ax[["names"]] <- ax[["names"]][-gn]
#       return(setAttributes(TRAlC(x[-gn],.Call(C_fmeanl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE),g[[2L]],TRA, ...), ax))
#     }
#   } else return(TRAlC(x,.Call(C_fmeanl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE),g[[2L]],TRA, ...))
# }
