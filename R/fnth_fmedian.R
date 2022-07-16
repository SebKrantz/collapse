# Note: for foundational changes to this code see fsum.R

fnth <- function(x, n = 0.5, ...) UseMethod("fnth") # , x

fnth.default <- function(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ties = "mean", nthreads = 1L, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fnth.matrix(x, n, g, w, TRA, na.rm, use.g.names, ties = ties, nthreads = nthreads, ...))
  ret <- switch(ties, mean = 1L, min = 2L, max = 3L, stop("ties must be 'mean', 'min' or 'max'"))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fnth,x,n,0L,0L,NULL,w,na.rm,ret,nthreads))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fnth,x,n,length(lev),g,NULL,w,na.rm,ret,nthreads), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fnth,x,n,fnlevels(g),g,NULL,w,na.rm,ret,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fnth,x,n,attr(g,"N.groups"),g,NULL,w,na.rm,ret,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(Cpp_fnth,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret,nthreads), GRPnames(g)))
    return(.Call(Cpp_fnth,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret,nthreads))
  }
  if(is.null(g)) return(TRAC(x,.Call(Cpp_fnth,x,n,0L,0L,NULL,w,na.rm,ret,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(Cpp_fnth,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret,nthreads),g[[2L]],TRA, ...)
}

fnth.matrix <- function(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ties = "mean", nthreads = 1L, ...) {
  ret <- switch(ties, mean = 1L, min = 2L, max = 3L, stop("ties must be 'mean', 'min' or 'max'"))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fnthm,x,n,0L,0L,NULL,w,na.rm,drop,ret,nthreads))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fnthm,x,n,length(lev),g,NULL,w,na.rm,FALSE,ret,nthreads), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fnthm,x,n,fnlevels(g),g,NULL,w,na.rm,FALSE,ret,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fnthm,x,n,attr(g,"N.groups"),g,NULL,w,na.rm,FALSE,ret,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(Cpp_fnthm,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(Cpp_fnthm,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads))
  }
  if(is.null(g)) return(TRAmC(x,.Call(Cpp_fnthm,x,n,0L,0L,NULL,w,na.rm,TRUE,ret,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(Cpp_fnthm,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads),g[[2L]],TRA, ...)
}

fnth.data.frame <- function(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ties = "mean", nthreads = 1L, ...) {
  ret <- switch(ties, mean = 1L, min = 2L, max = 3L, stop("ties must be 'mean', 'min' or 'max'"))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fnthl,x,n,0L,0L,NULL,w,na.rm,drop,ret,nthreads))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(Cpp_fnthl,x,n,length(lev),g,NULL,w,na.rm,FALSE,ret,nthreads), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fnthl,x,n,fnlevels(g),g,NULL,w,na.rm,FALSE,ret,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fnthl,x,n,attr(g,"N.groups"),g,NULL,w,na.rm,FALSE,ret,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(Cpp_fnthl,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads), groups))
    return(.Call(Cpp_fnthl,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads))
  }
  if(is.null(g)) return(TRAlC(x,.Call(Cpp_fnthl,x,n,0L,0L,NULL,w,na.rm,TRUE,ret,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(Cpp_fnthl,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads),g[[2L]],TRA, ...)
}

fnth.list <- function(x, ...) fnth.data.frame(x, ...)

fnth.grouped_df <- function(x, n = 0.5, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                             keep.group_vars = TRUE, keep.w = TRUE, ties = "mean", nthreads = 1L, ...) {
  ret <- switch(ties, mean = 1L, min = 2L, max = 3L, stop("ties must be 'mean', 'min' or 'max'"))
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
          return(setAttributes(c(g[[4L]], sumw, .Call(Cpp_fnthl,x[-gn],n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(Cpp_fnthl,x[-gn],n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(Cpp_fnthl,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads)), ax))
      } else return(setAttributes(.Call(Cpp_fnthl,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(Cpp_fnthl,x[-gn],n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(Cpp_fnthl,x[-gn],n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(Cpp_fnthl,x,n,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret,nthreads),g[[2L]],TRA, ...))
}




fmedian <- function(x, ...) UseMethod("fmedian") # , x

fmedian.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, nthreads = 1L, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fmedian.matrix(x, g, w, TRA, na.rm, use.g.names, nthreads = nthreads, ...))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fnth,x,0.5,0L,0L,NULL,w,na.rm,1L,nthreads))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fnth,x,0.5,length(lev),g,NULL,w,na.rm,1L,nthreads), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fnth,x,0.5,fnlevels(g),g,NULL,w,na.rm,1L,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fnth,x,0.5,attr(g,"N.groups"),g,NULL,w,na.rm,1L,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(Cpp_fnth,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,1L,nthreads), GRPnames(g)))
    return(.Call(Cpp_fnth,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,1L,nthreads))
  }
  if(is.null(g)) return(TRAC(x,.Call(Cpp_fnth,x,0.5,0L,0L,NULL,w,na.rm,1L,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(Cpp_fnth,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,1L,nthreads),g[[2L]],TRA, ...)
}

fmedian.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, nthreads = 1L, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fnthm,x,0.5,0L,0L,NULL,w,na.rm,drop,1L,nthreads))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fnthm,x,0.5,length(lev),g,NULL,w,na.rm,FALSE,1L,nthreads), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fnthm,x,0.5,fnlevels(g),g,NULL,w,na.rm,FALSE,1L,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fnthm,x,0.5,attr(g,"N.groups"),g,NULL,w,na.rm,FALSE,1L,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(Cpp_fnthm,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(Cpp_fnthm,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads))
  }
  if(is.null(g)) return(TRAmC(x,.Call(Cpp_fnthm,x,0.5,0L,0L,NULL,w,na.rm,TRUE,1L,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(Cpp_fnthm,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads),g[[2L]],TRA, ...)
}

fmedian.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, nthreads = 1L, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fnthl,x,0.5,0L,0L,NULL,w,na.rm,drop,1L,nthreads))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(Cpp_fnthl,x,0.5,length(lev),g,NULL,w,na.rm,FALSE,1L,nthreads), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fnthl,x,0.5,fnlevels(g),g,NULL,w,na.rm,FALSE,1L,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fnthl,x,0.5,attr(g,"N.groups"),g,NULL,w,na.rm,FALSE,1L,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(Cpp_fnthl,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads), groups))
    return(.Call(Cpp_fnthl,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads))
  }
  if(is.null(g)) return(TRAlC(x,.Call(Cpp_fnthl,x,0.5,0L,0L,NULL,w,na.rm,TRUE,1L,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(Cpp_fnthl,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads),g[[2L]],TRA, ...)
}

fmedian.list <- function(x, ...) fmedian.data.frame(x, ...)

fmedian.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
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
          return(setAttributes(c(g[[4L]], sumw, .Call(Cpp_fnthl,x[-gn],0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(Cpp_fnthl,x[-gn],0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(Cpp_fnthl,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads)), ax))
      } else return(setAttributes(.Call(Cpp_fnthl,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(Cpp_fnthl,x[-gn],0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(Cpp_fnthl,x[-gn],0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(Cpp_fnthl,x,0.5,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,1L,nthreads),g[[2L]],TRA, ...))
}
