# Note: Adapted from fmode.R

fnth <- function(x, n = 0.5, ...) UseMethod("fnth") # , x

fnth.default <- function(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, ties = "q7", nthreads = .op[["nthreads"]], o = NULL, check.o = is.null(attr(o, "sorted")), ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fnth.matrix(x, n, g, w, TRA, na.rm, use.g.names, ties = ties, nthreads = nthreads, ...))
  if(!is.null(g)) g <- GRP(g, return.groups = use.g.names && is.null(TRA), call = FALSE) # sort = FALSE for TRA: not faster here...
  res <- .Call(C_fnth, x, n, g, w, na.rm, ties, nthreads, o, check.o)
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(res)
    if(use.g.names) names(res) <- GRPnames(g, FALSE)
    return(res)
  }
  TRAC(x,res,g[[2L]],TRA, ...)
}

fnth.matrix <- function(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, drop = TRUE, ties = "q7", nthreads = .op[["nthreads"]], ...) {
  if(!is.null(g)) g <- GRP(g, return.groups = use.g.names && is.null(TRA), call = FALSE) # sort = FALSE for TRA: not faster here...
  res <- .Call(C_fnthm, x, n, g, w, na.rm, drop, ties, nthreads)
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(res)
    if(use.g.names) dimnames(res)[[1L]] <- GRPnames(g)
    return(res)
  }
  TRAmC(x,res,g[[2L]],TRA, ...)
}

fnth.data.frame <- function(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, drop = TRUE, ties = "q7", nthreads = .op[["nthreads"]], ...) {
  if(!is.null(g)) g <- GRP(g, return.groups = use.g.names && is.null(TRA), call = FALSE) # sort = FALSE for TRA: not faster here...
  res <- .Call(C_fnthl, x, n, g, w, na.rm, drop, ties, nthreads)
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(if(drop) unlist(res) else res)
    if(use.g.names && !inherits(x, "data.table") && length(gn <- GRPnames(g)))
      attr(res, "row.names") <- gn
    return(res)
  }
  TRAlC(x,res,g[[2L]],TRA, ...)
}

fnth.list <- function(x, ...) fnth.data.frame(x, ...)

fnth.grouped_df <- function(x, n = 0.5, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = FALSE,
                             keep.group_vars = TRUE, keep.w = TRUE, stub = .op[["stub"]], ties = "q7", nthreads = .op[["nthreads"]], ...) {

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
        if(nTRAl) sumw <- `names<-`(list(fsumC(w,g[[1L]],g[[2L]],NULL,na.rm)), do_stub(stub, if(length(wsym) == 1L) as.character(wsym) else deparse(wsym), "sum.")) else if(keep.group_vars)
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
          return(setAttributes(c(g[[4L]], sumw, .Call(C_fnthl,x[-gn],n,g,w,na.rm,FALSE,ties,nthreads)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(C_fnthl,x[-gn],n,g,w,na.rm,FALSE,ties,nthreads)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(C_fnthl,x,n,g,w,na.rm,FALSE,ties,nthreads)), ax))
      } else return(setAttributes(.Call(C_fnthl,x,n,g,w,na.rm,FALSE,ties,nthreads), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(C_fnthl,x[-gn],n,g,w,na.rm,FALSE,ties,nthreads),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_fnthl,x[-gn],n,g,w,na.rm,FALSE,ties,nthreads),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_fnthl,x,n,g,w,na.rm,FALSE,ties,nthreads),g[[2L]],TRA, ...))
}


fmedian <- function(x, ...) UseMethod("fmedian") # , x

fmedian.default <- function(x, ..., ties = "mean")
  fnth.default(x, 0.5, ..., ties = ties)

fmedian.matrix <- function(x, ..., ties = "mean")
  fnth.matrix(x, 0.5, ..., ties = ties)

fmedian.data.frame <- function(x, ..., ties = "mean")
  fnth.data.frame(x, 0.5, ..., ties = ties)

fmedian.list <- fmedian.data.frame

fmedian.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = FALSE,
                            keep.group_vars = TRUE, keep.w = TRUE, stub = .op[["stub"]], ties = "mean", nthreads = .op[["nthreads"]], ...) {

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
        if(nTRAl) sumw <- `names<-`(list(fsumC(w,g[[1L]],g[[2L]],NULL,na.rm)), do_stub(stub, if(length(wsym) == 1L) as.character(wsym) else deparse(wsym), "sum.")) else if(keep.group_vars)
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
          return(setAttributes(c(g[[4L]], sumw, .Call(C_fnthl,x[-gn],0.5,g,w,na.rm,FALSE,ties,nthreads)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(C_fnthl,x[-gn],0.5,g,w,na.rm,FALSE,ties,nthreads)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(C_fnthl,x,0.5,g,w,na.rm,FALSE,ties,nthreads)), ax))
      } else return(setAttributes(.Call(C_fnthl,x,0.5,g,w,na.rm,FALSE,ties,nthreads), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(C_fnthl,x[-gn],0.5,g,w,na.rm,FALSE,ties,nthreads),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_fnthl,x[-gn],0.5,g,w,na.rm,FALSE,ties,nthreads),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_fnthl,x,0.5,g,w,na.rm,FALSE,1L,nthreads),g[[2L]],TRA, ...))
}
