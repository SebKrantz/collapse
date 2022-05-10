
# Note: for principal innovations of this code see fsum.R

fmode <- function(x, ...) UseMethod("fmode") # , x

fmode.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ties = "first", nthreads = 1L, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fmode.matrix(x, g, w, TRA, na.rm, use.g.names, ties = ties, nthreads = nthreads, ...))
  r <- switch(ties, first = 0L, min = 1L, max = 2L, last = 3L, stop("Unknown ties option: ", ties))
  if(!is.null(g)) g <- GRP(g, return.groups = use.g.names && is.null(TRA), call = FALSE) # sort = FALSE for TRA: not faster here...
  res <- .Call(C_fmode,x,g,w,na.rm,r,nthreads)
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(res)
    if(use.g.names) names(res) <- GRPnames(g, FALSE)
    return(res)
  }
  TRAC(x,res,g[[2L]],TRA, ...)
}

fmode.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ties = "first", nthreads = 1L, ...) {
  r <- switch(ties, first = 0L, min = 1L, max = 2L, last = 3L, stop("Unknown ties option: ", ties))
  if(!is.null(g)) g <- GRP(g, return.groups = use.g.names && is.null(TRA), call = FALSE) # sort = FALSE for TRA: not faster here...
  res <- .Call(C_fmodem,x,g,w,na.rm,drop,r,nthreads)
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(res)
    if(use.g.names) dimnames(res)[[1L]] <- GRPnames(g)
    return(res)
  }
  TRAmC(x,res,g[[2L]],TRA, ...)
}

fmode.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ties = "first", nthreads = 1L, ...) {
  r <- switch(ties, first = 0L, min = 1L, max = 2L, last = 3L, stop("Unknown ties option: ", ties))
  if(!is.null(g)) g <- GRP(g, return.groups = use.g.names && is.null(TRA), call = FALSE) # sort = FALSE for TRA: not faster here...
  res <- .Call(C_fmodel,x,g,w,na.rm,r,nthreads)
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(if(drop) unlist(res) else res)
    if(use.g.names && !inherits(x, "data.table") && length(gn <- GRPnames(g)))
      attr(res, "row.names") <- gn
    return(res)
  }
  TRAlC(x,res,g[[2L]],TRA, ...)
}

fmode.list <- function(x, ...) fmode.data.frame(x, ...)

fmode.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                             keep.group_vars = TRUE, keep.w = TRUE, ties = "first", nthreads = 1L, ...) {
  r <- switch(ties, first = 0L, min = 1L, max = 2L, last = 3L, stop("Unknown ties option: ", ties))
  g <- GRP.grouped_df(x, call = FALSE)
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

  gl <- length(gn) > 0L

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
          return(setAttributes(c(g[[4L]], sumw, .Call(C_fmodel,x[-gn],g,w,na.rm,r,nthreads)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(C_fmodel,x[-gn],g,w,na.rm,r,nthreads)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(C_fmodel,x,g,w,na.rm,r,nthreads)), ax))
      } else return(setAttributes(.Call(C_fmodel,x,g,w,na.rm,r,nthreads), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(C_fmodel,x[-gn],g,w,na.rm,r,nthreads),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_fmodel,x[-gn],g,w,na.rm,r,nthreads),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_fmodel,x,g,w,na.rm,r,nthreads),g[[2L]],TRA, ...))
}
