

fsum <- function(x, ...) UseMethod("fsum") # , x

fsum.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, fill = FALSE, nthreads = .op[["nthreads"]], ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fsum.matrix(x, g, w, TRA, na.rm, use.g.names, fill = fill, nthreads = nthreads, ...))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fsum,x,0L,0L,w,na.rm,fill,nthreads))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(C_fsum,x,length(lev),g,w,na.rm,fill,nthreads), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fsum,x,fnlevels(g),g,w,na.rm,fill,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fsum,x,attr(g,"N.groups"),g,w,na.rm,fill,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(C_fsum,x,g[[1L]],g[[2L]],w,na.rm,fill,nthreads), GRPnames(g)))
    return(.Call(C_fsum,x,g[[1L]],g[[2L]],w,na.rm,fill,nthreads))
  }
  if(is.null(g)) return(TRAC(x,.Call(C_fsum,x,0L,0L,w,na.rm,fill,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(C_fsum,x,g[[1L]],g[[2L]],w,na.rm,fill,nthreads),g[[2L]],TRA, ...)
}

fsum.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, drop = TRUE, fill = FALSE, nthreads = .op[["nthreads"]], ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fsumm,x,0L,0L,w,na.rm,fill,drop,nthreads))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(C_fsumm,x,length(lev),g,w,na.rm,fill,FALSE,nthreads), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(C_fsumm,x,fnlevels(g),g,w,na.rm,fill,FALSE,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fsumm,x,attr(g,"N.groups"),g,w,na.rm,fill,FALSE,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(C_fsumm,x,g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(C_fsumm,x,g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads))
  }
  if(is.null(g)) return(TRAmC(x,.Call(C_fsumm,x,0L,0L,w,na.rm,fill,TRUE,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(C_fsumm,x,g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads),g[[2L]],TRA, ...)
}

fsum.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, drop = TRUE, fill = FALSE, nthreads = .op[["nthreads"]], ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fsuml,x,0L,0L,w,na.rm,fill,drop,nthreads))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(C_fsuml,x,length(lev),g,w,na.rm,fill,FALSE,nthreads), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fsuml,x,fnlevels(g),g,w,na.rm,fill,FALSE,nthreads))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fsuml,x,attr(g,"N.groups"),g,w,na.rm,fill,FALSE,nthreads))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads), groups))
    return(.Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads))
  }
  if(is.null(g)) return(TRAlC(x,.Call(C_fsuml,x,0L,0L,w,na.rm,fill,TRUE,nthreads),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads),g[[2L]],TRA, ...)
}

fsum.list <- function(x, ...) fsum.data.frame(x, ...)

fsum.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = FALSE,
                            keep.group_vars = TRUE, keep.w = TRUE, stub = .op[["stub"]], fill = FALSE,
                            nthreads = .op[["nthreads"]], ...) {
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
        if(nTRAl) sumw <- `names<-`(list(fsumC(w,g[[1L]],g[[2L]],NULL,na.rm,fill)), do_stub(stub, if(length(wsym) == 1L) as.character(wsym) else deparse(wsym), "sum.")) else if(keep.group_vars)
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
          return(setAttributes(c(g[[4L]], sumw, .Call(C_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(C_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads)), ax))
      } else return(setAttributes(.Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(C_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,fill,FALSE,nthreads),g[[2L]],TRA, ...))
}
