
# Note: for principal innovations of this code see fsum.R

fmode <- function(x, ...) UseMethod("fmode") # , x


fmode.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ties = "first", ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fmode.matrix(x, g, w, TRA, na.rm, use.g.names, ties = ties, ...))
  ret <- switch(ties, first = 0L, min = 1L, max = 2L, last = 3L, stop("Unknown ties option: ", ties))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fmode,x,0L,0L,NULL,w,na.rm,ret))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fmode,x,length(lev),g,NULL,w,na.rm,ret), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fmode,x,fnlevels(g),g,NULL,w,na.rm,ret))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fmode,x,attr(g,"N.groups"),g,NULL,w,na.rm,ret))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(Cpp_fmode,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret), GRPnames(g)))
    return(.Call(Cpp_fmode,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret))
  }
  if(is.null(g)) return(TRAC(x,.Call(Cpp_fmode,x,0L,0L,NULL,w,na.rm,ret),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(Cpp_fmode,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret),g[[2L]],TRA, ...)
}

fmode.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ties = "first", ...) {
  ret <- switch(ties, first = 0L, min = 1L, max = 2L, last = 3L, stop("Unknown ties option: ", ties))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fmodem,x,0L,0L,NULL,w,na.rm,drop,ret))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fmodem,x,length(lev),g,NULL,w,na.rm,FALSE,ret), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fmodem,x,fnlevels(g),g,NULL,w,na.rm,FALSE,ret))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fmodem,x,attr(g,"N.groups"),g,NULL,w,na.rm,FALSE,ret))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(Cpp_fmodem,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(Cpp_fmodem,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret))
  }
  if(is.null(g)) return(TRAmC(x,.Call(Cpp_fmodem,x,0L,0L,NULL,w,na.rm,TRUE,ret),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(Cpp_fmodem,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE,ret),g[[2L]],TRA, ...)
}

fmode.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ties = "first", ...) {
  ret <- switch(ties, first = 0L, min = 1L, max = 2L, last = 3L, stop("Unknown ties option: ", ties))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) if(drop) return(unlist(.Call(Cpp_fmodel,x,0L,0L,NULL,w,na.rm,ret))) else return(.Call(Cpp_fmodel,x,0L,0L,NULL,w,na.rm,ret))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(Cpp_fmodel,x,length(lev),g,NULL,w,na.rm,ret), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fmodel,x,fnlevels(g),g,NULL,w,na.rm,ret))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fmodel,x,attr(g,"N.groups"),g,NULL,w,na.rm,ret))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret), groups))
    return(.Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret))
  }
  if(is.null(g)) return(TRAlC(x,.Call(Cpp_fmodel,x,0L,0L,NULL,w,na.rm,ret),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret),g[[2L]],TRA, ...)
}

fmode.list <- function(x, ...) fmode.data.frame(x, ...)

fmode.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                             keep.group_vars = TRUE, keep.w = TRUE, ties = "first", ...) {
  ret <- switch(ties, first = 0L, min = 1L, max = 2L, last = 3L, stop("Unknown ties option: ", ties))
  g <- GRP.grouped_df(x, call = FALSE)
  wsym <- l1orn(as.character(substitute(w)), NULL)
  nam <- attr(x, "names")
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  sumw <- NULL

  if(length(wsym) && length(wn <- whichv(nam, wsym))) {
    w <- .subset2(x, wn)
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    gn <- c(gn, wn)
    if(keep.w) {
      if(nTRAl) sumw <- `names<-`(list(fsumC(w,g[[1L]],g[[2L]],NULL,na.rm)), paste0("sum.", wsym)) else if(keep.group_vars)
        gn2 <- gn else sumw <- gn2 <- wn
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
          return(setAttributes(c(g[[4L]], sumw, .Call(Cpp_fmodel,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(Cpp_fmodel,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret)), ax))
      } else return(setAttributes(.Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(Cpp_fmodel,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(Cpp_fmodel,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ret),g[[2L]],TRA, ...))
}
