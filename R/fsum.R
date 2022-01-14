

fsum <- function(x, ...) UseMethod("fsum") # , x

fsum.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fsum.matrix(x, g, w, TRA, na.rm, use.g.names, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(C_fsum,x,0L,0L,w,na.rm))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(C_fsum,x,length(lev),g,w,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fsum,x,fnlevels(g),g,w,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fsum,x,attr(g,"N.groups"),g,w,na.rm))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(C_fsum,x,g[[1L]],g[[2L]],w,na.rm), GRPnames(g)))
    return(.Call(C_fsum,x,g[[1L]],g[[2L]],w,na.rm))
  }
  if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(C_fsum,x,0L,0L,w,na.rm),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(C_fsum,x,fnlevels(g),g,w,na.rm),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRA,x,.Call(C_fsum,x,attr(g,"N.groups"),g,w,na.rm),g,TtI(TRA)))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRA,x,.Call(C_fsum,x,g[[1L]],g[[2L]],w,na.rm),g[[2L]],TtI(TRA))
}

fsum.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(C_fsumm,x,0L,0L,w,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(C_fsumm,x,length(lev),g,w,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(C_fsumm,x,fnlevels(g),g,w,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fsumm,x,attr(g,"N.groups"),g,w,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(C_fsumm,x,g[[1L]],g[[2L]],w,na.rm,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(C_fsumm,x,g[[1L]],g[[2L]],w,na.rm,FALSE))
  }
  if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(C_fsumm,x,0L,0L,w,na.rm,TRUE),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(C_fsumm,x,fnlevels(g),g,w,na.rm,FALSE),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRAm,x,.Call(C_fsumm,x,attr(g,"N.groups"),g,w,na.rm,FALSE),g,TtI(TRA)))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRAm,x,.Call(C_fsumm,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TtI(TRA))
}

fsum.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(C_fsuml,x,0L,0L,w,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(C_fsuml,x,length(lev),g,w,na.rm,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fsuml,x,fnlevels(g),g,w,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fsuml,x,attr(g,"N.groups"),g,w,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE), groups))
    return(.Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE))
  }
  if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(C_fsuml,x,0L,0L,w,na.rm,TRUE),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(C_fsuml,x,fnlevels(g),g,w,na.rm,FALSE),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRAl,x,.Call(C_fsuml,x,attr(g,"N.groups"),g,w,na.rm,FALSE),g,TtI(TRA)))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRAl,x,.Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TtI(TRA))
}

fsum.list <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
  fsum.data.frame(x, g, w, TRA, na.rm, use.g.names, drop, ...)

fsum.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                             keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  wsym <- l1orn(as.character(substitute(w)), NULL)
  nam <- attr(x, "names")
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  sumw <- NULL

  if(length(wsym) && length(wn <- whichv(nam, wsym))) {
    w <- .subset2(x, wn) # faster using unclass?
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    gn <- c(gn, wn)
    if(keep.w) {
      if(nTRAl) sumw <- `names<-`(list(fsumC(w,g[[1L]],g[[2L]],NULL,na.rm)), paste0("sum.", wsym)) else if(keep.group_vars)
        gn2 <- gn else sumw <- gn2 <- wn
    }
  }

  gl <- length(gn) > 0L # necessary here, not before !

  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- fsetdiff(ax[["class"]], c("GRP_df", "grouped_df"))
      ax[["row.names"]] <- if(use.g.names) GRPnames(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], names(sumw), nam[-gn])
          return(setAttributes(c(g[[4L]], sumw, .Call(C_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(C_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
      } else return(setAttributes(.Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],.Call(Cpp_TRAl,x[-gn],.Call(C_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TtI(TRA))), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(C_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TtI(TRA)), ax))
  } else return(.Call(Cpp_TRAl,x,.Call(C_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TtI(TRA)))
}
