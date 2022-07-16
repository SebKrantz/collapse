
# TODO: w.type - Implement reliability weights?
# Note: for principal innovations of this code see fsum.R

fsd <- function(x, ...) UseMethod("fsd") # , x

fsd.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fsd.matrix(x, g, w, TRA, na.rm, use.g.names, stable.algo = stable.algo, ...))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fvarsd,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fvarsd,x,length(lev),g,NULL,w,na.rm,stable.algo,TRUE), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fvarsd,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,TRUE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fvarsd,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,TRUE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE), GRPnames(g)))
    return(.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE))
  }
  if(is.null(g)) return(TRAC(x,.Call(Cpp_fvarsd,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE),g[[2L]],TRA, ...)
}

fsd.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fvarsdm,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fvarsdm,x,length(lev),g,NULL,w,na.rm,stable.algo,TRUE,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fvarsdm,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,TRUE,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fvarsdm,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,TRUE,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE))
  }
  if(is.null(g)) return(TRAmC(x,.Call(Cpp_fvarsdm,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE),g[[2L]],TRA, ...)
}

fsd.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fvarsdl,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,drop))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(Cpp_fvarsdl,x,length(lev),g,NULL,w,na.rm,stable.algo,TRUE,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fvarsdl,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,TRUE,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fvarsdl,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,TRUE,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE), groups))
    return(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE))
  }
  if(is.null(g)) return(TRAlC(x,.Call(Cpp_fvarsdl,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE),g[[2L]],TRA, ...)
}

fsd.list <- function(x, ...) fsd.data.frame(x, ...)

fsd.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                             keep.group_vars = TRUE, keep.w = TRUE, stable.algo = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]], sumw, .Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE)), ax))
      } else return(setAttributes(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE),g[[2L]],TRA, ...))
}



fvar <- function(x, ...) UseMethod("fvar") # , x

fvar.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fvar.matrix(x, g, w, TRA, na.rm, use.g.names, stable.algo = stable.algo, ...))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fvarsd,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fvarsd,x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fvarsd,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fvarsd,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE), GRPnames(g)))
    return(.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE))
  }
  if(is.null(g)) return(TRAC(x,.Call(Cpp_fvarsd,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE),g[[2L]],TRA, ...)
}

fvar.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fvarsdm,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fvarsdm,x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fvarsdm,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fvarsdm,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE))
  }
  if(is.null(g)) return(TRAmC(x,.Call(Cpp_fvarsdm,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE),g[[2L]],TRA, ...)
}

fvar.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(Cpp_fvarsdl,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,drop))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(Cpp_fvarsdl,x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fvarsdl,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fvarsdl,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE), groups))
    return(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE))
  }
  if(is.null(g)) return(TRAlC(x,.Call(Cpp_fvarsdl,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE),g[[2L]],TRA, ...)
}

fvar.list <- function(x, ...) fvar.data.frame(x, ...)

fvar.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                           keep.group_vars = TRUE, keep.w = TRUE, stable.algo = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]], sumw, .Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE)), ax))
        }
        ax[["names"]] <- c(names(sumw), nam[-gn])
        return(setAttributes(c(sumw, .Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE)), ax))
      } else return(setAttributes(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE),g[[2L]],TRA, ...))
}


