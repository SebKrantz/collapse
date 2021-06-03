
# For foundational changes to this code see fsum.R

fprod <- function(x, ...) UseMethod("fprod") # , x

fprod.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fprod.matrix(x, g, w, TRA, na.rm, use.g.names, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fprod,x,0L,0L,w,na.rm))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fprod,x,length(lev),g,w,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fprod,x,fnlevels(g),g,w,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fprod,x,attr(g,"N.groups"),g,w,na.rm))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(Cpp_fprod,x,g[[1L]],g[[2L]],w,na.rm), GRPnames(g)))
    return(.Call(Cpp_fprod,x,g[[1L]],g[[2L]],w,na.rm))
  }
  if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fprod,x,0L,0L,w,na.rm),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fprod,x,fnlevels(g),g,w,na.rm),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRA,x,.Call(Cpp_fprod,x,attr(g,"N.groups"),g,w,na.rm),g,TtI(TRA)))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRA,x,.Call(Cpp_fprod,x,g[[1L]],g[[2L]],w,na.rm),g[[2L]],TtI(TRA))
}

fprod.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fprodm,x,0L,0L,w,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fprodm,x,length(lev),g,w,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fprodm,x,fnlevels(g),g,w,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fprodm,x,attr(g,"N.groups"),g,w,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(Cpp_fprodm,x,g[[1L]],g[[2L]],w,na.rm,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(Cpp_fprodm,x,g[[1L]],g[[2L]],w,na.rm,FALSE))
  }
  if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fprodm,x,0L,0L,w,na.rm,TRUE),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fprodm,x,fnlevels(g),g,w,na.rm,FALSE),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRAm,x,.Call(Cpp_fprodm,x,attr(g,"N.groups"),g,w,na.rm,FALSE),g,TtI(TRA)))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRAm,x,.Call(Cpp_fprodm,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TtI(TRA))
}

fprod.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fprodl,x,0L,0L,w,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(Cpp_fprodl,x,length(lev),g,w,na.rm,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fprodl,x,fnlevels(g),g,w,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fprodl,x,attr(g,"N.groups"),g,w,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(Cpp_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE), groups))
    return(.Call(Cpp_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE))
  }
  if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fprodl,x,0L,0L,w,na.rm,TRUE),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fprodl,x,fnlevels(g),g,w,na.rm,FALSE),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRAl,x,.Call(Cpp_fprodl,x,attr(g,"N.groups"),g,w,na.rm,FALSE),g,TtI(TRA)))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRAl,x,.Call(Cpp_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TtI(TRA))
}

fprod.list <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
  fprod.data.frame(x, g, w, TRA, na.rm, use.g.names, drop, ...)

fprod.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                            keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  wsym <- l1orn(as.character(substitute(w)), NULL)
  nam <- attr(x, "names")
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  prodw <- NULL

  if(length(wsym) && length(wn <- which(wsym == nam))) {
    w <- .subset2(x, wn)
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    gn <- c(gn, wn)
    if(keep.w) {
      if(nTRAl) prodw <- `names<-`(list(fprodCpp(w,g[[1L]],g[[2L]],NULL,na.rm)), paste0("prod.", wsym)) else if(keep.group_vars)
        gn2 <- gn else prodw <- gn2 <- wn
    }
  }

  gl <- length(gn) > 0L

  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- fsetdiff(ax[["class"]], c("GRP_df", "grouped_df"))
      ax[["row.names"]] <- if(use.g.names) GRPnames(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], names(prodw), nam[-gn])
          return(setAttributes(c(g[[4L]], prodw, .Call(Cpp_fprodl,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
        }
        ax[["names"]] <- c(names(prodw), nam[-gn])
        return(setAttributes(c(prodw, .Call(Cpp_fprodl,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(Cpp_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
      } else return(setAttributes(.Call(Cpp_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE), ax))
    } else if(keep.group_vars || (keep.w && length(prodw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fprodl,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TtI(TRA))), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fprodl,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TtI(TRA)), ax))
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TtI(TRA)))
}

