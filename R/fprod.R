
# For foundational changes to this code see fsum.R

fprod <- function(x, ...) UseMethod("fprod") # , x

fprod.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fprod.matrix(x, g, w, TRA, na.rm, use.g.names, ...))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fprod,x,0L,0L,w,na.rm))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(C_fprod,x,length(lev),g,w,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fprod,x,fnlevels(g),g,w,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fprod,x,attr(g,"N.groups"),g,w,na.rm))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(C_fprod,x,g[[1L]],g[[2L]],w,na.rm), GRPnames(g)))
    return(.Call(C_fprod,x,g[[1L]],g[[2L]],w,na.rm))
  }
  if(is.null(g)) return(TRAC(x,.Call(C_fprod,x,0L,0L,w,na.rm),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(C_fprod,x,g[[1L]],g[[2L]],w,na.rm),g[[2L]],TRA, ...)
}

fprod.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, drop = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fprodm,x,0L,0L,w,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(C_fprodm,x,length(lev),g,w,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(C_fprodm,x,fnlevels(g),g,w,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fprodm,x,attr(g,"N.groups"),g,w,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(C_fprodm,x,g[[1L]],g[[2L]],w,na.rm,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(C_fprodm,x,g[[1L]],g[[2L]],w,na.rm,FALSE))
  }
  if(is.null(g)) return(TRAmC(x,.Call(C_fprodm,x,0L,0L,w,na.rm,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(C_fprodm,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TRA, ...)
}

fprod.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, drop = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fprodl,x,0L,0L,w,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(C_fprodl,x,length(lev),g,w,na.rm,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fprodl,x,fnlevels(g),g,w,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fprodl,x,attr(g,"N.groups"),g,w,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(C_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE), groups))
    return(.Call(C_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE))
  }
  if(is.null(g)) return(TRAlC(x,.Call(C_fprodl,x,0L,0L,w,na.rm,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(C_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TRA, ...)
}

fprod.list <- function(x, ...) fprod.data.frame(x, ...)

fprod.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = FALSE,
                            keep.group_vars = TRUE, keep.w = TRUE, stub = .op[["stub"]], ...) {
  g <- GRP.grouped_df(x, call = FALSE)
  if(is.null(g[[4L]])) keep.group_vars <- FALSE
  wsym <- substitute(w)
  nam <- attr(x, "names")
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  prodw <- NULL

  if(!is.null(wsym)) {
    w <- eval(wsym, x, parent.frame())
    if(length(wn <- which(nam %in% all.vars(wsym)))) {
      if(any(gn %in% wn)) stop("Weights coincide with grouping variables!")
      gn <- c(gn, wn)
      if(keep.w) {
        if(nTRAl) prodw <- `names<-`(list(.Call(C_fprod,w,g[[1L]],g[[2L]],NULL,na.rm)), do_stub(stub, if(length(wsym) == 1L) as.character(wsym) else deparse(wsym), "prod.")) else if(keep.group_vars)
          gn2 <- gn else prodw <- gn2 <- wn
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
          ax[["names"]] <- c(g[[5L]], names(prodw), nam[-gn])
          return(setAttributes(c(g[[4L]], prodw, .Call(C_fprodl,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
        }
        ax[["names"]] <- c(names(prodw), nam[-gn])
        return(setAttributes(c(prodw, .Call(C_fprodl,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(C_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
      } else return(setAttributes(.Call(C_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE), ax))
    } else if(keep.group_vars || (keep.w && length(prodw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],TRAlC(x[-gn],.Call(C_fprodl,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_fprodl,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_fprodl,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TRA, ...))
}

