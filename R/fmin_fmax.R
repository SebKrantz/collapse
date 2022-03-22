
# For foundational changes to this code see fsum.R !!

fmin <- function(x, ...) UseMethod("fmin") # , x

fmin.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fmin.matrix(x, g, TRA, na.rm, use.g.names, ...))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fmin,x,0L,0L,na.rm))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(C_fmin,x,length(lev),g,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fmin,x,fnlevels(g),g,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fmin,x,attr(g,"N.groups"),g,na.rm))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(C_fmin,x,g[[1L]],g[[2L]],na.rm), GRPnames(g)))
    return(.Call(C_fmin,x,g[[1L]],g[[2L]],na.rm))
  }
  if(is.null(g)) return(TRAC(x,.Call(C_fmin,x,0L,0L,na.rm),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(C_fmin,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRA, ...)
}

fmin.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fminm,x,0L,0L,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(C_fminm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(C_fminm,x,fnlevels(g),g,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fminm,x,attr(g,"N.groups"),g,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(C_fminm,x,g[[1L]],g[[2L]],na.rm,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(C_fminm,x,g[[1L]],g[[2L]],na.rm,FALSE))
  }
  if(is.null(g)) return(TRAmC(x,.Call(C_fminm,x,0L,0L,na.rm,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(C_fminm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRA, ...)
}

fmin.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fminl,x,0L,0L,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(C_fminl,x,length(lev),g,na.rm,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fminl,x,fnlevels(g),g,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fminl,x,attr(g,"N.groups"),g,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(C_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE), groups))
    return(.Call(C_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE))
  }
  if(is.null(g)) return(TRAlC(x,.Call(C_fminl,x,0L,0L,na.rm,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(C_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRA, ...)
}

fmin.list <- function(x, ...) fmin.data.frame(x, ...)

fmin.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  g <- GRP.grouped_df(x, call = FALSE)
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
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
          ax[["names"]] <- c(g[[5L]], nam[-gn])
          return(setAttributes(c(g[[4L]],.Call(C_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE)), ax))
        }
        ax[["names"]] <- nam[-gn]
        return(setAttributes(.Call(C_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(C_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE)), ax))
      } else return(setAttributes(.Call(C_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],TRAlC(x[-gn],.Call(C_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRA, ...))
}


fmax <- function(x, ...) UseMethod("fmax") # , x

fmax.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fmax.matrix(x, g, TRA, na.rm, use.g.names, ...))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fmax,x,0L,0L,na.rm))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(C_fmax,x,length(lev),g,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fmax,x,fnlevels(g),g,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fmax,x,attr(g,"N.groups"),g,na.rm))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(C_fmax,x,g[[1L]],g[[2L]],na.rm), GRPnames(g)))
    return(.Call(C_fmax,x,g[[1L]],g[[2L]],na.rm))
  }
  if(is.null(g)) return(TRAC(x,.Call(C_fmax,x,0L,0L,na.rm),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(C_fmax,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRA, ...)
}

fmax.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fmaxm,x,0L,0L,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(C_fmaxm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(C_fmaxm,x,fnlevels(g),g,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fmaxm,x,attr(g,"N.groups"),g,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(C_fmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(C_fmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE))
  }
  if(is.null(g)) return(TRAmC(x,.Call(C_fmaxm,x,0L,0L,na.rm,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(C_fmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRA, ...)
}

fmax.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fmaxl,x,0L,0L,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(C_fmaxl,x,length(lev),g,na.rm,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fmaxl,x,fnlevels(g),g,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fmaxl,x,attr(g,"N.groups"),g,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(C_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE), groups))
    return(.Call(C_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE))
  }
  if(is.null(g)) return(TRAlC(x,.Call(C_fmaxl,x,0L,0L,na.rm,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(C_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRA, ...)
}

fmax.list <- function(x, ...) fmax.data.frame(x, ...)

fmax.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  g <- GRP.grouped_df(x, call = FALSE)
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
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
          ax[["names"]] <- c(g[[5L]], nam[-gn])
          return(setAttributes(c(g[[4L]],.Call(C_fmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE)), ax))
        }
        ax[["names"]] <- nam[-gn]
        return(setAttributes(.Call(C_fmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(C_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE)), ax))
      } else return(setAttributes(.Call(C_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],TRAlC(x[-gn],.Call(C_fmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_fmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRA, ...))
}

