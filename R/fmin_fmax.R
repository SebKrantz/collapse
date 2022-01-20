
# For foundational changes to this code see fsum.R !!

fmin <- function(x, ...) UseMethod("fmin") # , x

fmin.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fmin.matrix(x, g, TRA, na.rm, use.g.names, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
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
  if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(C_fmin,x,0L,0L,na.rm),0L,TtI(TRA)))
  g <- G_guo(g)
  .Call(Cpp_TRA,x,.Call(C_fmin,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA))
}

fmin.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
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
  if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(C_fminm,x,0L,0L,na.rm,TRUE),0L,TtI(TRA)))
  g <- G_guo(g)
  .Call(Cpp_TRAm,x,.Call(C_fminm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA))
}

fmin.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
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
  if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(C_fminl,x,0L,0L,na.rm,TRUE),0L,TtI(TRA)))
  g <- G_guo(g)
  .Call(Cpp_TRAl,x,.Call(C_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA))
}

fmin.list <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
  fmin.data.frame(x, g, TRA, na.rm, use.g.names, drop, ...)

fmin.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
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
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(C_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA))), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(C_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA)), ax))
  } else return(.Call(Cpp_TRAl,x,.Call(C_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA)))
}


fmax <- function(x, ...) UseMethod("fmax") # , x

fmax.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fmax.matrix(x, g, TRA, na.rm, use.g.names, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
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
  if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(C_fmax,x,0L,0L,na.rm),0L,TtI(TRA)))
  g <- G_guo(g)
  .Call(Cpp_TRA,x,.Call(C_fmax,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA))
}

fmax.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
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
  if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(C_fmaxm,x,0L,0L,na.rm,TRUE),0L,TtI(TRA)))
  g <- G_guo(g)
  .Call(Cpp_TRAm,x,.Call(C_fmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA))
}

fmax.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
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
  if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(C_fmaxl,x,0L,0L,na.rm,TRUE),0L,TtI(TRA)))
  g <- G_guo(g)
  .Call(Cpp_TRAl,x,.Call(C_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA))
}

fmax.list <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
  fmax.data.frame(x, g, TRA, na.rm, use.g.names, drop, ...)

fmax.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
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
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(C_fmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA))), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(C_fmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA)), ax))
  } else return(.Call(Cpp_TRAl,x,.Call(C_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA)))
}

