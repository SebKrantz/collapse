
# For foundational changes to this code see fsum.R !!

fmin <- function(x, ...) UseMethod("fmin") # , x

fmin.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fminmax,x,0L,0L,na.rm,1L)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fminmax,x,length(lev),g,na.rm,1L), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fminmax,x,fnlevels(g),g,na.rm,1L)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fminmax,x,attr(g,"N.groups"),g,na.rm,1L))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fminmax,x,g[[1L]],g[[2L]],na.rm,1L), group_names.GRP(g))) else
        return(.Call(Cpp_fminmax,x,g[[1L]],g[[2L]],na.rm,1L))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fminmax,x,0L,0L,na.rm,1L),0L,TtI(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fminmax,x,fnlevels(g),g,na.rm,1L),g,TtI(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fminmax,x,attr(g,"N.groups"),g,na.rm,1L),g,TtI(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fminmax,x,g[[1L]],g[[2L]],na.rm,1L),g[[2L]],TtI(TRA)))
    }
  }
}
fmin.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fminmaxm,x,0L,0L,na.rm,drop,1L)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fminmaxm,x,length(lev),g,na.rm,FALSE,1L), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fminmaxm,x,fnlevels(g),g,na.rm,FALSE,1L)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fminmaxm,x,attr(g,"N.groups"),g,na.rm,FALSE,1L))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fminmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE,1L), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fminmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE,1L))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fminmaxm,x,0L,0L,na.rm,TRUE,1L),0L,TtI(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fminmaxm,x,fnlevels(g),g,na.rm,FALSE,1L),g,TtI(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fminmaxm,x,attr(g,"N.groups"),g,na.rm,FALSE,1L),g,TtI(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fminmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE,1L),g[[2L]],TtI(TRA)))
    }
  }
}
fmin.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fminmaxl,x,0L,0L,na.rm,drop,1L)) else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fminmaxl,x,length(lev),g,na.rm,FALSE,1L), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fminmaxl,x,fnlevels(g),g,na.rm,FALSE,1L)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fminmaxl,x,attr(g,"N.groups"),g,na.rm,FALSE,1L))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,1L), groups)) else
          return(.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,1L))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fminmaxl,x,0L,0L,na.rm,TRUE,1L),0L,TtI(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fminmaxl,x,fnlevels(g),g,na.rm,FALSE,1L),g,TtI(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fminmaxl,x,attr(g,"N.groups"),g,na.rm,FALSE,1L),g,TtI(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,1L),g[[2L]],TtI(TRA)))
    }
  }
}
fmin.list <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
  fmin.data.frame(x, g, TRA, na.rm, use.g.names, drop, ...)
fmin.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x)
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  gl <- length(gn) > 0L
  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], nam[-gn])
          return(setAttributes(c(g[[4L]],.Call(Cpp_fminmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE,1L)), ax))
        } else {
          ax[["names"]] <- nam[-gn]
          return(setAttributes(.Call(Cpp_fminmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE,1L), ax))
        }
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,1L)), ax))
      } else return(setAttributes(.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,1L), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fminmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE,1L),g[[2L]],TtI(TRA))), ax))
    } else {
      ax[["names"]] <- nam[-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fminmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE,1L),g[[2L]],TtI(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,1L),g[[2L]],TtI(TRA)))
}


fmax <- function(x, ...) UseMethod("fmax") # , x

fmax.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fminmax,x,0L,0L,na.rm,2L)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fminmax,x,length(lev),g,na.rm,2L), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fminmax,x,fnlevels(g),g,na.rm,2L)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fminmax,x,attr(g,"N.groups"),g,na.rm,2L))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fminmax,x,g[[1L]],g[[2L]],na.rm,2L), group_names.GRP(g))) else
        return(.Call(Cpp_fminmax,x,g[[1L]],g[[2L]],na.rm,2L))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fminmax,x,0L,0L,na.rm,2L),0L,TtI(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fminmax,x,fnlevels(g),g,na.rm,2L),g,TtI(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fminmax,x,attr(g,"N.groups"),g,na.rm,2L),g,TtI(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fminmax,x,g[[1L]],g[[2L]],na.rm,2L),g[[2L]],TtI(TRA)))
    }
  }
}
fmax.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fminmaxm,x,0L,0L,na.rm,drop,2L)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fminmaxm,x,length(lev),g,na.rm,FALSE,2L), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fminmaxm,x,fnlevels(g),g,na.rm,FALSE,2L)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fminmaxm,x,attr(g,"N.groups"),g,na.rm,FALSE,2L))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fminmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE,2L), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fminmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE,2L))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fminmaxm,x,0L,0L,na.rm,TRUE,2L),0L,TtI(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fminmaxm,x,fnlevels(g),g,na.rm,FALSE,2L),g,TtI(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fminmaxm,x,attr(g,"N.groups"),g,na.rm,FALSE,2L),g,TtI(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fminmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE,2L),g[[2L]],TtI(TRA)))
    }
  }
}
fmax.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fminmaxl,x,0L,0L,na.rm,drop,2L)) else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fminmaxl,x,length(lev),g,na.rm,FALSE,2L), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fminmaxl,x,fnlevels(g),g,na.rm,FALSE,2L)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fminmaxl,x,attr(g,"N.groups"),g,na.rm,FALSE,2L))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,2L), groups)) else
          return(.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,2L))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fminmaxl,x,0L,0L,na.rm,TRUE,2L),0L,TtI(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fminmaxl,x,fnlevels(g),g,na.rm,FALSE,2L),g,TtI(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fminmaxl,x,attr(g,"N.groups"),g,na.rm,FALSE,2L),g,TtI(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,2L),g[[2L]],TtI(TRA)))
    }
  }
}
fmax.list <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
  fmax.data.frame(x, g, TRA, na.rm, use.g.names, drop, ...)
fmax.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x)
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  gl <- length(gn) > 0L
  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], nam[-gn])
          return(setAttributes(c(g[[4L]],.Call(Cpp_fminmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE,2L)), ax))
        } else {
          ax[["names"]] <- nam[-gn]
          return(setAttributes(.Call(Cpp_fminmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE,2L), ax))
        }
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,2L)), ax))
      } else return(setAttributes(.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,2L), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fminmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE,2L),g[[2L]],TtI(TRA))), ax))
    } else {
      ax[["names"]] <- nam[-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fminmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE,2L),g[[2L]],TtI(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fminmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE,2L),g[[2L]],TtI(TRA)))
}

