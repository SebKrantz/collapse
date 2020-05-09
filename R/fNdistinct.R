

# For foundational changes to this code see fsum.R
# Note: matrix method needs memory equal to size of the object, while data.frame method does not need any memory ?!

fNdistinct <- function(x, ...) UseMethod("fNdistinct") # , x

fNdistinct.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNdistinct,x,0L,0L,NULL,na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fNdistinct,x,length(lev),g,NULL,na.rm), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fNdistinct,x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fNdistinct,x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fNdistinct,x,g[[1L]],g[[2L]],g[[3L]],na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_fNdistinct,x,g[[1L]],g[[2L]],g[[3L]],na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fNdistinct,x,0L,0L,NULL,na.rm),0L,TtI(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fNdistinct,x,fnlevels(g),g,NULL,na.rm),g,TtI(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fNdistinct,x,attr(g,"N.groups"),g,NULL,na.rm),g,TtI(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fNdistinct,x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TtI(TRA)))
    }
  }
}
fNdistinct.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNdistinctm,x,0L,0L,NULL,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fNdistinctm,x,length(lev),g,NULL,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fNdistinctm,x,fnlevels(g),g,NULL,na.rm,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fNdistinctm,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fNdistinctm,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fNdistinctm,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fNdistinctm,x,0L,0L,NULL,na.rm,TRUE),0L,TtI(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fNdistinctm,x,fnlevels(g),g,NULL,na.rm,FALSE),g,TtI(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fNdistinctm,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE),g,TtI(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fNdistinctm,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TtI(TRA)))
    }
  }
}
fNdistinct.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNdistinctl,x,0L,0L,NULL,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fNdistinctl,x,length(lev),g,NULL,na.rm,FALSE), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fNdistinctl,x,fnlevels(g),g,NULL,na.rm,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fNdistinctl,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), groups)) else
          return(.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fNdistinctl,x,0L,0L,NULL,na.rm,TRUE),0L,TtI(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fNdistinctl,x,fnlevels(g),g,NULL,na.rm,FALSE),g,TtI(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fNdistinctl,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE),g,TtI(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TtI(TRA)))
    }
  }
}
fNdistinct.list <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
  fNdistinct.data.frame(x, g, TRA, na.rm, use.g.names, drop, ...)
fNdistinct.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]],.Call(Cpp_fNdistinctl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE)), ax))
        } else {
          ax[["names"]] <- nam[-gn]
          return(setAttributes(.Call(Cpp_fNdistinctl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), ax))
        }
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE)), ax))
      } else return(setAttributes(.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fNdistinctl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TtI(TRA))), ax))
    } else {
      ax[["names"]] <- nam[-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fNdistinctl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TtI(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TtI(TRA)))
}

