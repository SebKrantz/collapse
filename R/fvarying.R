
# TODO: Test!! and make fpdim !!

fvarying <- function(x, ...) UseMethod("fvarying") # , x
fvarying.default <- function(x, g = NULL, TRA = NULL, any_group = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fvarying,x,0L,0L,any_group)) else if(is.atomic(g)) {
      if(use.g.names && !any_group) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fvarying,x,length(lev),g,any_group), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fvarying,x,fnlevels(g),g,any_group)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fvarying,x,attr(g,"N.groups"),g,any_group))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names && !any_group) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !any_group) return(`names<-`(.Call(Cpp_fvarying,x,g[[1L]],g[[2L]],any_group), group_names.GRP(g))) else
        return(.Call(Cpp_fvarying,x,g[[1L]],g[[2L]],any_group))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fvarying,x,0L,0L,any_group),0L,TRAtoInt(TRA))) else if(is.atomic(g)) {
      if(is.nmfactor(g)) ng <- fnlevels(g) else {
        g <- qG(g, na.exclude = FALSE)
        ng <- attr(g,"N.groups")
      }
      return(.Call(Cpp_TRA,x,.Call(Cpp_fvarying,x,ng,g,any_group),if(any_group) 0L else g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fvarying,x,g[[1L]],g[[2L]],any_group),if(any_group) 0L else g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fvarying.matrix <- function(x, g = NULL, TRA = NULL, any_group = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fvaryingm,x,0L,0L,any_group,drop)) else if(is.atomic(g)) {
      if(use.g.names && !any_group) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fvaryingm,x,length(lev),g,any_group,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fvaryingm,x,fnlevels(g),g,any_group,drop)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fvaryingm,x,attr(g,"N.groups"),g,any_group,drop))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names && !any_group) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !any_group) return(`dimnames<-`(.Call(Cpp_fvaryingm,x,g[[1L]],g[[2L]],any_group,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fvaryingm,x,g[[1L]],g[[2L]],any_group,drop))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fvaryingm,x,0L,0L,any_group,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) ng <- fnlevels(g) else {
        g <- qG(g, na.exclude = FALSE)
        ng <- attr(g,"N.groups")
      }
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fvaryingm,x,ng,g,any_group,drop),if(any_group) 0L else g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fvaryingm,x,g[[1L]],g[[2L]],any_group,drop),if(any_group) 0L else g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fvarying.data.frame <- function(x, g = NULL, TRA = NULL, any_group = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fvaryingl,x,0L,0L,any_group,drop)) else if(is.atomic(g)) {
      if(use.g.names && !any_group && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fvaryingl,x,length(lev),g,any_group,FALSE), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fvaryingl,x,fnlevels(g),g,any_group,drop)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fvaryingl,x,attr(g,"N.groups"),g,any_group,drop))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names && !any_group) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !any_group && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fvaryingl,x,g[[1L]],g[[2L]],any_group,FALSE), groups)) else
          return(.Call(Cpp_fvaryingl,x,g[[1L]],g[[2L]],any_group,drop))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fvaryingl,x,0L,0L,any_group,TRUE),0L,TRAtoInt(TRA))) else if(is.atomic(g)) {
      if(is.nmfactor(g)) ng <- fnlevels(g) else {
        g <- qG(g, na.exclude = FALSE)
        ng <- attr(g,"N.groups")
      }
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fvaryingl,x,ng,g,any_group,drop),if(any_group) 0L else g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fvaryingl,x,g[[1L]],g[[2L]],any_group,drop),if(any_group) 0L else g[[2L]],TRAtoInt(TRA)))
    }
  }
}
# Make better version ???
fvarying.grouped_df <- function(x, TRA = NULL, any_group = TRUE, use.g.names = FALSE, keep.group_vars = !any_group, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
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
      ax[["row.names"]] <- if(use.g.names && !any_group) group_names.GRP(g) else if(!any_group) .set_row_names(g[[1L]]) else 1L
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], nam[-gn])
          return(setAttributes(c(g[[4L]],.Call(Cpp_fvaryingl,x[-gn],g[[1L]],g[[2L]],any_group,FALSE)), ax))
        } else {
          ax[["names"]] <- nam[-gn]
          return(setAttributes(.Call(Cpp_fvaryingl,x[-gn],g[[1L]],g[[2L]],any_group,FALSE), ax))
        }
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(Cpp_fvaryingl,x,g[[1L]],g[[2L]],any_group,FALSE)), ax))
      } else return(setAttributes(.Call(Cpp_fvaryingl,x,g[[1L]],g[[2L]],any_group,FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fvaryingl,x[-gn],g[[1L]],g[[2L]],any_group,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- nam[-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fvaryingl,x[-gn],g[[1L]],g[[2L]],any_group,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fvaryingl,x,g[[1L]],g[[2L]],any_group,FALSE),g[[2L]],TRAtoInt(TRA)))
}

