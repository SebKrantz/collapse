# library(Rcpp)
# sourceCpp('src/fmin.cpp')
# sourceCpp('src/fmina.cpp')
# sourceCpp('src/fminl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

fmin <- function(x, ...) UseMethod("fmin") # , x

fmin.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fmin,x,0L,0L,na.rm)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fmin,x,length(lev),g,na.rm), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fmin,x,fnlevels(g),g,na.rm)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fmin,x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fmin,x,g[[1L]],g[[2L]],na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_fmin,x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmin,x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmin,x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fmin,x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fmin,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmin.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fminm,x,0L,0L,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fminm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fminm,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fminm,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fminm,x,g[[1L]],g[[2L]],na.rm,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fminm,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fminm,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fminm,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fminm,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fminm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmin.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fminl,x,0L,0L,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fminl,x,length(lev),g,na.rm,FALSE), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fminl,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fminl,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE), groups)) else
          return(.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fminl,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fminl,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fminl,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmin.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
      ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], nam[-gn])
          return(setAttributes(c(g[[4L]],.Call(Cpp_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE)), ax))
        } else {
          ax[["names"]] <- nam[-gn]
          return(setAttributes(.Call(Cpp_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE), ax))
        }
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE)), ax))
      } else return(setAttributes(.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- nam[-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
}


