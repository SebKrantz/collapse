# library(Rcpp)
# sourceCpp('src/fnobs.cpp')
# sourceCpp('src/fnobsa.cpp')
# sourceCpp('src/fnobsl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

fNobs <- function(x, ...) { # g = NULL, TRA = NULL, use.g.names = TRUE, drop = TRUE, keep.group_vars = TRUE,
  UseMethod("fNobs", x)
}
fNobs.default <- function(x, g = NULL, TRA = NULL, use.g.names = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNobs,x,0L,0L)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fNobs,x,length(lev),g), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fNobs,x,fnlevels(g),g)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fNobs,x,attr(g,"N.groups"),g))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fNobs,x,g[[1L]],g[[2L]]), group_names.GRP(g))) else
        return(.Call(Cpp_fNobs,x,g[[1L]],g[[2L]]))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fNobs,x,0L,0L),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fNobs,x,fnlevels(g),g),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fNobs,x,attr(g,"N.groups"),g),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fNobs,x,g[[1L]],g[[2L]]),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fNobs.matrix <- function(x, g = NULL, TRA = NULL, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNobsm,x,0L,0L,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fNobsm,x,length(lev),g,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fNobsm,x,fnlevels(g),g,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fNobsm,x,attr(g,"N.groups"),g,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fNobsm,x,g[[1L]],g[[2L]],FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fNobsm,x,g[[1L]],g[[2L]],FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fNobsm,x,0L,0L,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fNobsm,x,fnlevels(g),g,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fNobsm,x,attr(g,"N.groups"),g,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fNobsm,x,g[[1L]],g[[2L]],FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fNobs.data.frame <- function(x, g = NULL, TRA = NULL, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNobsl,x,0L,0L,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fNobsl,x,length(lev),g,FALSE), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fNobsl,x,fnlevels(g),g,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fNobsl,x,attr(g,"N.groups"),g,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE), groups)) else
          return(.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fNobsl,x,0L,0L,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fNobsl,x,fnlevels(g),g,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fNobsl,x,attr(g,"N.groups"),g,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fNobs.grouped_df <- function(x, TRA = NULL, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  gn <- which(names(x) %in% g[[5L]])
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
          ax[["names"]] <- c(g[[5L]], ax[["names"]][-gn])
          return(setAttributes(c(g[[4L]],.Call(Cpp_fNobsl,x[-gn],g[[1L]],g[[2L]],FALSE)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(.Call(Cpp_fNobsl,x[-gn],g[[1L]],g[[2L]],FALSE), ax))
        }
      } else return(setAttributes(.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fNobsl,x[-gn],g[[1L]],g[[2L]],FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fNobsl,x[-gn],g[[1L]],g[[2L]],FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE),g[[2L]],TRAtoInt(TRA)))
}
