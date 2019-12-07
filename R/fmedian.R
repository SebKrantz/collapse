# library(Rcpp)
# sourceCpp('src/fmedian.cpp')
# sourceCpp('src/fmediana.cpp')
# sourceCpp('src/fmedianl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

fmedian <- function(x, ...) { # g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_vars = TRUE,
  UseMethod("fmedian", x)
}
fmedian.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fmedian,x,0L,0L,NULL,na.rm)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fmedian,x,length(lev),g,NULL,na.rm), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fmedian,x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fmedian,x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fmedian,x,g[[1L]],g[[2L]],g[[3L]],na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_fmedian,x,g[[1L]],g[[2L]],g[[3L]],na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmedian,x,0L,0L,NULL,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmedian,x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fmedian,x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fmedian,x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmedian.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fmedianm,x,0L,0L,NULL,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fmedianm,x,length(lev),g,NULL,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fmedianm,x,fnlevels(g),g,NULL,na.rm,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fmedianm,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fmedianm,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fmedianm,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fmedianm,x,0L,0L,NULL,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fmedianm,x,fnlevels(g),g,NULL,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fmedianm,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fmedianm,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmedian.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fmedianl,x,0L,0L,NULL,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fmedianl,x,length(lev),g,NULL,na.rm,FALSE), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fmedianl,x,fnlevels(g),g,NULL,na.rm,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fmedianl,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fmedianl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), groups)) else
          return(.Call(Cpp_fmedianl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fmedianl,x,0L,0L,NULL,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fmedianl,x,fnlevels(g),g,NULL,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fmedianl,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fmedianl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmedian.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]],.Call(Cpp_fmedianl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(.Call(Cpp_fmedianl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), ax))
        }
      } else return(setAttributes(.Call(Cpp_fmedianl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fmedianl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fmedianl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fmedianl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
}
