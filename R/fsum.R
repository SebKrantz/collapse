# library(Rcpp)
# sourceCpp('src/fsum.cpp')
# sourceCpp('src/fsuma.cpp')
# sourceCpp('src/fsuml.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')


fsum <- function(x, ...) { # g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_vars = TRUE
  UseMethod("fsum", x)
}
fsum.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fsum,x,0L,0L,na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fsum,x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fsum,x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(.Call(Cpp_fsum,x,attr(g,"N.groups"),g,na.rm))
        }
      }
      # } else if (.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names) {
      #     lev <- attr(g, "levels")
      #     return(`names<-`(.Call(Cpp_fsum,x,length(lev),g,na.rm), lev))
      #   } else return(.Call(Cpp_fsum,x,fnlevels(g),g,na.rm))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fsum,x,g[[1L]],g[[2L]],na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_fsum,x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
    #   g <- interaction(g)
    #   return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsum.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fsumm,x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fsumm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fsumm,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g)
          return(.Call(Cpp_fsumm,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
      # } else if (.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names) {
      #     lev <- attr(g, "levels")
      #     return(`dimnames<-`(.Call(Cpp_fsumm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      #   } else return(.Call(Cpp_fsumm,x,fnlevels(g),g,na.rm,FALSE))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fsumm,x,g[[1L]],g[[2L]],na.rm,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fsumm,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsum.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fsuml,x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fsuml,x,length(lev),g,na.rm,FALSE), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fsuml,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g)
          return(.Call(Cpp_fsuml,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names && !inherits(x, "data.table")) {
      #     lev <- attr(g, "levels")
      #     return(setRow.names(.Call(Cpp_fsuml,x,length(lev),g,na.rm,FALSE), lev))
      #   } else return(.Call(Cpp_fsuml,x,fnlevels(g),g,na.rm,FALSE))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],na.rm,FALSE), groups)) else
          return(.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsum.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) { # drop grouping columns argument ??
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  gn <- which(names(x) %in% g[[5L]]) # faster than na.rm(match(names(g[[4L]]), names(x)))
  nTRAl <- is.null(TRA)
  gl <- length(gn) > 0L
  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL # remove attributes is more memory efficient than unclass !!
    # if(gl) ax[["names"]] <- if(keep.group_vars) c(g[[5L]], ax[["names"]][-gn]) else ax[["names"]][-gn]
    if(nTRAl) {
      ax[["groups"]] <- NULL # best way ??? -> ye, not slower than following (tested GGDC) ngp <- name(ax) != "groups"
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], ax[["names"]][-gn])
          return(setAttributes(c(g[[4L]],.Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(.Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE), ax))
        }
      } else return(setAttributes(.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],na.rm,FALSE), ax))
    } else if(keep.group_vars) {
         ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
         return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
         ax[["names"]] <- ax[["names"]][-gn]
         return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
}
