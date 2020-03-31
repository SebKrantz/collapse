# library(Rcpp)
# sourceCpp('src/fsum.cpp')
# sourceCpp('src/fsuma.cpp')
# sourceCpp('src/fsuml.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')


fsum <- function(x, ...) {
  UseMethod("fsum", x)
}
fsum.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fsum,x,0L,0L,w,na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fsum,x,length(lev),g,w,na.rm), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fsum,x,fnlevels(g),g,w,na.rm)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fsum,x,attr(g,"N.groups"),g,w,na.rm))
        }
      }
      # } else if (.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names) {
      #     lev <- attr(g, "levels")
      #     return(`names<-`(.Call(Cpp_fsum,x,length(lev),g,w,na.rm), lev))
      #   } else return(.Call(Cpp_fsum,x,fnlevels(g),g,w,na.rm))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fsum,x,g[[1L]],g[[2L]],w,na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_fsum,x,g[[1L]],g[[2L]],w,na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,0L,0L,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,fnlevels(g),g,w,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,attr(g,"N.groups"),g,w,na.rm),g,TRAtoInt(TRA)))
      }
    # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
    #   g <- interaction(g)
    #   return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,fnlevels(g),g,w,na.rm),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,g[[1L]],g[[2L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsum.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fsumm,x,0L,0L,w,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fsumm,x,length(lev),g,w,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fsumm,x,fnlevels(g),g,w,na.rm,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fsumm,x,attr(g,"N.groups"),g,w,na.rm,FALSE))
        }
      }
      # } else if (.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names) {
      #     lev <- attr(g, "levels")
      #     return(`dimnames<-`(.Call(Cpp_fsumm,x,length(lev),g,w,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      #   } else return(.Call(Cpp_fsumm,x,fnlevels(g),g,w,na.rm,FALSE))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fsumm,x,g[[1L]],g[[2L]],w,na.rm,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fsumm,x,g[[1L]],g[[2L]],w,na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,0L,0L,w,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,fnlevels(g),g,w,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,attr(g,"N.groups"),g,w,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,fnlevels(g),g,w,na.rm,FALSE),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsum.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fsuml,x,0L,0L,w,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fsuml,x,length(lev),g,w,na.rm,FALSE), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fsuml,x,fnlevels(g),g,w,na.rm,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fsuml,x,attr(g,"N.groups"),g,w,na.rm,FALSE))
        }
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names && !inherits(x, "data.table")) {
      #     lev <- attr(g, "levels")
      #     return(setRow.names(.Call(Cpp_fsuml,x,length(lev),g,w,na.rm,FALSE), lev))
      #   } else return(.Call(Cpp_fsuml,x,fnlevels(g),g,w,na.rm,FALSE))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE), groups)) else
          return(.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,0L,0L,w,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,fnlevels(g),g,w,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,attr(g,"N.groups"),g,w,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,fnlevels(g),g,w,na.rm,FALSE),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsum.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                            keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- attr(x, "names")
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  sumw <- NULL

  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- unclass(x)[[wn]] # faster using unclass??
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    onlyw <- !length(gn)
    gn <- c(gn, wn)
    if(keep.w) {
      if(nTRAl) sumw <- `names<-`(list(fsumCpp(w,g[[1L]],g[[2L]],NULL,na.rm)), paste0("sum.", wsym)) else if(keep.group_vars)
        gn2 <- gn else sumw <- gn2 <- wn
    }
  } else onlyw <- FALSE

  gl <- length(gn) > 0L # necessary here, not before !!!

  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars && !onlyw) {
          ax[["names"]] <- c(g[[5L]], names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(g[[4L]], sumw, .Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
        } else {
          ax[["names"]] <- c(names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(sumw, .Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE)), ax))
        }
      } else return(setAttributes(.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(ax[["names"]][gn2], ax[["names"]][-gn])
      return(setAttributes(c(x[gn2],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
}
