# library(Rcpp)
# sourceCpp('src/fsum.cpp')
# sourceCpp('src/fvarsd.cpp')
# sourceCpp('src/fvarsda.cpp')
# sourceCpp('src/fvarsdl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# w.type = "frequency"
# Note: for principal innovations of this code see fsum.R !!

fsd <- function(x, ...) { # g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, drop = TRUE, keep.group_vars = TRUE, keep.w = TRUE,
  UseMethod("fsd", x)
}
fsd.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fvarsd,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fvarsd,x,length(lev),g,NULL,w,na.rm,stable.algo,TRUE), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fvarsd,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,TRUE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fvarsd,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,TRUE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE), group_names.GRP(g))) else
        return(.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fvarsd,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fvarsd,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,TRUE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fvarsd,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,TRUE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsd.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fvarsdm,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fvarsdm,x,length(lev),g,NULL,w,na.rm,stable.algo,TRUE,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fvarsdm,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,TRUE,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fvarsdm,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,TRUE,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fvarsdm,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fvarsdm,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,TRUE,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fvarsdm,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,TRUE,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsd.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fvarsdl,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fvarsdl,x,length(lev),g,NULL,w,na.rm,stable.algo,TRUE,FALSE), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fvarsdl,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,TRUE,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fvarsdl,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,TRUE,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE), groups)) else
          return(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fvarsdl,x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fvarsdl,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,TRUE,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fvarsdl,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,TRUE,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsd.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                           keep.group_vars = TRUE, keep.w = TRUE, stable.algo = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w)) # Speed up somehow ??
  nam <- attr(x, "names")
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  sumw <- NULL

  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- unclass(x)[[wn]]
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    onlyw <- !length(gn)
    gn <- c(gn, wn)
    if(keep.w) {
     if(nTRAl) sumw <- `names<-`(list(fsumCpp(w,g[[1L]],g[[2L]],na.rm)), paste0("sum.", wsym)) else if(keep.group_vars)
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
          return(setAttributes(c(g[[4L]], sumw, .Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE)), ax))
        } else {
          ax[["names"]] <- c(names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(sumw, .Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE)), ax))
        }
      } else return(setAttributes(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(ax[["names"]][gn2], ax[["names"]][-gn])
      return(setAttributes(c(x[gn2],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,TRUE,FALSE),g[[2L]],TRAtoInt(TRA)))
}



fvar <- function(x, ...) { # g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, drop = TRUE, keep.group_vars = TRUE, keep.w = TRUE,
  UseMethod("fvar", x)
}
fvar.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fvarsd,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fvarsd,x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fvarsd,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fvarsd,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE), group_names.GRP(g))) else
        return(.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fvarsd,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fvarsd,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fvarsd,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fvarsd,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fvar.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fvarsdm,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fvarsdm,x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fvarsdm,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fvarsdm,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fvarsdm,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fvarsdm,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fvarsdm,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fvarsdm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fvar.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fvarsdl,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,drop)) else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fvarsdl,x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE,FALSE), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fvarsdl,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fvarsdl,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE), groups)) else
          return(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fvarsdl,x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fvarsdl,x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fvarsdl,x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fvar.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                            keep.group_vars = TRUE, keep.w = TRUE, stable.algo = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- attr(x, "names")
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  sumw <- NULL

  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- unclass(x)[[wn]]
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    onlyw <- !length(gn)
    gn <- c(gn, wn)
    if(keep.w) {
      if(nTRAl) sumw <- `names<-`(list(fsumCpp(w,g[[1L]],g[[2L]],na.rm)), paste0("sum.", wsym)) else if(keep.group_vars)
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
          return(setAttributes(c(g[[4L]], sumw, .Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE)), ax))
        } else {
          ax[["names"]] <- c(names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(sumw, .Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE)), ax))
        }
      } else return(setAttributes(.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(ax[["names"]][gn2], ax[["names"]][-gn])
      return(setAttributes(c(x[gn2],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fvarsdl,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fvarsdl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE,FALSE),g[[2L]],TRAtoInt(TRA)))
}



