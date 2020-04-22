# library(Rcpp)
# sourceCpp('src/fmax.cpp')
# sourceCpp('src/fmode.cpp') # Still test thoroughly !!, make better fmodelCpp !!!
# sourceCpp('src/fmodea.cpp') # Still test thoroughly !!
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# Note: for principal innovations of this code see fsum.R !!

fmode <- function(x, ...) UseMethod("fmode") # , x

fmode.default <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fmode,x,0L,0L,NULL,w,na.rm)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fmode,x,length(lev),g,NULL,w,na.rm), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fmode,x,fnlevels(g),g,NULL,w,na.rm)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fmode,x,attr(g,"N.groups"),g,NULL,w,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fmode,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_fmode,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmode,x,0L,0L,NULL,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmode,x,fnlevels(g),g,NULL,w,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fmode,x,attr(g,"N.groups"),g,NULL,w,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fmode,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmode.matrix <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fmodem,x,0L,0L,NULL,w,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fmodem,x,length(lev),g,NULL,w,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fmodem,x,fnlevels(g),g,NULL,w,na.rm,FALSE)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fmodem,x,attr(g,"N.groups"),g,NULL,w,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fmodem,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fmodem,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fmodem,x,0L,0L,NULL,w,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fmodem,x,fnlevels(g),g,NULL,w,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fmodem,x,attr(g,"N.groups"),g,NULL,w,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fmodem,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmode.data.frame <- function(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) {
      if(drop) return(unlist(.Call(Cpp_fmodel,x,0L,0L,NULL,w,na.rm))) else return(.Call(Cpp_fmodel,x,0L,0L,NULL,w,na.rm))
    } else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fmodel,x,length(lev),g,NULL,w,na.rm), lev))
      } else {
        if(is.nmfactor(g)) return(.Call(Cpp_fmodel,x,fnlevels(g),g,NULL,w,na.rm)) else {
          g <- qG(g, na.exclude = FALSE)
          return(.Call(Cpp_fmodel,x,attr(g,"N.groups"),g,NULL,w,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP.default(g) else GRP.default(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), groups)) else
          return(.Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fmodel,x,0L,0L,NULL,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fmodel,x,fnlevels(g),g,NULL,w,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g, na.exclude = FALSE)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fmodel,x,attr(g,"N.groups"),g,NULL,w,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmode.grouped_df <- function(x, w = NULL, TRA = NULL, na.rm = TRUE, use.g.names = FALSE,
                             keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  g <- GRP.grouped_df(x)
  wsym <- l1orn(all.vars(substitute(w)), "NULL")
  nam <- attr(x, "names")
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  maxw <- NULL

  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- unclass(x)[[wn]]
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    gn <- c(gn, wn)
    if(keep.w) {
      if(nTRAl) maxw <- `names<-`(list(fmaxCpp(w,g[[1L]],g[[2L]],na.rm)), paste0("max.", wsym)) else if(keep.group_vars)
        gn2 <- gn else maxw <- gn2 <- wn
    }
  }

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
          ax[["names"]] <- c(g[[5L]], names(maxw), nam[-gn])
          return(setAttributes(c(g[[4L]], maxw, .Call(Cpp_fmodel,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm)), ax))
        } else {
          ax[["names"]] <- c(names(maxw), nam[-gn])
          return(setAttributes(c(maxw, .Call(Cpp_fmodel,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm)), ax))
        }
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]], .Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm)), ax))
      } else return(setAttributes(.Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), ax))
    } else if(keep.group_vars || (keep.w && length(maxw))) {
      ax[["names"]] <- c(nam[gn2], nam[-gn])
      return(setAttributes(c(x[gn2],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fmodel,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- nam[-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fmodel,x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fmodel,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
}
