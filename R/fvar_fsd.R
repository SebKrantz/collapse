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

fsd <- function(x, ...) { # g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, drop = TRUE, keep.group_vars = TRUE, keep.w = TRUE,
  UseMethod("fsd", x)
}
fsd.default <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdCpp(x,0L,0L,NULL,w,na.rm,stable.algo)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fvarsdCpp(x,length(lev),g,NULL,w,na.rm,stable.algo), lev))
      } else {
        if(is.factor(g)) return(fvarsdCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo)) else {
          g <- qG(g)
          return(fvarsdCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fvarsdCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo), group_names.GRP(g))) else
        return(fvarsdCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fvarsdCpp(x,0L,0L,NULL,w,na.rm,stable.algo),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fvarsdCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fvarsdCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fvarsdCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsd.matrix <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdmCpp(x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fvarsdmCpp(x,length(lev),g,NULL,w,na.rm,stable.algo), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(fvarsdmCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo)) else {
          g <- qG(g)
          return(fvarsdmCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fvarsdmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(fvarsdmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fvarsdmCpp(x,0L,0L,NULL,w,na.rm,stable.algo),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fvarsdmCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fvarsdmCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fvarsdmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsd.data.frame <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdlCpp(x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fvarsdlCpp(x,length(lev),g,NULL,w,na.rm,stable.algo), lev))
      } else {
        if(is.factor(g)) return(fvarsdlCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo)) else {
          g <- qG(g)
          return(fvarsdlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(fvarsdlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo), groups)) else
          return(fvarsdlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo))
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fvarsdlCpp(x,0L,0L,NULL,w,na.rm,stable.algo),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fvarsdlCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fvarsdlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fvarsdlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsd.grouped_df <- function(x, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE, stable.algo = TRUE,
                           keep.group_vars = TRUE, keep.w = TRUE, ...) {
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- names(x)
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- TRA == FALSE
  sumw <- NULL

  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- x[[wn]]
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    gn <- c(gn, wn)
    if(keep.w) {
     if(nTRAl) sumw <- `names<-`(list(fsumCpp(w,g[[1L]],g[[2L]],na.rm)), paste0("sum.", wsym)) else if(keep.group_vars)
     gn2 <- gn else sumw <- gn2 <- wn
    }
  }

  gl <- length(gn) > 0L # necessary here, not before !!!

  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(g[[4L]], sumw, fvarsdlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo)), ax))
        } else {
          ax[["names"]] <- c(names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(sumw, fvarsdlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo)), ax))
        }
      } else return(setAttributes(fvarsdlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(ax[["names"]][gn2], ax[["names"]][-gn])
      return(setAttributes(c(x[gn2],TRAlCpp(x[-gn],fvarsdlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(TRAlCpp(x[-gn],fvarsdlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(TRAlCpp(x,fvarsdlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo),g[[2L]],TRAtoInt(TRA)))
}


fvar <- function(x, ...) { # g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, drop = TRUE, keep.group_vars = TRUE, keep.w = TRUE,
  UseMethod("fvar", x)
}
fvar.default <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fvarsdCpp(x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE), lev))
      } else {
        if(is.factor(g)) return(fvarsdCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE)) else {
          g <- qG(g)
          return(fvarsdCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fvarsdCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE), group_names.GRP(g))) else
        return(fvarsdCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fvarsdCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fvarsdCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fvarsdCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fvarsdCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fvar.matrix <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdmCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fvarsdmCpp(x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(fvarsdmCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE)) else {
          g <- qG(g)
          return(fvarsdmCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fvarsdmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(fvarsdmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fvarsdmCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fvarsdmCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fvarsdmCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fvarsdmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fvar.data.frame <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdlCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,drop)) else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fvarsdlCpp(x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE), lev))
      } else {
        if(is.factor(g)) return(fvarsdlCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE)) else {
          g <- qG(g)
          return(fvarsdlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(fvarsdlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE), groups)) else
          return(fvarsdlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE))
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fvarsdlCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fvarsdlCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fvarsdlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fvarsdlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fvar.grouped_df <- function(x, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE, stable.algo = TRUE,
                            keep.group_vars = TRUE, keep.w = TRUE, ...) {
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- names(x)
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- TRA == FALSE
  sumw <- NULL

  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- x[[wn]]
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    gn <- c(gn, wn)
    if(keep.w) {
      if(nTRAl) sumw <- `names<-`(list(fsumCpp(w,g[[1L]],g[[2L]],na.rm)), paste0("sum.", wsym)) else if(keep.group_vars)
        gn2 <- gn else sumw <- gn2 <- wn
    }
  }

  gl <- length(gn) > 0L # necessary here, not before !!!

  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(g[[4L]], sumw, fvarsdlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE)), ax))
        } else {
          ax[["names"]] <- c(names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(sumw, fvarsdlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE)), ax))
        }
      } else return(setAttributes(fvarsdlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE), ax))
    } else if(keep.group_vars || (keep.w && length(sumw))) {
      ax[["names"]] <- c(ax[["names"]][gn2], ax[["names"]][-gn])
      return(setAttributes(c(x[gn2],TRAlCpp(x[-gn],fvarsdlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(TRAlCpp(x[-gn],fvarsdlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(TRAlCpp(x,fvarsdlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,stable.algo,FALSE),g[[2L]],TRAtoInt(TRA)))
}

