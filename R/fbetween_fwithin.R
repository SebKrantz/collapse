
# Note: for principal innovations of this code see fsum.R and fscale.R. Old code is commented out below and was innovated in flag.R.
#  replaced give.names = TRUE with stub

ckm <- function(x) if(is.double(x)) x else if(is.character(x) && x == "overall.mean") -Inf else stop("mean must be a number or 'overall.mean'") # better than switch !!

fwithin <- function(x, ...) { # g = NULL, w = NULL, na.rm = TRUE, mean = 0,
  UseMethod("fwithin", x)
}
fwithin.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BW,x,0L,0L,NULL,w,na.rm,ckm(mean),FALSE,FALSE)) else if (is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BW,x,attr(g,"N.groups"),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BW,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ckm(mean),FALSE,FALSE))
  }
}
fwithin.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE)
}
fwithin.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BWm,x,0L,0L,NULL,w,na.rm,ckm(mean),FALSE,FALSE)) else if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BWm,x,fnlevels(g),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWm,x,attr(g,"N.groups"),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BWm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ckm(mean),FALSE,FALSE))
  }
}
fwithin.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,ckm(mean),FALSE,FALSE)) else if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWl,x,attr(g,"N.groups"),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ckm(mean),FALSE,FALSE))
  }
}
fwithin.pdata.frame <- function(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE)
}
fwithin.grouped_df <- function(x, w = NULL, na.rm = TRUE, mean = 0,
                               keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- attr(x, "names")
  gn2 <- which(nam %in% g[[5L]])
  gn <- if(keep.group_vars) gn2 else NULL
  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- unclass(x)[[wn]]
    if(any(gn2 == wn)) stop("Weights coincide with grouping variables!")
    gn2 <- c(gn2,wn)
    if(keep.w) gn <- c(gn,wn)
  }
  if(length(gn2)) {
    if(!length(gn))
      return(.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,ckm(mean),FALSE,FALSE)) else {
        ax <- attributes(x)
        attributes(x) <- NULL
        ax[["names"]] <- c(nam[gn], nam[-gn2])
        return(setAttributes(c(x[gn],.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,ckm(mean),FALSE,FALSE)), ax))
      }
  } else return(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ckm(mean),FALSE,FALSE))
}


W <- function(x, ...) { # g = NULL, w = NULL, na.rm = TRUE, mean = 0,
  UseMethod("W", x)
}
W.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, ...) {
  fwithin.default(x, g, w, na.rm, mean, ...)
}
W.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE)
}
W.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, stub = "W.", ...) {
  add_stub(fwithin.matrix(x, g, w, na.rm, mean, ...), stub)
}
W.grouped_df <- function(x, w = NULL, na.rm = TRUE, mean = 0,
                         stub = "W.", keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- attr(x, "names")
  gn2 <- which(nam %in% g[[5L]])
  gn <- if(keep.group_vars) gn2 else NULL
  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- unclass(x)[[wn]]
    if(any(gn2 == wn)) stop("Weights coincide with grouping variables!")
    gn2 <- c(gn2,wn)
    if(keep.w) gn <- c(gn, wn)
  }
  if(length(gn2)) {
    if(!length(gn))
      return(add_stub(.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,ckm(mean),FALSE,FALSE), stub)) else {
        ax <- attributes(x)
        attributes(x) <- NULL
        ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[-gn2]) else nam[-gn2])
        return(setAttributes(c(x[gn],.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,ckm(mean),FALSE,FALSE)), ax))
      }
  } else return(add_stub(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,ckm(mean),FALSE,FALSE), stub))
}
W.pdata.frame <- function(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE, mean = 0,
                          stub = "W.", keep.ids = TRUE, keep.w = TRUE, ...) {

  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  ax <- attributes(x)
  class(x) <- NULL
  nam <- names(x)
  g <- if(length(effect) == 1L) unclass(ax[["index"]])[[effect]] else
    finteraction(unclass(ax[["index"]])[effect])

  if(keep.ids) {
    gn <- which(nam %in% attr(ax[["index"]], "names"))
    if(length(gn) && is.null(cols)) cols <- seq_along(x)[-gn]
  } else gn <- NULL

  if(!is.null(cols)) cols <- cols2int(cols, x, nam)

  if(is.call(w)) {
    w <- all.vars(w)
    wn <- ckmatch(w, nam)
    w <- x[[wn]]
    cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
    if(keep.w) gn <- c(gn, wn)
  }

  if(length(gn) && !is.null(cols)) {
    ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols])
    return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE)), ax))
  } else if(!length(gn)) {
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
    return(setAttributes(.Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE), ax))
  } else {
    if(is.character(stub)) {
      ax[["names"]] <- paste0(stub, nam)
      return(setAttributes(.Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE), ax))
    } else
      return(.Call(Cpp_BWl,`oldClass<-`(x, ax[["class"]]),fnlevels(g),g,NULL,w,na.rm,ckm(mean),FALSE,FALSE))
  }
}
W.data.frame <- function(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
                         mean = 0, stub = "W.", keep.by = TRUE, keep.w = TRUE, ...) {

  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.call(by) || is.call(w)) {
    ax <- attributes(x)
    class(x) <- NULL
    nam <- names(x)

    if(is.call(by)) {
      if(length(by) == 3L) {
        cols <- ckmatch(all.vars(by[[2L]]), nam)
        gn <- ckmatch(all.vars(by[[3L]]), nam)
      } else {
        gn <- ckmatch(all.vars(by), nam)
        cols <- if(is.null(cols)) seq_along(x)[-gn] else cols2int(cols, x, nam)
      }
      by <- if(length(gn) == 1L) at2GRP(x[[gn]]) else GRP.default(x, gn, return.groups = FALSE)
      if(!keep.by) gn <- NULL
    } else {
      gn <- NULL
      if(!is.null(cols)) cols <- cols2int(cols, x, nam)
      if(!is.GRP(by)) by <- if(is.null(by)) list(0L, 0L, NULL) else if(is.atomic(by)) # Necessary for if by is passed externally !!
        at2GRP(by) else GRP.default(by, return.groups = FALSE)
    }

    if(is.call(w)) {
      w <- all.vars(w)
      wn <- ckmatch(w, nam)
      w <- x[[wn]]
      cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
      if(keep.w) gn <- c(gn, wn)
    }

    if(length(gn)) {
      ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols])
      return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,ckm(mean),FALSE,FALSE)), ax))
    } else {
      ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
      return(setAttributes(.Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,ckm(mean),FALSE,FALSE), ax))
    }
  } else if(!is.null(cols)) { # Need to do like this, otherwise list-subsetting drops attributes !!
    ax <- attributes(x)
    x <- unclass(x)[cols2int(cols, x, ax[["names"]])]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }
  if(is.character(stub)) attr(x, "names") <- paste0(stub, attr(x, "names"))

  if(is.null(by)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,ckm(mean),FALSE,FALSE)) else if (is.atomic(by)) {
    if(is.nmfactor(by)) return(.Call(Cpp_BWl,x,fnlevels(by),by,NULL,w,na.rm,ckm(mean),FALSE,FALSE)) else {
      by <- qG(by, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWl,x,attr(by,"N.groups"),by,NULL,w,na.rm,ckm(mean),FALSE,FALSE))
    }
  } else {
    if(!is.GRP(by)) by <- GRP.default(by, return.groups = FALSE)
    return(.Call(Cpp_BWl,x,by[[1L]],by[[2L]],by[[3L]],w,na.rm,ckm(mean),FALSE,FALSE))
  }
}


fbetween <- function(x, ...) { # g = NULL, w = NULL, na.rm = TRUE, fill = FALSE,
  UseMethod("fbetween", x)
}
fbetween.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BW,x,0L,0L,NULL,w,na.rm,0,TRUE,fill)) else if (is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,0,TRUE,fill)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BW,x,attr(g,"N.groups"),g,NULL,w,na.rm,0,TRUE,fill))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BW,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,0,TRUE,fill))
  }
}
fbetween.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,0,TRUE,fill)
}
fbetween.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BWm,x,0L,0L,NULL,w,na.rm,0,TRUE,fill)) else if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BWm,x,fnlevels(g),g,NULL,w,na.rm,0,TRUE,fill)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWm,x,attr(g,"N.groups"),g,NULL,w,na.rm,0,TRUE,fill))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BWm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,0,TRUE,fill))
  }
}
fbetween.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,0,TRUE,fill)) else if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,0,TRUE,fill)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWl,x,attr(g,"N.groups"),g,NULL,w,na.rm,0,TRUE,fill))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,0,TRUE,fill))
  }
}
fbetween.pdata.frame <- function(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,0,TRUE,fill)
}
fbetween.grouped_df <- function(x, w = NULL, na.rm = TRUE, fill = FALSE,
                                keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- attr(x, "names")
  gn2 <- which(nam %in% g[[5L]])
  gn <- if(keep.group_vars) gn2 else NULL
  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- unclass(x)[[wn]]
    if(any(gn2 == wn)) stop("Weights coincide with grouping variables!")
    gn2 <- c(gn2,wn)
    if(keep.w) gn <- c(gn,wn)
  }
  if(length(gn2)) {
    if(!length(gn))
      return(.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,0,TRUE,fill)) else {
        ax <- attributes(x)
        attributes(x) <- NULL
        ax[["names"]] <- c(nam[gn], nam[-gn2])
        return(setAttributes(c(x[gn],.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,0,TRUE,fill)), ax))
      }
  } else return(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,0,TRUE,fill))
}


B <- function(x, ...) { # g = NULL, w = NULL, na.rm = TRUE, fill = FALSE,
  UseMethod("B", x)
}
B.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  fbetween.default(x, g, w, na.rm, fill, ...)
}
B.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,0,TRUE,fill)
}
B.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, stub = "B.", ...) {
  add_stub(fbetween.matrix(x, g, w, na.rm, fill, ...), stub)
}
B.grouped_df <- function(x, w = NULL, na.rm = TRUE, fill = FALSE,
                         stub = "B.", keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- attr(x, "names")
  gn2 <- which(nam %in% g[[5L]])
  gn <- if(keep.group_vars) gn2 else NULL
  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- unclass(x)[[wn]]
    if(any(gn2 == wn)) stop("Weights coincide with grouping variables!")
    gn2 <- c(gn2,wn)
    if(keep.w) gn <- c(gn, wn)
  }
  if(length(gn2)) {
    if(!length(gn))
      return(add_stub(.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,0,TRUE,fill), stub)) else {
        ax <- attributes(x)
        attributes(x) <- NULL
        ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[-gn2]) else nam[-gn2])
        return(setAttributes(c(x[gn],.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,0,TRUE,fill)), ax))
      }
  } else return(add_stub(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,0,TRUE,fill), stub))
}
B.pdata.frame <- function(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE,
                          stub = "B.", keep.ids = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  ax <- attributes(x)
  class(x) <- NULL
  nam <- names(x)
  g <- if(length(effect) == 1L) unclass(ax[["index"]])[[effect]] else
    finteraction(unclass(ax[["index"]])[effect])

  if(keep.ids) {
    gn <- which(nam %in% attr(ax[["index"]], "names"))
    if(length(gn) && is.null(cols)) cols <- seq_along(x)[-gn]
  } else gn <- NULL

  if(!is.null(cols)) cols <- cols2int(cols, x, nam)

  if(is.call(w)) {
    w <- all.vars(w)
    wn <- ckmatch(w, nam)
    w <- x[[wn]]
    cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
    if(keep.w) gn <- c(gn, wn)
  }

  if(length(gn) && !is.null(cols)) {
    ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols])
    return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,0,TRUE,fill)), ax))
  } else if(!length(gn)) {
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
    return(setAttributes(.Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,0,TRUE,fill), ax))
  } else {
    if(is.character(stub)) {
      ax[["names"]] <- paste0(stub, nam)
      return(setAttributes(.Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,0,TRUE,fill), ax))
    } else
      return(.Call(Cpp_BWl,`oldClass<-`(x, ax[["class"]]),fnlevels(g),g,NULL,w,na.rm,0,TRUE,fill))
  }
}
B.data.frame <- function(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
                         fill = FALSE, stub = "B.", keep.by = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.call(by) || is.call(w)) {
    ax <- attributes(x)
    class(x) <- NULL
    nam <- names(x)

    if(is.call(by)) {
      if(length(by) == 3L) {
        cols <- ckmatch(all.vars(by[[2L]]), nam)
        gn <- ckmatch(all.vars(by[[3L]]), nam)
      } else {
        gn <- ckmatch(all.vars(by), nam)
        cols <- if(is.null(cols)) seq_along(x)[-gn] else cols2int(cols, x, nam)
      }
      by <- if(length(gn) == 1L) at2GRP(x[[gn]]) else GRP.default(x, gn, return.groups = FALSE)
      if(!keep.by) gn <- NULL
    } else {
      gn <- NULL
      if(!is.null(cols)) cols <- cols2int(cols, x, nam)
      if(!is.GRP(by)) by <- if(is.null(by)) list(0L, 0L, NULL) else if(is.atomic(by)) # Necessary if by is passed externally !!
        at2GRP(by) else GRP.default(by, return.groups = FALSE)
    }

    if(is.call(w)) {
      w <- all.vars(w)
      wn <- ckmatch(w, nam)
      w <- x[[wn]]
      cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
      if(keep.w) gn <- c(gn, wn)
    }

    if(length(gn)) {
      ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols])
      return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,0,TRUE,fill)), ax))
    } else {
      ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
      return(setAttributes(.Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,0,TRUE,fill), ax))
    }
  } else if(!is.null(cols)) { # Necessary, else attributes are dropped by list-subsetting !!
    ax <- attributes(x)
    x <- unclass(x)[cols2int(cols, x, ax[["names"]])]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }
  if(is.character(stub)) attr(x, "names") <- paste0(stub, attr(x, "names"))

  if(is.null(by)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,0,TRUE,fill)) else if (is.atomic(by)) {
    if(is.nmfactor(by)) return(.Call(Cpp_BWl,x,fnlevels(by),by,NULL,w,na.rm,0,TRUE,fill)) else {
      by <- qG(by, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWl,x,attr(by,"N.groups"),by,NULL,w,na.rm,0,TRUE,fill))
    }
  } else {
    if(!is.GRP(by)) by <- GRP.default(by, return.groups = FALSE)
    return(.Call(Cpp_BWl,x,by[[1L]],by[[2L]],by[[3L]],w,na.rm,0,TRUE,fill))
  }
}
