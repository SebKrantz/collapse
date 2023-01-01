# library(Rcpp)
# sourceCpp('src/BW.cpp')
# sourceCpp('src/BWa.cpp')
# sourceCpp('src/BWl.cpp')
# source("R/GRP.R")
# source("R/small_helper.R")
# source("R/quick_conversion.R")

# Note: for principal innovations of this code see fsum.R and fscale.R. Old code is commented out below and was innovated in flag.R.
#  replaced give.names = TRUE with stub
fwithin <- function(x, ...) { # g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE,
  UseMethod("fwithin", x)
}
fwithin.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BW,x,0L,0L,NULL,w,na.rm,add.global.mean,FALSE)) else if (is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,add.global.mean,FALSE)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BW,x,attr(g,"N.groups"),g,NULL,w,na.rm,add.global.mean,FALSE))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BW,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,add.global.mean,FALSE))
  }
}
fwithin.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, add.global.mean = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,add.global.mean,FALSE)
}
fwithin.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BWm,x,0L,0L,NULL,w,na.rm,add.global.mean,FALSE)) else if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BWm,x,fnlevels(g),g,NULL,w,na.rm,add.global.mean,FALSE)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWm,x,attr(g,"N.groups"),g,NULL,w,na.rm,add.global.mean,FALSE))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BWm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,add.global.mean,FALSE))
  }
}
fwithin.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,add.global.mean,FALSE)) else if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,add.global.mean,FALSE)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWl,x,attr(g,"N.groups"),g,NULL,w,na.rm,add.global.mean,FALSE))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,add.global.mean,FALSE))
  }
}
fwithin.pdata.frame <- function(x, effect = 1L, w = NULL, na.rm = TRUE, add.global.mean = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,add.global.mean,FALSE)
}
fwithin.grouped_df <- function(x, w = NULL, na.rm = TRUE, add.global.mean = FALSE,
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
      return(.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,add.global.mean,FALSE)) else {
        ax <- attributes(x)
        attributes(x) <- NULL
        ax[["names"]] <- c(nam[gn], nam[-gn2])
        return(setAttributes(c(x[gn],.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,add.global.mean,FALSE)), ax))
      }
  } else return(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,add.global.mean,FALSE))
}


W <- function(x, ...) { # g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE,
  UseMethod("W", x)
}
W.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, ...) {
  fwithin.default(x, g, w, na.rm, add.global.mean, ...)
}
W.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, add.global.mean = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,add.global.mean,FALSE)
}
W.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, stub = "W.", ...) {
  add_stub(fwithin.matrix(x, g, w, na.rm, add.global.mean, ...), stub)
}
W.grouped_df <- function(x, w = NULL, na.rm = TRUE, add.global.mean = FALSE,
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
      return(add_stub(.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,add.global.mean,FALSE), stub)) else {
        ax <- attributes(x)
        attributes(x) <- NULL
        ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[-gn2]) else nam[-gn2])
        return(setAttributes(c(x[gn],.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,add.global.mean,FALSE)), ax))
      }
  } else return(add_stub(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,add.global.mean,FALSE), stub))
}
W.pdata.frame <- function(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE, add.global.mean = FALSE,
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
    return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,add.global.mean,FALSE)), ax))
  } else if(!length(gn)) {
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
    return(setAttributes(.Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,add.global.mean,FALSE), ax))
  } else {
    if(is.character(stub)) {
      ax[["names"]] <- paste0(stub, nam)
      return(setAttributes(.Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,add.global.mean,FALSE), ax))
    } else
      return(.Call(Cpp_BWl,`oldClass<-`(x, ax[["class"]]),fnlevels(g),g,NULL,w,na.rm,add.global.mean,FALSE))
  }
}
W.data.frame <- function(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
                         add.global.mean = FALSE, stub = "W.", keep.by = TRUE, keep.w = TRUE, ...) {
  
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
      return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,add.global.mean,FALSE)), ax))
    } else {
      ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
      return(setAttributes(.Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,add.global.mean,FALSE), ax))
    }
  } else if(!is.null(cols)) { # Need to do like this, otherwise list-subsetting drops attributes !!
    ax <- attributes(x)
    x <- unclass(x)[cols2int(cols, x, ax[["names"]])]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }
  if(is.character(stub)) attr(x, "names") <- paste0(stub, attr(x, "names"))
  
  if(is.null(by)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,add.global.mean,FALSE)) else if (is.atomic(by)) {
    if(is.nmfactor(by)) return(.Call(Cpp_BWl,x,fnlevels(by),by,NULL,w,na.rm,add.global.mean,FALSE)) else {
      by <- qG(by, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWl,x,attr(by,"N.groups"),by,NULL,w,na.rm,add.global.mean,FALSE))
    }
  } else {
    if(!is.GRP(by)) by <- GRP.default(by, return.groups = FALSE)
    return(.Call(Cpp_BWl,x,by[[1L]],by[[2L]],by[[3L]],w,na.rm,add.global.mean,FALSE))
  }
}


fbetween <- function(x, ...) { # g = NULL, w = NULL, na.rm = TRUE, fill = FALSE,
  UseMethod("fbetween", x)
}
fbetween.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BW,x,0L,0L,NULL,w,na.rm,fill,TRUE)) else if (is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,fill,TRUE)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BW,x,attr(g,"N.groups"),g,NULL,w,na.rm,fill,TRUE))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BW,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,fill,TRUE))
  }
}
fbetween.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,fill,TRUE)
}
fbetween.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BWm,x,0L,0L,NULL,w,na.rm,fill,TRUE)) else if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BWm,x,fnlevels(g),g,NULL,w,na.rm,fill,TRUE)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWm,x,attr(g,"N.groups"),g,NULL,w,na.rm,fill,TRUE))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BWm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,fill,TRUE))
  }
}
fbetween.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,fill,TRUE)) else if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,fill,TRUE)) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWl,x,attr(g,"N.groups"),g,NULL,w,na.rm,fill,TRUE))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,fill,TRUE))
  }
}
fbetween.pdata.frame <- function(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,fill,TRUE)
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
      return(.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,fill,TRUE)) else {
        ax <- attributes(x)
        attributes(x) <- NULL
        ax[["names"]] <- c(nam[gn], nam[-gn2])
        return(setAttributes(c(x[gn],.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,fill,TRUE)), ax))
      }
  } else return(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,fill,TRUE))
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
  .Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,fill,TRUE)
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
      return(add_stub(.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,fill,TRUE), stub)) else {
        ax <- attributes(x)
        attributes(x) <- NULL
        ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[-gn2]) else nam[-gn2])
        return(setAttributes(c(x[gn],.Call(Cpp_BWl,x[-gn2],g[[1L]],g[[2L]],g[[3L]],w,na.rm,fill,TRUE)), ax))
      }
  } else return(add_stub(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,fill,TRUE), stub))
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
    return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,fill,TRUE)), ax))
  } else if(!length(gn)) {
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
    return(setAttributes(.Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,fill,TRUE), ax))
  } else {
    if(is.character(stub)) {
      ax[["names"]] <- paste0(stub, nam)
      return(setAttributes(.Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,fill,TRUE), ax))
    } else
      return(.Call(Cpp_BWl,`oldClass<-`(x, ax[["class"]]),fnlevels(g),g,NULL,w,na.rm,fill,TRUE))
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
      return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,fill,TRUE)), ax))
    } else {
      ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
      return(setAttributes(.Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,fill,TRUE), ax))
    }
  } else if(!is.null(cols)) { # Necessary, else attributes are dropped by list-subsetting !!
    ax <- attributes(x)
    x <- unclass(x)[cols2int(cols, x, ax[["names"]])]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }
  if(is.character(stub)) attr(x, "names") <- paste0(stub, attr(x, "names"))
  
  if(is.null(by)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,fill,TRUE)) else if (is.atomic(by)) {
    if(is.nmfactor(by)) return(.Call(Cpp_BWl,x,fnlevels(by),by,NULL,w,na.rm,fill,TRUE)) else {
      by <- qG(by, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_BWl,x,attr(by,"N.groups"),by,NULL,w,na.rm,fill,TRUE))
    }
  } else {
    if(!is.GRP(by)) by <- GRP.default(by, return.groups = FALSE)
    return(.Call(Cpp_BWl,x,by[[1L]],by[[2L]],by[[3L]],w,na.rm,fill,TRUE))
  }
}
