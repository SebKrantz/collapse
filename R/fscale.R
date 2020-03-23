# Make faster ???
cm <- function(x) if(is.double(x)) x else if(is.character(x) && x == "overall.mean") -Inf else if(isFALSE(x)) Inf else stop("mean must be a number, 'overall.mean' or FALSE")
csd <- function(x) if(is.double(x)) x else if(is.character(x) && x == "within.sd") -Inf else stop("sd must be a number or 'within.sd'")

# w.type = "frequency"
# Todo: center and scale arguments -> link to fmean and fsd with TRA ??

fscale <- function(x, ...) { # g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1,
  UseMethod("fscale", x)
}
fscale.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_fscale,x,0L,0L,w,na.rm,cm(mean),csd(sd))) else if (is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_fscale,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd))) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_fscale,x,attr(g,"N.groups"),g,w,na.rm,cm(mean),csd(sd)))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_fscale,x,g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd)))
  }
}
fscale.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, sd = 1, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_fscale,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd))
}
fscale.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_fscalem,x,0L,0L,w,na.rm,cm(mean),csd(sd))) else if (is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_fscalem,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd))) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_fscalem,x,attr(g,"N.groups"),g,w,na.rm,cm(mean),csd(sd)))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_fscalem,x,g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd)))
  }
}
fscale.grouped_df <- function(x, w = NULL, na.rm = TRUE, mean = 0, sd = 1, keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w)) # faster way ???
  nam <- attr(x, "names")
  gn2 <- which(nam %in% g[[5L]])
  gn <- if(keep.group_vars) gn2 else NULL
  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- unclass(x)[[wn]]
    if(any(gn2 == wn)) stop("Weights coincide with grouping variables!")
    gn2 <- c(gn2, wn)
    if(keep.w) gn <- c(gn, wn)
  }
  if(length(gn2)) {
    if(!length(gn))
      return(.Call(Cpp_fscalel,x[-gn2],g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd))) else {
        ax <- attributes(x)
        attributes(x) <- NULL
        ax[["names"]] <- c(nam[gn], nam[-gn2])
        return(setAttributes(c(x[gn],.Call(Cpp_fscalel,x[-gn2],g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd))), ax))
      }
  } else return(.Call(Cpp_fscalel,x,g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd)))
}
fscale.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) return(.Call(Cpp_fscalel,x,0L,0L,w,na.rm,cm(mean),csd(sd))) else if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_fscalel,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd))) else {
      g <- qG(g, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_fscalel,x,attr(g,"N.groups"),g,w,na.rm,cm(mean),csd(sd)))
    }
  } else {
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
    return(.Call(Cpp_fscalel,x,g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd)))
  }
}
fscale.pdata.frame <- function(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, sd = 1, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_fscale,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd))
}


# Todo: Could still implement return = c("cols","all","add")

STD <- function(x, ...) { # g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1,
  UseMethod("STD", x)
}
STD.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1, ...) {
  fscale.default(x, g, w, na.rm, mean, sd)
}
STD.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, sd = 1, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- if(length(effect) == 1L) unclass(attr(x, "index"))[[effect]] else finteraction(unclass(attr(x, "index"))[effect])
  .Call(Cpp_fscale,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd))
}
STD.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1, stub = "STD.", ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  add_stub(fscale.matrix(x, g, w, na.rm, mean, sd), stub)
}
STD.grouped_df <- function(x, w = NULL, na.rm = TRUE, mean = 0, sd = 1, stub = "STD.", keep.group_vars = TRUE, keep.w = TRUE, ...) {
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
      return(add_stub(.Call(Cpp_fscalel,x[-gn2],g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd)), stub)) else {
        ax <- attributes(x)
        attributes(x) <- NULL
        ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[-gn2]) else nam[-gn2])
        return(setAttributes(c(x[gn],.Call(Cpp_fscalel,x[-gn2],g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd))), ax))
      }
  } else return(add_stub(.Call(Cpp_fscalel,x,g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd)), stub))
}

# updated (best) version !!
STD.pdata.frame <- function(x, effect = 1L, w = NULL, cols = is.numeric,
                            na.rm = TRUE, mean = 0, sd = 1, stub = "STD.", keep.ids = TRUE,
                            keep.w = TRUE, ...) {
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
    wn <- ckmatch(w, nam, "Unknown weight variable:")
    w <- x[[wn]]
    cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
    if(keep.w) gn <- c(gn, wn)
  }

  if(length(gn) && !is.null(cols)) {
    ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols])
    return(setAttributes(c(x[gn], .Call(Cpp_fscalel,x[cols],fnlevels(g),g,w,na.rm,cm(mean),csd(sd))), ax))
  } else if(!length(gn)) {
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
    return(setAttributes(.Call(Cpp_fscalel,x[cols],fnlevels(g),g,w,na.rm,cm(mean),csd(sd)), ax))
  } else {
    if(is.character(stub)) {
      ax[["names"]] <- paste0(stub, nam)
      return(setAttributes(.Call(Cpp_fscalel,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd)), ax))
    } else
      return(.Call(Cpp_fscalel,`oldClass<-`(x, ax[["class"]]),fnlevels(g),g,w,na.rm,cm(mean),csd(sd)))
  }
}
# updated, fast and data.table proof version !!!
STD.data.frame <- function(x, by = NULL, w = NULL, cols = is.numeric,
                           na.rm = TRUE, mean = 0, sd = 1, stub = "STD.", keep.by = TRUE,
                           keep.w = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  # fastest solution??
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
      wn <- ckmatch(w, nam, "Unknown weight variable:")
      w <- x[[wn]]
      cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
      if(keep.w) gn <- c(gn, wn)
    }

    if(length(gn)) {
      ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols])
      return(setAttributes(c(x[gn], .Call(Cpp_fscalel,x[cols],by[[1L]],by[[2L]],w,na.rm,cm(mean),csd(sd))), ax))
    } else {
      ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
      return(setAttributes(.Call(Cpp_fscalel,x[cols],by[[1L]],by[[2L]],w,na.rm,cm(mean),csd(sd)), ax))
    }
  } else if(!is.null(cols)) { # Needs to be like this, otherwise subsetting dropps the attributes !!
    ax <- attributes(x)
    x <- unclass(x)[cols2int(cols, x, ax[["names"]])]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }
  if(is.character(stub)) attr(x, "names") <- paste0(stub, attr(x, "names"))

  if(is.null(by)) return(.Call(Cpp_fscalel,x,0L,0L,w,na.rm,cm(mean),csd(sd))) else if (is.atomic(by)) {
    if(is.nmfactor(by)) return(.Call(Cpp_fscalel,x,fnlevels(by),by,w,na.rm,cm(mean),csd(sd))) else {
      by <- qG(by, ordered = FALSE, na.exclude = FALSE)
      return(.Call(Cpp_fscalel,x,attr(by,"N.groups"),by,w,na.rm,cm(mean),csd(sd)))
    }
  } else {
    if(!is.GRP(by)) by <- GRP.default(by, return.groups = FALSE)
    return(.Call(Cpp_fscalel,x,by[[1L]],by[[2L]],w,na.rm,cm(mean),csd(sd)))
  }
}
