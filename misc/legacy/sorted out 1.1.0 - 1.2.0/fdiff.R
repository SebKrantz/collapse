# library(Rcpp)
# sourceCpp('src/fdiff.cpp', rebuild = TRUE)
# sourceCpp('src/fdiffa.cpp', rebuild = TRUE)
# sourceCpp('src/fdiffl.cpp', rebuild = TRUE) # On large test data (10 mio obs), unordered panel-difference still gave error !!!!!
# source("R/GRP.R")
# source("R/small_helper.R")
# source("R/quick_conversion.R")

# For principle innovations of this code see flag.R and flag.cpp # stubs instead of give.names !!

fdiff <- function(x, n = 1, diff = 1, ...) UseMethod("fdiff") # , x

fdiff.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fdiff.default", dotstostr(...)))
  if(is.null(g))
    return(.Call(Cpp_fdiff,x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(g)) {
      if(is.nmfactor(g)) nl <- fnlevels(g) else {
        g <- qG(g, na.exclude = FALSE)
        nl <- attr(g, "N.groups")
      }
      return(.Call(Cpp_fdiff,x,n,diff,fill,nl,g,NULL,G_t(t,wm=2L),stubs))
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_fdiff,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=2L),stubs))
    }
}
fdiff.pseries <- function(x, n = 1, diff = 1, fill = NA, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fdiff.pseries", dotstostr(...)))
  index <- unclass(attr(x, "index"))
  if(length(index) > 2L) index <- c(finteraction(index[-length(index)]), index[length(index)])
  if(is.matrix(x))
    .Call(Cpp_fdiffm,x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs) else
      .Call(Cpp_fdiff,x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)
}
fdiff.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fdiff.matrix", dotstostr(...)))
  if(is.null(g))
    return(.Call(Cpp_fdiffm,x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(g)) {
      if(is.nmfactor(g)) nl <- fnlevels(g) else {
        g <- qG(g, na.exclude = FALSE)
        nl <- attr(g, "N.groups")
      }
      .Call(Cpp_fdiffm,x,n,diff,fill,nl,g,NULL,G_t(t,wm=2L),stubs)
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      .Call(Cpp_fdiffm,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=2L),stubs)
    }
}
fdiff.grouped_df <- function(x, n = 1, diff = 1, t = NULL, fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fdiff.grouped_df", dotstostr(...)))
  g <- GRP.grouped_df(x)
  tsym <- l1orn(all.vars(substitute(t)), "NULL")
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  if(!(tsym == "NULL" || is.na(tn <- match(tsym, nam)))) {
    if(any(gn == tn)) stop("timevar coincides with grouping variables!")
    t <- unclass(x)[[tn]]
    gn <- c(gn, tn)
  }
  if(length(gn)) {
    if(!keep.ids)
      return(.Call(Cpp_fdiffl,x[-gn],n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=2L),stubs)) else {
        ax <- attributes(x)
        class(x) <- NULL # Works for multiple lags !!
        res <- c(x[gn],.Call(Cpp_fdiffl,x[-gn],n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=2L),stubs))
        ax[["names"]] <- names(res)
        return(setAttributes(res, ax))
      }
  } else return(.Call(Cpp_fdiffl,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=2L),stubs))
}
fdiff.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fdiff.data.frame", dotstostr(...)))
  if(is.null(g))
    return(.Call(Cpp_fdiffl,x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(g)) {
      if(is.nmfactor(g)) nl <- fnlevels(g) else {
        g <- qG(g, na.exclude = FALSE)
        nl <- attr(g, "N.groups")
      }
      .Call(Cpp_fdiffl,x,n,diff,fill,nl,g,NULL,G_t(t,wm=2L),stubs)
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      .Call(Cpp_fdiffl,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=2L),stubs)
    }
}
fdiff.pdata.frame <- function(x, n = 1, diff = 1, fill = NA, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fdiff.pdata.frame", dotstostr(...)))
  index <- unclass(attr(x, "index"))
  if(length(index) > 2L) index <- c(finteraction(index[-length(index)]), index[length(index)])
  .Call(Cpp_fdiffl,x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)
}

# Difference Operator (masks stats::D)
# use xt instead of by ???
D <- function(x, n = 1, diff = 1, ...) UseMethod("D") # , x

D.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...)
  fdiff.default(x, n, diff, g, t, fill, stubs, ...)

D.pseries <- function(x, n = 1, diff = 1, fill = NA, stubs = TRUE, ...)
  fdiff.pseries(x, n, diff, fill, stubs, ...)

D.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...)
  fdiff.matrix(x, n, diff, g, t, fill, stubs, ...)

D.grouped_df <- fdiff.grouped_df

D.data.frame <- function(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {

  if(!missing(...)) stop(sprintf("Unknown argument %s passed to D.data.frame", dotstostr(...)))
  if(is.call(by) || is.call(t)) {
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
      if(!keep.ids) gn <- NULL
    } else {
      gn <- NULL
      if(!is.null(cols)) cols <- cols2int(cols, x, nam)
      if(!is.GRP(by)) by <- if(is.null(by)) list(0L, 0L, NULL) else if(is.atomic(by)) # Necessary for if by is passed externally !!
        at2GRP(by) else GRP.default(by, return.groups = FALSE)
    }

    if(is.call(t)) {
      t <- all.vars(t)
      tn <- ckmatch(t, nam)
      t1 <- length(tn) == 1L
      t <- if(t1) x[[tn]] else GRP.default(x[tn], return.groups = FALSE)[[2L]]
      cols <- if(is.null(cols)) seq_along(x)[-tn] else if(t1) cols[cols != tn] else fsetdiff(cols, tn)
      if(keep.ids) gn <- c(gn, tn)
    }

    res <- if(length(gn))
      c(x[gn], .Call(Cpp_fdiffl,x[cols],n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t,wm=2L),stubs)) else
        .Call(Cpp_fdiffl,x[cols],n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t,wm=2L),stubs)
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!is.null(cols)) { # Needs to be done like this, otherwise list-subsetting drops attributes !!
    ax <- attributes(x)
    class(x) <- NULL
    x <- x[cols2int(cols, x, names(x))]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }

  if(is.null(by))
    return(.Call(Cpp_fdiffl,x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(by)) {
      if(is.nmfactor(by)) nl <- fnlevels(by) else {
        by <- qG(by, na.exclude = FALSE)
        nl <- attr(by, "N.groups")
      }
      .Call(Cpp_fdiffl,x,n,diff,fill,nl,by,NULL,G_t(t,wm=2L),stubs)
    } else {
      if(!is.GRP(by)) by <- GRP.default(by, return.groups = FALSE)
      .Call(Cpp_fdiffl,x,n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t,wm=2L),stubs)
    }
}

D.pdata.frame <- function(x, n = 1, diff = 1, cols = is.numeric, fill = NA, stubs = TRUE,
                          keep.ids = TRUE, ...) {

  if(!missing(...)) stop(sprintf("Unknown argument %s passed to D.pdata.frame", dotstostr(...)))
  ax <- attributes(x)
  nam <- ax[["names"]]
  index <- unclass(ax[["index"]])

  if(keep.ids) {
    gn <- which(nam %in% names(index))
    if(length(gn) && is.null(cols)) cols <- seq_along(unclass(x))[-gn]
  } else gn <- NULL

  if(length(index) > 2L) index <- c(finteraction(index[-length(index)]), index[length(index)])

  if(!is.null(cols)) cols <- cols2int(cols, x, nam)

  if(length(gn) && !is.null(cols)) {
    class(x) <- NULL # Works for multiple lags !!
    res <- c(x[gn], .Call(Cpp_fdiffl,x[cols],n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs))
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!length(gn)) # could speed up ??
    return(.Call(Cpp_fdiffl,fcolsubset(x, cols),n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)) else
      return(.Call(Cpp_fdiffl,x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs))
}
