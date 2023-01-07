# library(Rcpp)
# sourceCpp('src/fgrowth.cpp', rebuild = TRUE) # Todo: Thoroughly test these functions !!!
# sourceCpp('src/fgrowtha.cpp', rebuild = TRUE)
# sourceCpp('src/fgrowthl.cpp', rebuild = TRUE)
# source("R/GRP.R")
# source("R/small_helper.R")
# source("R/quick_conversion.R")

# For principle innovations of this code see flag.R and flag.cpp # stubs instead of give.names !!

# ,logdiff,stubs  # TRY   x, n = 1, diff = 1, ..., t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE
fgrowth <- function(x, n = 1, diff = 1, ...) UseMethod("fgrowth") # , x

fgrowth.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.default", dotstostr(...)))
  if(is.null(g))
    return(.Call(Cpp_fgrowth,x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),logdiff,stubs)) else if(is.atomic(g)) {
      if(is.nmfactor(g)) nl <- fnlevels(g) else {
        g <- qG(g, na.exclude = FALSE)
        nl <- attr(g, "N.groups")
      }
      return(.Call(Cpp_fgrowth,x,n,diff,fill,nl,g,NULL,G_t(t,wm=3L),logdiff,stubs))
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      return(.Call(Cpp_fgrowth,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=3L),logdiff,stubs))
    }
}
fgrowth.pseries <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.pseries", dotstostr(...)))
  index <- unclass(attr(x, "index"))
  if(length(index) > 2L) index <- c(finteraction(index[-length(index)]), index[length(index)])
  if(is.matrix(x))
  .Call(Cpp_fgrowthm,x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs) else
  .Call(Cpp_fgrowth,x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs)
}
fgrowth.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.matrix", dotstostr(...)))
  if(is.null(g))
    return(.Call(Cpp_fgrowthm,x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),logdiff,stubs)) else if(is.atomic(g)) {
      if(is.nmfactor(g)) nl <- fnlevels(g) else {
        g <- qG(g, na.exclude = FALSE)
        nl <- attr(g, "N.groups")
      }
      .Call(Cpp_fgrowthm,x,n,diff,fill,nl,g,NULL,G_t(t,wm=3L),logdiff,stubs)
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      .Call(Cpp_fgrowthm,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=3L),logdiff,stubs)
    }
}
fgrowth.grouped_df <- function(x, n = 1, diff = 1, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, keep.ids = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.grouped_df", dotstostr(...)))
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
      return(.Call(Cpp_fgrowthl,x[-gn],n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=3L),logdiff,stubs)) else {
        ax <- attributes(x)
        class(x) <- NULL # Works for multiple lags !!
        res <- c(x[gn],.Call(Cpp_fgrowthl,x[-gn],n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=3L),logdiff,stubs))
        ax[["names"]] <- names(res)
        return(setAttributes(res, ax))
      }
  } else return(.Call(Cpp_fgrowthl,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=3L),logdiff,stubs))
}
fgrowth.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.data.frame", dotstostr(...)))
  if(is.null(g))
    return(.Call(Cpp_fgrowthl,x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),logdiff,stubs)) else if(is.atomic(g)) {
      if(is.nmfactor(g)) nl <- fnlevels(g) else {
        g <- qG(g, na.exclude = FALSE)
        nl <- attr(g, "N.groups")
      }
      .Call(Cpp_fgrowthl,x,n,diff,fill,nl,g,NULL,G_t(t,wm=3L),logdiff,stubs)
    } else {
      if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE)
      .Call(Cpp_fgrowthl,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t,wm=3L),logdiff,stubs)
    }
}
fgrowth.pdata.frame <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.pdata.frame", dotstostr(...)))
  index <- unclass(attr(x, "index"))
  if(length(index) > 2L) index <- c(finteraction(index[-length(index)]), index[length(index)])
  .Call(Cpp_fgrowthl,x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs)
}

# Growth Operator
G <- function(x, n = 1, diff = 1, ...) UseMethod("G") # , x

G.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...)
  fgrowth.default(x, n, diff, g, t, fill, logdiff, stubs, ...)

G.pseries <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, ...)
  fgrowth.pseries(x, n, diff, fill, logdiff, stubs, ...)

G.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...)
  fgrowth.matrix(x, n, diff, g, t, fill, logdiff, stubs, ...)

G.grouped_df <- fgrowth.grouped_df

G.data.frame <- function(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, logdiff = FALSE, stubs = TRUE, keep.ids = TRUE, ...) {

  if(!missing(...)) stop(sprintf("Unknown argument %s passed to G.data.frame", dotstostr(...)))
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
      c(x[gn], .Call(Cpp_fgrowthl,x[cols],n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t,wm=3L),logdiff,stubs)) else
        .Call(Cpp_fgrowthl,x[cols],n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t,wm=3L),logdiff,stubs)
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!is.null(cols)) { # Needs to be like this, otherwise list-subsetting removes attributes !!
    ax <- attributes(x)
    class(x) <- NULL
    x <- x[cols2int(cols, x, names(x))]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }

  if(is.null(by))
    return(.Call(Cpp_fgrowthl,x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),logdiff,stubs)) else if(is.atomic(by)) {
      if(is.nmfactor(by)) nl <- fnlevels(by) else {
        by <- qG(by, na.exclude = FALSE)
        nl <- attr(by, "N.groups")
      }
      .Call(Cpp_fgrowthl,x,n,diff,fill,nl,by,NULL,G_t(t,wm=3L),logdiff,stubs)
    } else {
      if(!is.GRP(by)) by <- GRP.default(by, return.groups = FALSE)
      .Call(Cpp_fgrowthl,x,n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t,wm=3L),logdiff,stubs)
    }
}

G.pdata.frame <- function(x, n = 1, diff = 1, cols = is.numeric, fill = NA, logdiff = FALSE, stubs = TRUE, keep.ids = TRUE, ...) {

  if(!missing(...)) stop(sprintf("Unknown argument %s passed to G.pdata.frame", dotstostr(...)))
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
    res <- c(x[gn], .Call(Cpp_fgrowthl,x[cols],n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs))
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!length(gn)) # could speed up ??
    return(.Call(Cpp_fgrowthl,fcolsubset(x, cols),n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs)) else
      return(.Call(Cpp_fgrowthl,x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs))
}
