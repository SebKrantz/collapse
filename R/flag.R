
flag <- function(x, n = 1, ...) UseMethod("flag") # , x

flag.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(UseMethod("flag", unclass(x)))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_flaglead,x,n,fill,0L,0L,G_t(t,0L),stubs))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) nl <- fnlevels(g) else {
      g <- qG(g, na.exclude = FALSE)
      nl <- attr(g, "N.groups")
    }
    return(.Call(Cpp_flaglead,x,n,fill,nl,g,G_t(t),stubs))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_flaglead,x,n,fill,g[[1L]],g[[2L]],G_t(t),stubs)
}

flag.pseries <- function(x, n = 1, fill = NA, stubs = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- unclass(attr(x, "index"))
  if(length(index) > 2L) index <- list(finteraction(index[-length(index)]), index[[length(index)]])
  if(is.matrix(x))
  .Call(Cpp_flagleadm,x,n,fill,fnlevels(index[[1L]]),index[[1L]],index[[2L]],stubs) else
  .Call(Cpp_flaglead,x,n,fill,fnlevels(index[[1L]]),index[[1L]],index[[2L]],stubs)
}

flag.matrix <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = length(n) > 1L, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_flagleadm,x,n,fill,0L,0L,G_t(t,0L),stubs))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) nl <- fnlevels(g) else {
      g <- qG(g, na.exclude = FALSE)
      nl <- attr(g, "N.groups")
    }
    return(.Call(Cpp_flagleadm,x,n,fill,nl,g,G_t(t),stubs))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_flagleadm,x,n,fill,g[[1L]],g[[2L]],G_t(t),stubs)
}

flag.grouped_df <- function(x, n = 1, t = NULL, fill = NA, stubs = length(n) > 1L, keep.ids = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  tsym <- all.vars(substitute(t))
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  if(length(tsym) && !anyNA(tn <- match(tsym, nam))) {
    if(length(tn) == 1L) {
      if(any(gn == tn)) stop("timevar coincides with grouping variables!")
      t <- .subset2(x, tn)
    } else {
      if(any(gn %in% tn)) stop("timevar coincides with grouping variables!")
      t <- .subset(x, tn)
    }
    gn <- c(gn, tn)
  }
  if(length(gn)) {
    ax <- attributes(x)
    res <- .Call(Cpp_flagleadl, .subset(x, -gn), n,fill,g[[1L]],g[[2L]],G_t(t),stubs)
    if(keep.ids) res <- c(.subset(x, gn), res)
    ax[["names"]] <- names(res)  # Works for multiple lags !
    return(setAttributes(res, ax))
  }
  .Call(Cpp_flagleadl,x,n,fill,g[[1L]],g[[2L]],G_t(t),stubs)
}
flag.data.frame <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = length(n) > 1L, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_flagleadl,x,n,fill,0L,0L,G_t(t,0L),stubs))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) nl <- fnlevels(g) else {
      g <- qG(g, na.exclude = FALSE)
      nl <- attr(g, "N.groups")
    }
    return(.Call(Cpp_flagleadl,x,n,fill,nl,g,G_t(t),stubs))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_flagleadl,x,n,fill,g[[1L]],g[[2L]],G_t(t),stubs)
}
flag.list <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = length(n) > 1L, ...)
  flag.data.frame(x, n, g, t, fill, stubs, ...)

flag.pdata.frame <- function(x, n = 1, fill = NA, stubs = length(n) > 1L, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- unclass(attr(x, "index"))
  if(length(index) > 2L) index <- list(finteraction(index[-length(index)]), index[[length(index)]])
  .Call(Cpp_flagleadl,x,n,fill,fnlevels(index[[1L]]),index[[1L]],index[[2L]],stubs)
}

# Lag Operator   # use xt instead of by ?
L <- function(x, n = 1, ...) UseMethod("L") # , x

L.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(flag.matrix(x, n, g, t, fill, stubs, ...))
  flag.default(x, n, g, t, fill, stubs, ...)
}

L.pseries <- function(x, n = 1, fill = NA, stubs = TRUE, ...)
  flag.pseries(x, n, fill, stubs, ...)

L.matrix <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...)
  flag.matrix(x, n, g, t, fill, stubs, ...)

L.grouped_df <- function(x, n = 1, t = NULL, fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {
  x <- x
  eval(substitute(flag.grouped_df(x, n, t, fill, stubs, keep.ids, ...)))
}

L.data.frame <- function(x, n = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.call(by) || is.call(t)) {
    ax <- attributes(x)
    class(x) <- NULL
    nam <- names(x)

    if(is.call(by)) {
      if(length(by) == 3L) {
        cols <- ckmatch(all.vars(by[[2L]]), nam, "Unknown variables:")
        gn <- ckmatch(all.vars(by[[3L]]), nam, "Unknown variables:")
      } else {
        gn <- ckmatch(all.vars(by), nam, "Unknown variables:")
        cols <- if(is.null(cols)) seq_along(x)[-gn] else cols2int(cols, x, nam)
      }
      by <- if(length(gn) == 1L) at2GRP(x[[gn]]) else GRP.default(x, gn, return.groups = FALSE, call = FALSE)
      if(!keep.ids) gn <- NULL
    } else {
      gn <- NULL
      if(length(cols)) cols <- cols2int(cols, x, nam)
      if(!is.GRP(by)) by <- if(is.null(by)) list(0L, 0L, NULL) else if(is.atomic(by)) # Necessary for if by is passed externally !
        at2GRP(by) else GRP.default(by, return.groups = FALSE, call = FALSE)
    }

    if(is.call(t)) {
      tn <- ckmatch(all.vars(t), nam, "Unknown variables:")
      t1 <- length(tn) == 1L
      t <- if(t1) x[[tn]] else GRP.default(x[tn], return.groups = FALSE, call = FALSE)[[2L]]
      cols <- if(is.null(cols)) seq_along(x)[-tn] else if(t1) cols[cols != tn] else fsetdiff(cols, tn)
      if(keep.ids) gn <- c(gn, tn)
    }

    res <- if(length(gn))
    c(x[gn], .Call(Cpp_flagleadl,x[cols],n,fill,by[[1L]],by[[2L]],G_t(t),stubs)) else
    .Call(Cpp_flagleadl,x[cols],n,fill,by[[1L]],by[[2L]],G_t(t),stubs)
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(length(cols)) { # Needs to be like this, otherwise subsetting dropps the attributes !
    ax <- attributes(x)
    class(x) <- NULL
    x <- x[cols2int(cols, x, names(x), FALSE)]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }

  if(is.null(by)) return(.Call(Cpp_flagleadl,x,n,fill,0L,0L,G_t(t,0L),stubs))
  if(is.atomic(by)) {
    if(is.nmfactor(by)) nl <- fnlevels(by) else {
      by <- qG(by, na.exclude = FALSE)
      nl <- attr(by, "N.groups")
    }
    return(.Call(Cpp_flagleadl,x,n,fill,nl,by,G_t(t),stubs))
  }
  if(!is.GRP(by)) by <- GRP.default(by, return.groups = FALSE, call = FALSE)
  .Call(Cpp_flagleadl,x,n,fill,by[[1L]],by[[2L]],G_t(t),stubs)
}

L.list <- function(x, n = 1, by = NULL, t = NULL, cols = is.numeric, fill = NA, stubs = TRUE, keep.ids = TRUE, ...)
  L.data.frame(x, n, by, t, cols, fill, stubs, keep.ids, ...)

L.pdata.frame <- function(x, n = 1, cols = is.numeric, fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  ax <- attributes(x)
  nam <- ax[["names"]]
  index <- unclass(ax[["index"]])

  if(keep.ids) {
    gn <- which(nam %in% names(index))
    if(length(gn) && is.null(cols)) cols <- seq_along(unclass(x))[-gn]
  } else gn <- NULL

  if(length(index) > 2L) index <- list(finteraction(index[-length(index)]), index[[length(index)]])

  if(length(cols)) cols <- cols2int(cols, x, nam, FALSE)

  if(length(gn) && length(cols)) {
    class(x) <- NULL # Works for multiple lags !
    res <- c(x[gn], .Call(Cpp_flagleadl,x[cols],n,fill,fnlevels(index[[1L]]),index[[1L]],index[[2L]],stubs))
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!length(gn)) # could speed up ?
    return(.Call(Cpp_flagleadl,fcolsubset(x, cols),n,fill,fnlevels(index[[1L]]),index[[1L]],index[[2L]],stubs))
  .Call(Cpp_flagleadl,x,n,fill,fnlevels(index[[1L]]),index[[1L]],index[[2L]],stubs)
}


# Lead Operator
F <- function(x, n = 1, ...) UseMethod("F") # , x

F.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(flag.matrix(x, -n, g, t, fill, stubs, ...))
  flag.default(x, -n, g, t, fill, stubs, ...)
}

F.pseries <- function(x, n = 1, fill = NA, stubs = TRUE, ...)
  flag.pseries(x, -n, fill, stubs, ...)

F.matrix <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...)
  flag.matrix(x, -n, g, t, fill, stubs, ...)

F.grouped_df <- function(x, n = 1, t = NULL, fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {
  x <- x
  eval(substitute(flag.grouped_df(x, -n, t, fill, stubs, keep.ids, ...)))
}

F.data.frame <- function(x, n = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, stubs = TRUE, keep.ids = TRUE, ...)
  L.data.frame(x, -n, by, t, cols, fill, stubs, keep.ids, ...)

F.list <- F.data.frame

F.pdata.frame <- function(x, n = 1, cols = is.numeric, fill = NA, stubs = TRUE, keep.ids = TRUE, ...)
  L.pdata.frame(x, -n, cols, fill, stubs, keep.ids, ...)


