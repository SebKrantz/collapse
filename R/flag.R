
flag <- function(x, n = 1, ...) UseMethod("flag") # , x

flag.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(UseMethod("flag", unclass(x)))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_flaglead,x,n,fill,0L,0L,G_t(t),stubs))
  g <- G_guo(g)
  .Call(Cpp_flaglead,x,n,fill,g[[1L]],g[[2L]],G_t(t),stubs)
}

flag.pseries <- function(x, n = 1, fill = NA, stubs = length(n) > 1L, shift = "time", ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- uncl2pix(x)
  g <- index[[1L]]
  t <- switch(shift, time = index[[2L]], row = NULL, stop("'shift' must be either 'time' or 'row'"))
  if(length(t) && !inherits(x, "indexed_series")) t <- plm_check_time(t)
  if(is.matrix(x))
  .Call(Cpp_flagleadm,x,n,fill,fnlevels(g),g,t,stubs) else
  .Call(Cpp_flaglead,x,n,fill,fnlevels(g),g,t,stubs)
}

flag.matrix <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = length(n) > 1L, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_flagleadm,x,n,fill,0L,0L,G_t(t),stubs))
  g <- G_guo(g)
  .Call(Cpp_flagleadm,x,n,fill,g[[1L]],g[[2L]],G_t(t),stubs)
}

flag.grouped_df <- function(x, n = 1, t = NULL, fill = NA, stubs = length(n) > 1L, keep.ids = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  tsym <- substitute(t)
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  if(!is.null(tsym)) {
    t <- eval(tsym, x, parent.frame())
    if(!anyNA(tn <- match(all.vars(tsym), nam))) {
      gn <- c(gn, tn)
      if(anyDuplicated.default(gn)) stop("timevar coincides with grouping variables!")
    }
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
  if(is.null(g)) return(.Call(Cpp_flagleadl,x,n,fill,0L,0L,G_t(t),stubs))
  g <- G_guo(g)
  .Call(Cpp_flagleadl,x,n,fill,g[[1L]],g[[2L]],G_t(t),stubs)
}

flag.list <- function(x, ...) flag.data.frame(x, ...)

flag.pdata.frame <- function(x, n = 1, fill = NA, stubs = length(n) > 1L, shift = "time", ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- uncl2pix(x)
  g <- index[[1L]]
  t <- switch(shift, time = index[[2L]], row = NULL, stop("'shift' must be either 'time' or 'row'"))
  if(length(t) && !inherits(x, "indexed_frame")) t <- plm_check_time(t)
  .Call(Cpp_flagleadl,x,n,fill,fnlevels(g),g,t,stubs)
}

# Lag Operator   # use xt instead of by ?
L <- function(x, n = 1, ...) UseMethod("L") # , x

L.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(flag.matrix(x, n, g, t, fill, stubs, ...))
  flag.default(x, n, g, t, fill, stubs, ...)
}

L.pseries <- function(x, n = 1, fill = NA, stubs = TRUE, shift = "time", ...)
  flag.pseries(x, n, fill, stubs, shift, ...)

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
        cols <- cols2intrmgn(gn, cols, x)
      }
      by <- G_guo(if(length(gn) == 1L) x[[gn]] else x[gn])
      if(!keep.ids) gn <- NULL
    } else {
      gn <- NULL
      if(length(cols)) cols <- cols2int(cols, x, nam)
      by <- if(is.null(by)) list(0L, 0L) else G_guo(by)
    }

    if(is.call(t)) {
      tn <- ckmatch(all.vars(t), nam, "Unknown variables:")
      t1 <- length(tn) == 1L
      t <- eval(if(t1) t[[2L]] else attr(terms.formula(t), "variables"), x, attr(t, ".Environment")) # if(t1) x[[tn]] else x[tn]
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

  if(is.null(by)) return(.Call(Cpp_flagleadl,x,n,fill,0L,0L,G_t(t),stubs))
  by <- G_guo(by)
  .Call(Cpp_flagleadl,x,n,fill,by[[1L]],by[[2L]],G_t(t),stubs)
}

L.list <- function(x, ...) L.data.frame(x, ...)

L.pdata.frame <- function(x, n = 1, cols = is.numeric, fill = NA, stubs = TRUE, shift = "time", keep.ids = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  ax <- attributes(x)
  nam <- ax[["names"]]
  index <- uncl2pix(x)
  cols_fun <- is.function(cols)

  if(cols_fun && identical(cols, is.numeric)) cols <- which(.Call(C_vtypes, x, 1L))
  else if(length(cols)) cols <- cols2int(cols, x, nam, FALSE)

  if(cols_fun || keep.ids) {
    gn <- which(nam %in% attr(index, "nam")) # Needed for 1 or 3+ index variables
    if(length(gn)) {
      if(cols_fun) cols <- fsetdiff(cols, gn)
      else if(is.null(cols)) cols <- seq_along(unclass(x))[-gn]
    }
    if(!keep.ids) gn <- NULL
  } else gn <- NULL

  g <- index[[1L]]
  t <- switch(shift, time = index[[2L]], row = NULL, stop("'shift' must be either 'time' or 'row'"))
  if(length(t) && !any(ax$class == "indexed_frame")) t <- plm_check_time(t)

  if(length(gn) && length(cols)) {
    class(x) <- NULL # Works for multiple lags !
    res <- c(x[gn], .Call(Cpp_flagleadl,x[cols],n,fill,fnlevels(g),g,t,stubs))
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!length(gn)) # could speed up ?
    return(.Call(Cpp_flagleadl,fcolsubset(x, cols),n,fill,fnlevels(g),g,t,stubs))
  .Call(Cpp_flagleadl,x,n,fill,fnlevels(g),g,t,stubs)
}


# Lead Operator
F <- function(x, n = 1, ...) UseMethod("F") # , x

F.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(flag.matrix(x, -n, g, t, fill, stubs, ...))
  flag.default(x, -n, g, t, fill, stubs, ...)
}

F.pseries <- function(x, n = 1, fill = NA, stubs = TRUE, shift = "time", ...)
  flag.pseries(x, -n, fill, stubs, shift, ...)

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

F.pdata.frame <- function(x, n = 1, cols = is.numeric, fill = NA, stubs = TRUE, shift = "time", keep.ids = TRUE, ...)
  L.pdata.frame(x, -n, cols, fill, stubs, shift, keep.ids, ...)


