
# For principle innovations of this code see flag.R and flag.cpp

# Helper functions
checkld <- function(...) {
  if(any(names(list(...)) == "logdiff")) {
    warning("argument 'logdiff' was renamed to 'log'")
    TRUE
  } else FALSE
}
baselog <- base::log


fdiff <- function(x, n = 1, diff = 1, ...) UseMethod("fdiff") # , x

fdiff.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, log = FALSE, rho = 1, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(UseMethod("fdiff", unclass(x)))
  if(!missing(...)) if(checkld(...)) log <- list(...)[["logdiff"]] else unused_arg_action(match.call(), ...)
  if(log) x <- baselog(x)
  if(is.null(g)) return(.Call(Cpp_fdiffgrowth,x,n,diff,fill,0L,0L,NULL,G_t(t),1L+log,rho,stubs,1))
  g <- G_guo(g)
  .Call(Cpp_fdiffgrowth,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),1L+log,rho,stubs,1)
}

fdiff.pseries <- function(x, n = 1, diff = 1, fill = NA, log = FALSE, rho = 1, stubs = length(n) + length(diff) > 2L, shift = "time", ...) {
  if(!missing(...)) if(checkld(...)) log <- list(...)[["logdiff"]] else unused_arg_action(match.call(), ...)
  index <- uncl2pix(x)
  if(log) x <- baselog(x)
  g <- index[[1L]]
  t <- switch(shift, time = index[[2L]], row = NULL, stop("'shift' must be either 'time' or 'row'"))
  if(length(t) && !inherits(x, "indexed_series")) t <- plm_check_time(t)
  if(is.matrix(x))
    .Call(Cpp_fdiffgrowthm,x,n,diff,fill,fnlevels(g),g,NULL,t,1L+log,rho,stubs,1) else
      .Call(Cpp_fdiffgrowth,x,n,diff,fill,fnlevels(g),g,NULL,t,1L+log,rho,stubs,1)
}

fdiff.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, log = FALSE, rho = 1, stubs = length(n) + length(diff) > 2L, ...) {
  if(!missing(...)) if(checkld(...)) log <- list(...)[["logdiff"]] else unused_arg_action(match.call(), ...)
  if(log) x <- baselog(x)
  if(is.null(g)) return(.Call(Cpp_fdiffgrowthm,x,n,diff,fill,0L,0L,NULL,G_t(t),1L+log,rho,stubs,1))
  g <- G_guo(g)
  .Call(Cpp_fdiffgrowthm,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),1L+log,rho,stubs,1)
}

fdiff.grouped_df <- function(x, n = 1, diff = 1, t = NULL, fill = NA, log = FALSE, rho = 1, stubs = length(n) + length(diff) > 2L, keep.ids = TRUE, ...) {
  if(!missing(...)) if(checkld(...)) log <- list(...)[["logdiff"]] else unused_arg_action(match.call(), ...)
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
  cld <- function(x) if(log) fdapply(x, baselog) else x
  if(length(gn)) {
    ax <- attributes(x)
    res <- .Call(Cpp_fdiffgrowthl,cld(.subset(x, -gn)),n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),1L+log,rho,stubs,1)
    if(keep.ids) res <- c(.subset(x, gn), res)
    ax[["names"]] <- names(res)  # Works for multiple lags / differences !
    return(setAttributes(res, ax))
  }
  .Call(Cpp_fdiffgrowthl,cld(x),n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),1L+log,rho,stubs,1)
}

fdiff.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, log = FALSE, rho = 1, stubs = length(n) + length(diff) > 2L, ...) {
  if(!missing(...)) if(checkld(...)) log <- list(...)[["logdiff"]] else unused_arg_action(match.call(), ...)
  if(log) x <- fdapply(x, baselog)
  if(is.null(g)) return(.Call(Cpp_fdiffgrowthl,x,n,diff,fill,0L,0L,NULL,G_t(t),1L+log,rho,stubs,1))
  g <- G_guo(g)
  .Call(Cpp_fdiffgrowthl,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),1L+log,rho,stubs,1)
}

fdiff.list <- function(x, ...) fdiff.data.frame(x, ...)

fdiff.pdata.frame <- function(x, n = 1, diff = 1, fill = NA, log = FALSE, rho = 1, stubs = length(n) + length(diff) > 2L, shift = "time", ...) {
  if(!missing(...)) if(checkld(...)) log <- list(...)[["logdiff"]] else unused_arg_action(match.call(), ...)
  index <- uncl2pix(x)
  if(log) x <- fdapply(x, baselog)
  g <- index[[1L]]
  t <- switch(shift, time = index[[2L]], row = NULL, stop("'shift' must be either 'time' or 'row'"))
  if(length(t) && !inherits(x, "indexed_frame")) t <- plm_check_time(t)
  .Call(Cpp_fdiffgrowthl,x,n,diff,fill,fnlevels(g),g,NULL,t,1L+log,rho,stubs,1)
}




fgrowth <- function(x, n = 1, diff = 1, ...) UseMethod("fgrowth") # , x

fgrowth.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(UseMethod("fgrowth", unclass(x)))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(logdiff) x <- if(scale == 1) baselog(x) else scale * baselog(x)
  if(is.null(g)) return(.Call(Cpp_fdiffgrowth,x,n,diff,fill,0L,0L,NULL,G_t(t),4L-logdiff,scale,stubs,power))
  g <- G_guo(g)
  .Call(Cpp_fdiffgrowth,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),4L-logdiff,scale,stubs,power)
}

fgrowth.pseries <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = length(n) + length(diff) > 2L, shift = "time", ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- uncl2pix(x)
  if(logdiff) x <- if(scale == 1) baselog(x) else scale * baselog(x)
  g <- index[[1L]]
  t <- switch(shift, time = index[[2L]], row = NULL, stop("'shift' must be either 'time' or 'row'"))
  if(length(t) && !inherits(x, "indexed_series")) t <- plm_check_time(t)
  if(is.matrix(x))
    .Call(Cpp_fdiffgrowthm,x,n,diff,fill,fnlevels(g),g,NULL,t,4L-logdiff,scale,stubs,power) else
      .Call(Cpp_fdiffgrowth,x,n,diff,fill,fnlevels(g),g,NULL,t,4L-logdiff,scale,stubs,power)
}

fgrowth.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = length(n) + length(diff) > 2L, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(logdiff) x <- if(scale == 1) baselog(x) else scale * baselog(x)
  if(is.null(g)) return(.Call(Cpp_fdiffgrowthm,x,n,diff,fill,0L,0L,NULL,G_t(t),4L-logdiff,scale,stubs,power))
  g <- G_guo(g)
  .Call(Cpp_fdiffgrowthm,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),4L-logdiff,scale,stubs,power)
}

fgrowth.grouped_df <- function(x, n = 1, diff = 1, t = NULL, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = length(n) + length(diff) > 2L, keep.ids = TRUE, ...) {
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
  cld <- function(x) if(!logdiff) x else if(scale != 1) fdapply(x, function(y) scale * baselog(y)) else fdapply(x, baselog)
  if(length(gn)) {
    ax <- attributes(x)
    res <- .Call(Cpp_fdiffgrowthl,cld(.subset(x, -gn)),n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),4L-logdiff,scale,stubs,power)
    if(keep.ids) res <- c(.subset(x, gn), res)
    ax[["names"]] <- names(res)  # Works for multiple lags / differences !
    return(setAttributes(res, ax))
  }
  .Call(Cpp_fdiffgrowthl,cld(x),n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),4L-logdiff,scale,stubs,power)
}

fgrowth.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = length(n) + length(diff) > 2L, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(logdiff) x <- if(scale == 1) fdapply(x, baselog) else fdapply(x, function(y) scale * baselog(y))
  if(is.null(g)) return(.Call(Cpp_fdiffgrowthl,x,n,diff,fill,0L,0L,NULL,G_t(t),4L-logdiff,scale,stubs,power))
  g <- G_guo(g)
  .Call(Cpp_fdiffgrowthl,x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),4L-logdiff,scale,stubs,power)
}

fgrowth.list <- function(x, ...) fgrowth.data.frame(x, ...)

fgrowth.pdata.frame <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = length(n) + length(diff) > 2L, shift = "time", ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- uncl2pix(x)
  if(logdiff) x <- if(scale == 1) fdapply(x, baselog) else fdapply(x, function(y) scale * baselog(y))
  g <- index[[1L]]
  t <- switch(shift, time = index[[2L]], row = NULL, stop("'shift' must be either 'time' or 'row'"))
  if(length(t) && !inherits(x, "indexed_frame")) t <- plm_check_time(t)
  .Call(Cpp_fdiffgrowthl,x,n,diff,fill,fnlevels(g),g,NULL,t,4L-logdiff,scale,stubs,power)
}

# Operator data frame methods templates

DG_data_frame_template <- function(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, return = 1L, rho = 1, stubs = TRUE, keep.ids = TRUE, power = 1, ...) { # , message = 2L, power = 1

  if(!missing(...)) unused_arg_action(match.call(), ...)

  cld <- function(y) switch(return, y, fdapply(y, baselog), if(rho == 1) fdapply(y, baselog) else fdapply(y, function(k) rho * baselog(k)), y)

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
        cols <- cols2intrmgn(gn, cols, x)
      }
      by <- G_guo(if(length(gn) == 1L) x[[gn]] else x[gn])
      if(!keep.ids) gn <- NULL
    } else {
      gn <- NULL
      if(length(cols)) cols <- cols2int(cols, x, nam)
      by <- if(is.null(by)) list(0L, 0L, NULL) else G_guo(by)
    }

    if(is.call(t)) {
      tn <- ckmatch(all.vars(t), nam)
      t1 <- length(tn) == 1L
      t <- eval(if(t1) t[[2L]] else attr(terms.formula(t), "variables"), x, attr(t, ".Environment")) # if(t1) x[[tn]] else x[tn]
      cols <- if(is.null(cols)) seq_along(x)[-tn] else if(t1) cols[cols != tn] else fsetdiff(cols, tn)
      if(keep.ids) gn <- c(gn, tn)
    }

    res <- if(length(gn))
      c(x[gn], .Call(Cpp_fdiffgrowthl,cld(x[cols]),n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),return,rho,stubs,power)) else
        .Call(Cpp_fdiffgrowthl,cld(x[cols]),n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),return,rho,stubs,power)
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(length(cols)) { # Needs to be done like this, otherwise list-subsetting drops attributes !
    ax <- attributes(x)
    class(x) <- NULL
    x <- x[cols2int(cols, x, names(x), FALSE)]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }

  if(is.null(by)) return(.Call(Cpp_fdiffgrowthl,cld(x),n,diff,fill,0L,0L,NULL,G_t(t),return,rho,stubs,power))
  by <- G_guo(by)
  .Call(Cpp_fdiffgrowthl,cld(x),n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),return,rho,stubs,power)
}

DG_pdata_frame_template <- function(x, n = 1, diff = 1, cols = is.numeric, fill = NA, return = 1L, rho = 1, stubs = TRUE, shift = "time",
                          keep.ids = TRUE, power = 1, ...) {

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

  cld <- function(y) switch(return, y, fdapply(y, baselog), if(rho == 1) fdapply(y, baselog) else fdapply(y, function(k) rho * baselog(k)), y)

  if(length(gn) && length(cols)) {
    class(x) <- NULL # Works for multiple lags !
    res <- c(x[gn], .Call(Cpp_fdiffgrowthl,cld(x[cols]),n,diff,fill,fnlevels(g),g,NULL,t,return,rho,stubs,power))
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!length(gn)) # could speed up ?
    return(.Call(Cpp_fdiffgrowthl,cld(fcolsubset(x, cols)),n,diff,fill,fnlevels(g),g,NULL,t,return,rho,stubs,power))
  .Call(Cpp_fdiffgrowthl,cld(x),n,diff,fill,fnlevels(g),g,NULL,t,return,rho,stubs,power)
}

# Difference Operator (masks stats::D)  # use xt instead of by ?

# setGeneric("D")

D <- function(x, n = 1, diff = 1, ...) UseMethod("D") # , x

D.expression <- function(x, ...) if(missing(x)) stats::D(...) else stats::D(x, ...)
D.call <- function(x, ...) if(missing(x)) stats::D(...) else stats::D(x, ...)
D.name <- function(x, ...) if(missing(x)) stats::D(...) else stats::D(x, ...)

D.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, rho = 1, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fdiff.matrix(x, n, diff, g, t, fill, FALSE, rho, stubs, ...))
  fdiff.default(x, n, diff, g, t, fill, FALSE, rho, stubs, ...)
}

D.pseries <- function(x, n = 1, diff = 1, fill = NA, rho = 1, stubs = TRUE, shift = "time", ...)
  fdiff.pseries(x, n, diff, fill, FALSE, rho, stubs, shift, ...)

# setOldClass("pseries")
# setMethod("D", signature(expr = "pseries"), D.pseries)

D.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, rho = 1, stubs = TRUE, ...)
  fdiff.matrix(x, n, diff, g, t, fill, FALSE, rho, stubs, ...)

# setMethod("D", "matrix")

D.grouped_df <- function(x, n = 1, diff = 1, t = NULL, fill = NA, rho = 1, stubs = TRUE, keep.ids = TRUE, ...) {
  x <- x # because of piped calls -> "." is not in global environment ...
  eval(substitute(fdiff.grouped_df(x, n, diff, t, fill, FALSE, rho, stubs, keep.ids, ...)))
}

D.data.frame <- function(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, rho = 1, stubs = TRUE, keep.ids = TRUE, ...)
  DG_data_frame_template(x, n, diff, by, t, cols, fill, 1L, rho, stubs, keep.ids, ...)

D.list <- function(x, ...) D.data.frame(x, ...)

D.pdata.frame <- function(x, n = 1, diff = 1, cols = is.numeric, fill = NA, rho = 1, stubs = TRUE, shift = "time",
                          keep.ids = TRUE, ...)
  DG_pdata_frame_template(x, n, diff, cols, fill, 1L, rho, stubs, shift, keep.ids, ...)

# Log-Difference Operator

Dlog <- function(x, n = 1, diff = 1, ...) UseMethod("Dlog") # , x

Dlog.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, rho = 1, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fdiff.matrix(x, n, diff, g, t, fill, TRUE, rho, stubs, ...))
  fdiff.default(x, n, diff, g, t, fill, TRUE, rho, stubs, ...)
}

Dlog.pseries <- function(x, n = 1, diff = 1, fill = NA, rho = 1, stubs = TRUE, shift = "time", ...)
  fdiff.pseries(x, n, diff, fill, TRUE, rho, stubs, shift, ...)

Dlog.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, rho = 1, stubs = TRUE, ...)
  fdiff.matrix(x, n, diff, g, t, fill, TRUE, rho, stubs, ...)

Dlog.grouped_df <- function(x, n = 1, diff = 1, t = NULL, fill = NA, rho = 1, stubs = TRUE, keep.ids = TRUE, ...) {
  x <- x
  eval(substitute(fdiff.grouped_df(x, n, diff, t, fill, TRUE, rho, stubs, keep.ids, ...)))
}

Dlog.data.frame <- function(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, rho = 1, stubs = TRUE, keep.ids = TRUE, ...)
  DG_data_frame_template(x, n, diff, by, t, cols, fill, 2L, rho, stubs, keep.ids, ...)

Dlog.list <- function(x, ...) Dlog.data.frame(x, ...)

Dlog.pdata.frame <- function(x, n = 1, diff = 1, cols = is.numeric, fill = NA, rho = 1, stubs = TRUE, shift = "time",
                          keep.ids = TRUE, ...)
  DG_pdata_frame_template(x, n, diff, cols, fill, 2L, rho, stubs, shift, keep.ids, ...)


# Growth Operator

G <- function(x, n = 1, diff = 1, ...) UseMethod("G") # , x

G.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fgrowth.matrix(x, n, diff, g, t, fill, logdiff, scale, power, stubs, ...))
  fgrowth.default(x, n, diff, g, t, fill, logdiff, scale, power, stubs, ...)
}

G.pseries <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = TRUE, shift = "time", ...)
  fgrowth.pseries(x, n, diff, fill, logdiff, scale, power, stubs, shift, ...)

G.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = TRUE, ...)
  fgrowth.matrix(x, n, diff, g, t, fill, logdiff, scale, power, stubs, ...)

G.grouped_df <- function(x, n = 1, diff = 1, t = NULL, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = TRUE, keep.ids = TRUE, ...) {
  x <- x
  eval(substitute(fgrowth.grouped_df(x, n, diff, t, fill, logdiff, scale, power, stubs, keep.ids, ...)))
}

G.data.frame <- function(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = TRUE, keep.ids = TRUE, ...)
  DG_data_frame_template(x, n, diff, by, t, cols, fill, 4L-logdiff, scale, stubs, keep.ids, power, ...)

G.list <- function(x, ...) G.data.frame(x, ...)

G.pdata.frame <- function(x, n = 1, diff = 1, cols = is.numeric, fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = TRUE, shift = "time", keep.ids = TRUE, ...)
  DG_pdata_frame_template(x, n, diff, cols, fill, 4L-logdiff, scale, stubs, shift, keep.ids, power, ...)
