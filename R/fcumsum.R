fcumsum <- function(x, ...) UseMethod("fcumsum") # , x

fcumsum.default <- function(x, g = NULL, o = NULL, na.rm = TRUE, fill = FALSE, check.o = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(UseMethod("fcumsum", unclass(x)))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(!is.null(o) && check.o) o <- ford(o, g)
  if(is.null(g)) return(.Call(C_fcumsum,x,0L,0L,o,na.rm,fill))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) nl <- fnlevels(g) else {
      g <- qG(g, na.exclude = FALSE)
      nl <- attr(g, "N.groups")
    }
    return(.Call(C_fcumsum,x,nl,g,o,na.rm,fill))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(C_fcumsum,x,g[[1L]],g[[2L]],o,na.rm,fill)
}

fcumsum.pseries <- function(x, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- unclass(attr(x, "index"))
  g <- if(length(index) > 2L) finteraction(index[-length(index)]) else index[[1L]]
  o <- ford(index[[length(index)]], g)
  if(is.matrix(x))
    .Call(C_fcumsumm,x,fnlevels(g),g,o,na.rm,fill) else
      .Call(C_fcumsum,x,fnlevels(g),g,o,na.rm,fill)
}

fcumsum.matrix <- function(x, g = NULL, o = NULL, na.rm = TRUE, fill = FALSE, check.o = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(!is.null(o) && check.o) o <- ford(o, g)
  if(is.null(g)) return(.Call(C_fcumsumm,x,0L,0L,o,na.rm,fill))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) nl <- fnlevels(g) else {
      g <- qG(g, na.exclude = FALSE)
      nl <- attr(g, "N.groups")
    }
    return(.Call(C_fcumsumm,x,nl,g,o,na.rm,fill))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(C_fcumsumm,x,g[[1L]],g[[2L]],o,na.rm,fill)
}

fcumsum.grouped_df <- function(x, o = NULL, na.rm = TRUE, fill = FALSE, check.o = TRUE, keep.ids = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  osym <- all.vars(substitute(o))
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  if(length(osym) && !anyNA(on <- match(osym, nam))) {
    if(length(on) == 1L) {
      if(any(gn == on)) stop("ordervar coincides with grouping variables!")
      o <- .subset2(x, on)
    } else {
      if(any(gn %in% on)) stop("ordervar coincides with grouping variables!")
      o <- .subset(x, on)
    }
    if(check.o) o <- ford(o, g)
    gn <- c(gn, on)
  }
  if(length(gn)) {
    ax <- attributes(x)
    res <- .Call(C_fcumsuml,.subset(x,-gn),g[[1L]],g[[2L]],o,na.rm,fill)
    if(keep.ids) res <- c(.subset(x, gn), res)
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  }
  .Call(C_fcumsuml,x,g[[1L]],g[[2L]],o,na.rm,fill)
}

fcumsum.data.frame <- function(x, g = NULL, o = NULL, na.rm = TRUE, fill = FALSE, check.o = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(!is.null(o) && check.o) o <- ford(o, g)
  if(is.null(g)) return(.Call(C_fcumsuml,x,0L,0L,o,na.rm,fill))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) nl <- fnlevels(g) else {
      g <- qG(g, na.exclude = FALSE)
      nl <- attr(g, "N.groups")
    }
    return(.Call(C_fcumsuml,x,nl,g,o,na.rm,fill))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(C_fcumsuml,x,g[[1L]],g[[2L]],o,na.rm,fill)
}

fcumsum.list <- function(x, g = NULL, o = NULL, na.rm = TRUE, fill = FALSE, check.o = TRUE, ...)
  fcumsum.data.frame(x, g, o, na.rm, fill, check.o, ...)

fcumsum.pdata.frame <- function(x, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- unclass(attr(x, "index"))
  g <- if(length(index) > 2L) finteraction(index[-length(index)]) else index[[1L]]
  o <- ford(index[[length(index)]], g)
  .Call(C_fcumsuml,x,fnlevels(g),g,o,na.rm,fill)
}
