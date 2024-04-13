
ford <- function(x, g = NULL) {
  if(is.null(x)) return(NULL)
  if(!is.null(g)) {
    x <- c(if(is.atomic(g)) list(g) else if(is_GRP(g)) g[2L] else g,
           if(is.atomic(x)) list(x) else x, list(method = "radix"))
    return(do.call(order, x))
  }
  if(is.list(x)) return(do.call(order, c(x, list(method = "radix"))))
  if(length(x) < 1000L) .Call(C_radixsort, TRUE, FALSE, FALSE, FALSE, TRUE, pairlist(x)) else order(x, method = "radix")
}

fcumsum <- function(x, ...) UseMethod("fcumsum") # , x

fcumsum.default <- function(x, g = NULL, o = NULL, na.rm = .op[["na.rm"]], fill = FALSE, check.o = TRUE, ...) {
  # if(is.matrix(x) && !inherits(x, "matrix")) return(UseMethod("fcumsum", unclass(x)))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(length(o) && check.o) o <- ford(o, g)
  if(is.null(g)) return(.Call(C_fcumsum,x,0L,0L,o,na.rm,fill))
  g <- G_guo(g)
  .Call(C_fcumsum,x,g[[1L]],g[[2L]],o,na.rm,fill)
}

fcumsum.pseries <- function(x, na.rm = .op[["na.rm"]], fill = FALSE, shift = "time", ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- uncl2pix(x)
  g <- index[[1L]]
  o <- switch(shift, time = ford(index[[2L]], g), row = NULL, stop("'shift' must be either 'time' or 'row'"))
  if(is.matrix(x))
    .Call(C_fcumsumm,x,fnlevels(g),g,o,na.rm,fill) else
      .Call(C_fcumsum,x,fnlevels(g),g,o,na.rm,fill)
}

fcumsum.matrix <- function(x, g = NULL, o = NULL, na.rm = .op[["na.rm"]], fill = FALSE, check.o = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(length(o) && check.o) o <- ford(o, g)
  if(is.null(g)) return(.Call(C_fcumsumm,x,0L,0L,o,na.rm,fill))
  g <- G_guo(g)
  .Call(C_fcumsumm,x,g[[1L]],g[[2L]],o,na.rm,fill)
}

fcumsum.zoo <- function(x, ...) if(is.matrix(x)) fcumsum.matrix(x, ...) else fcumsum.default(x, ...)
fcumsum.units <- fcumsum.zoo

fcumsum.grouped_df <- function(x, o = NULL, na.rm = .op[["na.rm"]], fill = FALSE, check.o = TRUE, keep.ids = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  osym <- substitute(o)
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  if(!is.null(osym)) {
    o <- eval(osym, x, parent.frame())
    if(!anyNA(on <- match(all.vars(osym), nam))) {
      gn <- c(gn, on)
      if(anyDuplicated.default(gn)) stop("timevar coincides with grouping variables!")
    }
    if(check.o) o <- ford(o, g)
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

fcumsum.data.frame <- function(x, g = NULL, o = NULL, na.rm = .op[["na.rm"]], fill = FALSE, check.o = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(length(o) && check.o) o <- ford(o, g)
  if(is.null(g)) return(.Call(C_fcumsuml,x,0L,0L,o,na.rm,fill))
  g <- G_guo(g)
  .Call(C_fcumsuml,x,g[[1L]],g[[2L]],o,na.rm,fill)
}

fcumsum.list <- function(x, ...) fcumsum.data.frame(x, ...)

fcumsum.pdata.frame <- function(x, na.rm = .op[["na.rm"]], fill = FALSE, shift = "time", ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- uncl2pix(x)
  g <- index[[1L]]
  o <- switch(shift, time = ford(index[[2L]], g), row = NULL, stop("'shift' must be either 'time' or 'row'"))
  .Call(C_fcumsuml,x,fnlevels(g),g,o,na.rm,fill)
}
