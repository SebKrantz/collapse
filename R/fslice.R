

fslice <- function(x, ..., n = 1, how = "first", order.by = NULL,
                   na.rm = .op[["na.rm"]], sort = FALSE, with.ties = FALSE) {

  # handle grouping
  if(!missing(...)) {
    g <- GRP.default(if(is.list(x)) fselect(unclass(x), ...) else list(...), sort = sort, return.groups = FALSE, return.order = sort, call = FALSE)
  } else if(is.list(x) && inherits(x, "grouped_df")) {
    g <- GRP.grouped_df(x, return.groups = FALSE, call = FALSE)
    x <- fungroup2(x, oldClass(x))
  } else g <- NULL

  # resolve values to order by
  if(switch(how, min = TRUE, max = TRUE, FALSE)) {
    if(is.list(x)) order.by <- eval(substitute(order.by), x, parent.frame())
    if(is.character(order.by) && length(order.by) == 1L && anyv(attr(x, "names"), order.by))
      order.by <- .subset2(x, order.by)
    if(length(order.by) != fnrow(x)) stop("order.by must be a numeric vector of the same length as the number of rows in x, or the name of a column in x.")
  }

  fslice_core(x, g, n, how, order.by, na.rm, with.ties, sort)
}

fslicev <- function(x, cols = NULL, n = 1, how = "first", order.by = NULL,
                   na.rm = .op[["na.rm"]], sort = FALSE, with.ties = FALSE, ...) {

  # handle grouping
  if(!is.null(cols)) {
    cond <- is.list(cols) || is.atomic(x)
    g <- GRP.default(if(cond) cols else x,
                     by = if(cond) NULL else cols,
                     sort = sort, return.groups = FALSE, return.order = sort, call = FALSE, ...)
  } else if(is.list(x) && inherits(x, "grouped_df")) {
    g <- GRP.grouped_df(x, return.groups = FALSE, call = FALSE)
    x <- fungroup2(x, oldClass(x))
  } else g <- NULL

  # resolve values to order by
  if(switch(how, min = TRUE, max = TRUE, FALSE)) {
    if(is.character(order.by) && length(order.by) == 1L && anyv(attr(x, "names"), order.by))
      order.by <- .subset2(x, order.by)
    if(length(order.by) != fnrow(x)) stop("order.by must be a numeric vector of the same length as the number of rows in x, or the name of a column in x.")
  }

  fslice_core(x, g, n, how, order.by, na.rm, with.ties, sort)
}


fslice_core <- function(x, g, n, how, order.by, na.rm, with.ties, sort) {

  # convert a proportion to a number if applicable
  if(n < 1) n <- if(is.null(g)) max(1L, as.integer(round(n * fnrow(x)))) else max(1L, as.integer(round(n * fnrow(x)/g[[1L]])))
  if(n > 1 && with.ties) stop("with.ties = TRUE is only supported for n = 1")

  if(is.null(g)) {
    ind <- switch(how,
      first = 1:n,
      last = (fnrow(x)-n+1L):fnrow(x),
      min = if(n > 1) radixorderv(order.by, decreasing = FALSE, na.last = na.rm)[1:n] else if(with.ties) order.by %==% fmin.default(order.by, na.rm = na.rm) else which.min(order.by),
      max = if(n > 1) radixorderv(order.by, decreasing = TRUE, na.last = na.rm)[1:n] else if(with.ties) order.by %==% fmax.default(order.by, na.rm = na.rm) else which.max(order.by),
      stop("Unknown 'how' option: ", how)
    )
    return(ss(x, ind, check = FALSE))
  }

  if(n == 1) {
    if(with.ties && sort) warning("sorting with ties is currently not supported")
    return(switch(how,
      first = condalc(ffirst(x, g, na.rm = FALSE), inherits(x, "data.table")),
      last = condalc(flast(x, g, na.rm = FALSE), inherits(x, "data.table")),
      # TODO: sort with ties?
      min = if(with.ties) ss(x, order.by %==% fmin(order.by, g, TRA = "fill", na.rm = na.rm, use.g.names = FALSE), check = FALSE) else
            ss(x, .Call(C_gwhich_first, order.by, g, fmin.default(order.by, g, na.rm = na.rm, use.g.names = FALSE)), check = FALSE),
      max = if(with.ties) ss(x, order.by %==% fmax(order.by, g, TRA = "fill", na.rm = na.rm, use.g.names = FALSE), check = FALSE) else
            ss(x, .Call(C_gwhich_first, order.by, g, fmax.default(order.by, g, na.rm = na.rm, use.g.names = FALSE)), check = FALSE),
      stop("Unknown 'how' option: ", how)
  ))
  }

  ind <- switch(how,
      first = .Call(C_gslice_multi, g, g$order, n, TRUE), # g$order is NULL if sort = FALSE
      last = .Call(C_gslice_multi, g, g$order, n, FALSE), # g$order is NULL if sort = FALSE
      min = .Call(C_gslice_multi, g, radixorder(g$group.id, order.by, decreasing = FALSE, na.last = na.rm), n, TRUE),
      max = .Call(C_gslice_multi, g, radixorder(g$group.id, order.by, decreasing = c(FALSE, TRUE), na.last = na.rm), n, TRUE),
      stop("Unknown 'how' option: ", how)
    )

  return(ss(x, ind, check = FALSE))
}
