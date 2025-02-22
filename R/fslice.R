fslice <- function(x, ..., n = 1, how = "first", order.by = NULL, na.rm = .op[["na.rm"]], with.ties = TRUE) {
  # handle grouping
  if (!missing(...)) {
    g <- GRP(fselect(x, ...), sort = FALSE, return.order = FALSE, call = FALSE)
  } else {
    g <- NULL
  }

  # convert a proportion to a number if applicable
  if (n < 1) {
    if (is.null(g)) {
      n <- max(1, round(n * fnrow(x)))
    } else {
      n <- pmax(1, round(n * g$group.sizes))
    }
  }

  # resolve values to order by
  if (!is.null(order.by)) {
    if (is.character(order.by) && order.by %in% names(x)) {
      order.by <- vec(fselect(x, order.by))
    }
    if (!is.vector(order.by) || length(order.by) != fnrow(x)) {
      stop("order.by must be a vector of the same length as the number of rows in x")
    }
  }

  # get the right function
  select_fn <- switch(how,
    "first" = function(idx) head(idx, n),
    "last" = function(idx) tail(idx, n),
    "min" = function(idx) {
      if (is.null(order.by)) stop("'order.by' must be provided for 'min'")
      values <- order.by[idx]
      indices <- radixorderv(values, decreasing = FALSE, na.last = !na.rm)
      if (with.ties) {
        threshold <- values[indices[n]]
        indices <- indices[values[indices] <= threshold]
      } else {
        indices <- head(indices, n)
      }
      idx[indices]
    },
    "max" = function(idx) {
      if (is.null(order.by)) stop("'order.by' must be provided for 'max'")
      values <- order.by[idx]
      indices <- radixorderv(values, decreasing = TRUE, na.last = !na.rm)
      if (with.ties) {
        threshold <- values[indices[n]]
        indices <- indices[values[indices] >= threshold]
      } else {
        indices <- head(indices, n)
      }
      idx[indices]
    },
    stop("Invalid 'how' argument: must be 'first', 'last', 'min', or 'max'")
  )

  indices <- seq_len(fnrow(x))
  if (is.null(g)) {
    selected <- select_fn(indices)
  } else {
    selected <- unlist(lapply(split(indices, g$group.id), select_fn))
  }
  fsubset(x, selected)
}
