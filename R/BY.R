
BY <- function(x, ...) UseMethod("BY")


BY.default <- function(x, g, FUN, ..., use.g.names = TRUE, sort = .op[["sort"]], reorder = TRUE,
                       expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                       return = c("same", "vector", "list")) {

  # If matrix, dispatch to matrix method
  if(is.matrix(x) && !inherits(x, "matrix")) return(UseMethod("BY", unclass(x)))
  if(!is.function(FUN)) FUN <- match.fun(FUN)

  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply
  simplify <- switch(return[1L], same = 1L, vector = 2L, list = 3L, stop("BY.default only supports same, vector and list output!"))
  g <- GRP(g, return.groups = use.g.names, sort = sort, call = FALSE)

  # Computing result: unsimplified
  if(!missing(...) && g[[1L]] > 1L && length(ln <- whichv(vlengths(dots <- list(...), FALSE), length(x)))) {
    asl <- lapply(dots[ln], gsplit, g)
    if(length(dots) > length(ln)) {
      mord <- dots[-ln]
      if(is.null(names(mord)) && is.null(names(asl))) warning("If some arguments have the same length as the data (vectors) while others have length 1 (scalars), please ensure that at least one of the two have keywords e.g. argname = value. This is because the latter are passed to the 'MoreArgs' argument of .mapply, and thus the order in which arguments are passed to the function might be different from your top-level call. In particular, .mapply will first pass the vector valued arguments followed by the scalar valued ones.")
    } else mord <- NULL
    res <- .mapply(FUN, c(list(gsplit(x, g)), asl), mord)
  } else res <- aplyfun(gsplit(x, g), FUN, ...)

  # Returning raw or wide result
  if(simplify == 3L || expand.wide) {
    if(use.g.names) names(res) <- GRPnames(g, FALSE)
    if(simplify == 3L) return(res)
    return(do.call(rbind, res))
  }

  # If using names and function also assigns names e.g. quantile()
  if(use.g.names && length(names(res[[1L]]))) {

    names(res) <- GRPnames(g, FALSE)
    res <- unlist(res, recursive = FALSE, use.names = TRUE)

    if(reorder && length(res) == length(g[[2L]]) && !isTRUE(g$ordered[2L]))
      warning("result is same length as x but the grouping is not sorted and the function used added names. Thus BY cannot decisively distinguish whether you are using a transformation function like scale() or a summary function like quantile() that computes a vector of statistics. The latter is assumed and the result is not reordered. To receive reordered output without constructed names set use.g.names = FALSE")

  } else { # Function does not assign names... or not using group names...

    res <- unlist(res, FALSE, FALSE)

    if(length(res) == g[[1L]]) {

      if(use.g.names) names(res) <- GRPnames(g, FALSE)

    } else if(length(res) == length(g[[2L]])) {

      if(reorder) res <- .Call(C_greorder, res, g)
      if(length(names(x)) && (reorder || isTRUE(g$ordered[2L]))) # Making sure we don't assign wrong names..
        names(res) <- names(x)

    }
  }

  if(simplify == 1L) return(copyMostAttributes(res, x)) # here needs to be copyMostAttributes... otherwise overwrites names
  res
}

# Experimental: But not really faster and also risky because vapply checks types and types may differ...
# copysplaplfun <- function(x, g, FUN, ...) {
#   sx <- gsplit(x, g)
#   if(length(sx) > 100000L && length(r1 <- match.fun(FUN)(sx[[1L]], ...)) == 1L)
#     return(copyMostAttributes(vapply(sx, FUN, r1, ..., USE.NAMES = FALSE), x))
#   copyMostAttributes(unlist(lapply(sx, FUN, ...), FALSE, FALSE), x)
# }

copysplaplfun <- function(x, g, FUN, ...) copyMostAttributes(unlist(lapply(gsplit(x, g), FUN, ...), FALSE, FALSE), x)
copysplmaplfun <- function(x, g, FUN, asl, mord) copyMostAttributes(unlist(.mapply(FUN, c(list(gsplit(x, g)), asl), mord), FALSE, FALSE), x)
splaplfun <- function(x, g, FUN, ...) unlist(lapply(gsplit(x, g), FUN, ...), FALSE, FALSE)
splmaplfun <- function(x, g, FUN, asl, mord) unlist(.mapply(FUN, c(list(gsplit(x, g)), asl), mord), FALSE, FALSE)

BY.data.frame <- function(x, g, FUN, ..., use.g.names = TRUE, sort = .op[["sort"]], reorder = TRUE,
                          expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                          return = c("same", "matrix", "data.frame", "list")) {

  if(!is.list(x)) stop("x needs to be a list")
  if(!is.function(FUN)) FUN <- match.fun(FUN)

  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply
  return <- switch(return[1L], same = 1L, matrix = 3L, data.frame = 2L, list = 0L,
                   stop("Unknown return option!"))
  g <- GRP(g, return.groups = use.g.names, sort = sort, call = FALSE)
  n <- length(g[[2L]])

  if(!missing(...) && g[[1L]] > 1L && length(ln <- whichv(vlengths(dots <- list(...), FALSE), n))) {
    asl <- lapply(dots[ln], gsplit, g)
    if(length(dots) > length(ln)) {
      mord <- dots[-ln]
      if(is.null(names(mord)) && is.null(names(asl))) warning("If some arguments have the same length as the data (vectors) while others have length 1 (scalars), please ensure that at least one of the two have keywords e.g. argname = value. This is because the latter are passed to the 'MoreArgs' argument of .mapply, and thus the order in which arguments are passed to the function might be different from your top-level call. In particular, .mapply will first pass the vector valued arguments followed by the scalar valued ones.")
    } else mord <- NULL
    multi <- TRUE
  } else multi <- FALSE

  # Just plain list output
  if(return == 0L) {
    if(multi) {
      if(expand.wide) return(aplyfun(x, function(y) do.call(rbind, .mapply(FUN, c(list(gsplit(y, g, use.g.names)), asl), mord))))
      return(aplyfun(x, function(y) .mapply(FUN, c(list(gsplit(y, g, use.g.names)), asl), mord)))
    }
    if(expand.wide) return(aplyfun(x, function(y) do.call(rbind, lapply(gsplit(y, g, use.g.names), FUN, ...))))
    return(aplyfun(x, function(y) lapply(gsplit(y, g, use.g.names), FUN, ...)))
  }

  ax <- attributes(x)

  # Wider output (for multiple summary statistics like quantile())
  if(expand.wide) {
    if(return < 3L) { # Return a data.frame
      splitfun <- if(multi) function(y) .Call(Cpp_mctl, do.call(rbind, .mapply(FUN, c(list(gsplit(y, g)), asl), mord)), TRUE, 0L) else
        function(y) .Call(Cpp_mctl, do.call(rbind, lapply(gsplit(y, g), FUN, ...)), TRUE, 0L)
      res <- unlist(aplyfun(x, splitfun), recursive = FALSE, use.names = TRUE)
      if(return == 1L) {
        isDTl <- inherits(x, "data.table")
        ax[["names"]] <- names(res)
        ax[["row.names"]] <- if(use.g.names && !isDTl && length(gn <- GRPnames(g))) gn else
          .set_row_names(.Call(C_fnrow, res))
      } else {
        isDTl <- FALSE
        ax <- list(names = names(res),
                   row.names = if(use.g.names && length(gn <- GRPnames(g))) gn else .set_row_names(.Call(C_fnrow, res)),
                   class = "data.frame")
      }
      return(condalcSA(res, ax, isDTl))
    } else { # Return a matrix
      attributes(x) <- NULL
      splitfun <- if(multi) function(y) do.call(rbind, .mapply(FUN, c(list(gsplit(y, g)), asl), mord)) else
        function(y) do.call(rbind, lapply(gsplit(y, g), FUN, ...))
      res <- do.call(cbind, aplyfun(x, splitfun))
      cn <- dimnames(res)[[2L]]
      namr <- rep(ax[["names"]], each = ncol(res)/length(x))
      dimnames(res) <- list(if(use.g.names) GRPnames(g) else NULL,
                            if(length(cn)) paste(namr, cn, sep = ".") else namr)
      return(res)
    }
  }

  # No expand wide (classical result)
  matl <- return == 3L
  isDTl <- !matl && return != 2L && inherits(x, "data.table") # is data table and return data.table
  rownam <- ax[["row.names"]]
  attributes(x) <- NULL

  # Returning plain data frame
  if(return == 2L) ax <- list(names = ax[["names"]], row.names = rownam, class = "data.frame")

  # Using group names...
  if(use.g.names && !isDTl && !is.null(g$groups)) {
    res <- vector("list", length(x))
    res1 <- if(multi) .mapply(FUN, c(list(gsplit(x[[1L]], g)), asl), mord) else lapply(gsplit(x[[1L]], g), FUN, ...)

    if(length(names(res1[[1L]]))) { # We apply a function that assigns names (e.g. quantile())
      names(res1) <- GRPnames(g)
      res1 <- unlist(res1, recursive = FALSE, use.names = TRUE)
      rn <- names(res1)
      names(res1) <- NULL
      if(reorder && length(res1) == n && !isTRUE(g$ordered[2L])) {
        warning("nrow(result) is same as nrow(x) but the grouping is not sorted and the function used added names. Thus BY cannot decisively distinguish whether you are using a transformation function like scale() or a summary function like quantile() that computes a vector of statistics. The latter is assumed and the result is not reordered. To receive reordered output without constructed names set use.g.names = FALSE")
        reorder <- FALSE
      }
    } else { # function doesn't assign names, different options.
      res1 <- unlist(res1, FALSE, FALSE)
      if(length(res1) == g[[1L]]) rn <- GRPnames(g)
      else if(matl) {
        rn <- if(length(res1) == n && is.character(rownam) && rownam[1L] != "1" && (reorder || isTRUE(g$ordered[2L]))) rownam else NULL
      } else {  # Important to check lenth(rn) below (simply keeps ax[["row.names"]])
        rn <- if(length(res1) != n || !(reorder || isTRUE(g$ordered[2L]))) .set_row_names(length(res1)) else NULL
      }
    }
    # Finish computing results...
    if(matl) {
      res[[1L]] <- res1
      if(length(res) > 1L) res[-1L] <- if(multi) aplyfun(x[-1L], splmaplfun, g, FUN, asl, mord) else aplyfun(x[-1L], splaplfun, g, FUN, ...)
      res <- do.call(cbind, res)
    } else {
      res[[1L]] <- copyMostAttributes(res1, x[[1L]])
      if(length(res) > 1L) res[-1L] <- if(multi) aplyfun(x[-1L], copysplmaplfun, g, FUN, asl, mord) else aplyfun(x[-1L], copysplaplfun, g, FUN, ...)
    }

  # Not using group names...
  } else {
    if(matl) {
      res <- if(multi) do.call(cbind, aplyfun(x, splmaplfun, g, FUN, asl, mord)) else do.call(cbind, aplyfun(x, splaplfun, g, FUN, ...))
      rn <- if(nrow(res) == n && is.character(rownam) && rownam[1L] != "1" && (reorder || isTRUE(g$ordered[2L]))) rownam else NULL
    } else {
      res <- if(multi) aplyfun(x, copysplmaplfun, g, FUN, asl, mord) else aplyfun(x, copysplaplfun, g, FUN, ...) # isDTL ? -> Not needed as data.tables cannot have character row-names anyway.
      rn <- if(.Call(C_fnrow, res) != n || !(reorder || isTRUE(g$ordered[2L]))) .set_row_names(.Call(C_fnrow, res)) else NULL
    }
  }

  # reorder result if necessary, without dimnames...
  if(reorder && fnrow(res) == n && !isTRUE(g$ordered[2L])) {
    ind <- .Call(C_greorder, seq_len(n), g)
    res <- if(matl) res[ind, , drop = FALSE] else .Call(C_subsetDT, res, ind, seq_along(res), FALSE)
  }

  if(matl) {
    dimnames(res) <- list(rn, ax[["names"]])
    return(res)
  }

  if(length(rn)) ax[["row.names"]] <- rn
  return(condalcSA(res, ax, isDTl))
}

BY.list <- function(x, ...) BY.data.frame(x, ...)

BY.matrix <- function(x, g, FUN, ..., use.g.names = TRUE, sort = .op[["sort"]], reorder = TRUE,
                      expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                      return = c("same", "matrix", "data.frame", "list")) {

  if(!is.matrix(x)) stop("x needs to be a matrix")
  if(!is.function(FUN)) FUN <- match.fun(FUN)

  aplyfun <- if(parallel) function(...) parallel::mclapply(..., mc.cores = mc.cores) else lapply
  return <- switch(return[1L], same = 3L, matrix = 2L, data.frame = 1L, list = 0L,
                   stop("Unknown return option!"))
  g <- GRP(g, return.groups = use.g.names, sort = sort, call = FALSE)
  n <- nrow(x)

  if(!missing(...) && g[[1L]] > 1L && length(ln <- whichv(vlengths(dots <- list(...), FALSE), n))) {
    asl <- lapply(dots[ln], gsplit, g)
    if(length(dots) > length(ln)) {
      mord <- dots[-ln]
      if(is.null(names(mord)) && is.null(names(asl))) warning("If some arguments have the same length as the data (vectors) while others have length 1 (scalars), please ensure that at least one of the two have keywords e.g. argname = value. This is because the latter are passed to the 'MoreArgs' argument of .mapply, and thus the order in which arguments are passed to the function might be different from your top-level call. In particular, .mapply will first pass the vector valued arguments followed by the scalar valued ones.")
    } else mord <- NULL
    multi <- TRUE
  } else multi <- FALSE

  # Just plain list output
  if(return == 0L) {
    xln <- .Call(Cpp_mctl, x, TRUE, 0L) # Named list from matrix
    if(multi) {
      if(expand.wide) return(aplyfun(xln, function(y) do.call(rbind, .mapply(FUN, c(list(gsplit(y, g, use.g.names)), asl), mord))))
      return(aplyfun(xln, function(y) .mapply(FUN, c(list(gsplit(y, g, use.g.names)), asl), mord)))
    }
    if(expand.wide) return(aplyfun(xln, function(y) do.call(rbind, lapply(gsplit(y, g, use.g.names), FUN, ...))))
    return(aplyfun(xln, function(y) lapply(gsplit(y, g, use.g.names), FUN, ...)))
  }

  # Wider output (for multiple summary statistics like quantile())
  if(expand.wide) {
    if(return == 1L) { # Return data frame
      splitfun <- if(multi) function(y) .Call(Cpp_mctl, do.call(rbind, .mapply(FUN, c(list(gsplit(y, g)), asl), mord)), TRUE, 0L) else
        function(y) .Call(Cpp_mctl, do.call(rbind, lapply(gsplit(y, g), FUN, ...)), TRUE, 0L)
      res <- unlist(aplyfun(.Call(Cpp_mctl, x, TRUE, 0L), splitfun), recursive = FALSE, use.names = TRUE)
      ax <- list(names = names(res),
                 row.names = if(use.g.names && length(gn <- GRPnames(g))) gn else .set_row_names(.Call(C_fnrow, res)),
                 class = "data.frame")
    } else { # Return a matrix
      splitfun2 <- if(multi) function(y) do.call(rbind, .mapply(FUN, c(list(gsplit(y, g)), asl), mord)) else
        function(y) do.call(rbind, lapply(gsplit(y, g), FUN, ...))
      res <- do.call(cbind, aplyfun(.Call(Cpp_mctl, x, FALSE, 0L), splitfun2))
      cn <- dimnames(res)[[2L]]
      namr <- rep(dimnames(x)[[2L]], each = ncol(res)/ncol(x))
      dn <- list(if(use.g.names) GRPnames(g) else NULL,
                 if(length(cn)) paste(namr, cn, sep = ".") else namr)
      if(return == 2L) return(`dimnames<-`(res, dn))
      ax <- attributes(x)
      ax[["dim"]] <- dim(res)
      ax[["dimnames"]] <- dn
    }
    return(setAttributes(res, ax))
  }

  dn <- dimnames(x)
  matl <- return > 1L
  xl <- .Call(Cpp_mctl, x, FALSE, 0L) # Plain list from matrix columns

  # No expand wide (classical result)
  if(use.g.names && !is.null(g$groups)) {
    res <- vector("list", length(xl))
    res1 <- if(multi) .mapply(FUN, c(list(gsplit(xl[[1L]], g)), asl), mord) else lapply(gsplit(xl[[1L]], g), FUN, ...)

    if(length(names(res1[[1L]]))) { # We apply a function that assigns names (e.g. quantile())
      names(res1) <- GRPnames(g)
      res1 <- unlist(res1, recursive = FALSE, use.names = TRUE)
      rn <- names(res1)
      names(res1) <- NULL
      if(reorder && length(res1) == n && !isTRUE(g$ordered[2L])) {
        warning("nrow(result) is same as nrow(x) but the grouping is not sorted and the function used added names. Thus BY cannot decisively distinguish whether you are using a transformation function like scale() or a summary function like quantile() that computes a vector of statistics. The latter is assumed and the result is not reordered. To receive reordered output without constructed names set use.g.names = FALSE")
        reorder <- FALSE
      }
    } else { # function doesn't assign names, different options.
      res1 <- unlist(res1, FALSE, FALSE)
      rn <- if(length(res1) == g[[1L]]) GRPnames(g) else
            if(length(res1) == n && (reorder || isTRUE(g$ordered[2L]))) dn[[1L]] else NULL
    }

    # Finish computing results...
    res[[1L]] <- res1
    if(length(res) > 1L) res[-1L] <- if(multi) aplyfun(xl[-1L], splmaplfun, g, FUN, asl, mord) else aplyfun(xl[-1L], splaplfun, g, FUN, ...)

    if(matl) { # Return a matrix
      res <- do.call(cbind, res)
      dn <- list(rn, dn[[2L]])
    }

  } else { # Not using group names

    res <- if(multi) aplyfun(xl, splmaplfun, g, FUN, asl, mord) else aplyfun(xl, splaplfun, g, FUN, ...)

    if(matl) { # Return matrix

      res <- do.call(cbind, res)
      if(nrow(res) != n || !(reorder || isTRUE(g$ordered[2L]))) dn <- list(NULL, dn[[2L]])

    } else { # Return data frame
      rn <- if(.Call(C_fnrow, res) == n && (reorder || isTRUE(g$ordered[2L]))) dn[[1L]] else NULL
    }

  }

  # reorder result if necessary, without dimnames...
  if(reorder && fnrow(res) == n && !isTRUE(g$ordered[2L])) {
    ind <- .Call(C_greorder, seq_len(n), g)
    res <- if(matl) res[ind, , drop = FALSE] else .Call(C_subsetDT, res, ind, seq_along(res), FALSE)
  }

  if(matl) {
    if(return == 2L) return(`dimnames<-`(res, dn))
    ax <- attributes(x)
    ax[["dim"]] <- dim(res)
    ax[["dimnames"]] <- dn
  } else { # Returning a data.frame
    ax <- list(names = dn[[2L]],
               row.names = if(length(rn)) rn else .set_row_names(.Call(C_fnrow, res)),
               class = "data.frame")
  }

  return(setAttributes(res, ax))
}

BY.grouped_df <- function(x, FUN, ..., reorder = TRUE, keep.group_vars = TRUE, use.g.names = FALSE) {
  g <- GRP.grouped_df(x, call = FALSE)
  gn <- which(attr(x, "names") %in% g[[5L]])

  res <- BY.data.frame(if(length(gn)) fcolsubset(x, -gn) else x, g, FUN, ..., reorder = reorder, use.g.names = use.g.names)

  # Other return options
  if(!is.data.frame(res)) return(res)

  n <- fnrow(res)
  same_size <- n == fnrow(x)
  if(!same_size && is.null(g[[4L]])) keep.group_vars <- FALSE

  # Not preserving grouping variable or same size and no grouping variables: return appropriate object
  if(!keep.group_vars || (same_size && length(gn) == 0L)) return(if(same_size && (reorder || isTRUE(g$ordered[2L]))) res else fungroup(res))

  # If same size, with grouping variables...
  if(same_size) {
    if(!(reorder || isTRUE(g$ordered[2L]))) return(fungroup(res))
    ar <- attributes(res)
    ar[["names"]] <- c(g[[5L]], ar[["names"]])
    return(condalcSA(c(.subset(x, gn), res), ar, any(ar$class == "data.table")))
  }

  # If other size or no groups
  if(n != g[[1L]]) {
    if(is.null(g[[4L]])) return(fungroup(res))
    len <- n / g[[1L]]
    if(len != as.integer(len)) stop("length of output (", n, ") is not a multiple of the number of groups: ", g[[1L]])
    g[[4L]] <- .Call(C_subsetDT, g[[4L]], rep(seq_len(g[[1L]]), each = len), seq_along(g[[5L]]), FALSE)
  }

  # Aggregation
  ar <- attributes(fungroup2(res, oldClass(res)))
  ar[["names"]] <- c(g[[5L]], ar[["names"]])
  condalcSA(c(g[[4L]], res), ar, any(ar$class == "data.table"))
}


