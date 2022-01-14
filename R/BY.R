

BY <- function(x, ...) UseMethod("BY")

BY.default <- function(x, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                       expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                       return = c("same", "vector", "list")) {
  if(!is.atomic(x)) stop("x needs to be an atomic vector") # redundant ?
  if(is.matrix(x) && !inherits(x, "matrix")) return(UseMethod("BY", unclass(x)))
  if(!(is.function(FUN) || is.character(FUN))) stop("FUN needs to be a function")
  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply
  simplify <- switch(return[1L], same = 1L, vector = 2L, list = 3L, stop("BY.default only supports same, vector and list output!"))
  g <- GRP(g, return.groups = use.g.names, sort = sort, call = FALSE)
  res <- aplyfun(gsplit(x, g), FUN, ...)
  if(use.g.names) names(res) <- GRPnames(g, FALSE)
  if(simplify == 3L) return(res)
  if(expand.wide) return(do.call(rbind, res))
  if(use.g.names) {
    res <- unlist(res, recursive = FALSE)
    if(simplify == 1L) return(copyMostAttributes(res, x)) # here needs to be copyMostAttributes... otherwise overwrites names
  } else {
    if(simplify == 1L) {
      res <- unlist(res, FALSE, FALSE)
      if(length(res) == length(x) && typeof(res) == typeof(x) && isTRUE(g$ordered[2L])) return(duplAttributes(res, x))
      return(copyMostAttributes(res, x))
    }
    # If we return a vector and do not use group names but a function like quantile(), we may still replicate the names given by that function...
    ll <- length(res)
    nr1 <- names(res[[1L]])
    res <- unlist(res, FALSE, FALSE)
    if(length(res) != ll && length(nr1) && length(res) == length(nr1)*ll)
      names(res) <- rep(nr1, ll)
  }
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
splaplfun <- function(x, g, FUN, ...) unlist(lapply(gsplit(x, g), FUN, ...), FALSE, FALSE)


BY.data.frame <- function(x, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                          expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                          return = c("same", "matrix", "data.frame", "list")) {
  if(!is.list(x)) stop("x needs to be a list")
  if(!(is.function(FUN) || is.character(FUN))) stop("FUN needs to be a function")
  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply
  return <- switch(return[1L], same = 1L, matrix = 3L, data.frame = 2L, list = 0L,
                   stop("Unknown return option!"))
  g <- GRP(g, return.groups = use.g.names, sort = sort, call = FALSE)
  if(return != 0L) {
    ax <- attributes(x)
    if(expand.wide) {
      if(return < 3L) { # Return a data.frame
        splitfun <- function(y) .Call(Cpp_mctl, do.call(rbind, lapply(gsplit(y, g), FUN, ...)), TRUE, 0L)
        res <- unlist(aplyfun(x, splitfun), recursive = FALSE, use.names = TRUE)
        if(return == 1L) {
          isDTl <- inherits(x, "data.table")
          ax[["names"]] <- names(res)
          ax[["row.names"]] <- if(use.g.names && !inherits(x, "data.table") && length(gn <- GRPnames(g))) gn else
            .set_row_names(length(res[[1L]]))
        } else {
          isDTl <- FALSE
          ax <- list(names = names(res),
                     row.names = if(use.g.names && length(gn <- GRPnames(g))) gn else .set_row_names(length(res[[1L]])),
                     class = "data.frame")
        }
      } else { # Return a matrix
        attributes(x) <- NULL
        splitfun <- function(y) do.call(rbind, lapply(gsplit(y, g), FUN, ...))
        res <- do.call(cbind, aplyfun(x, splitfun))
        cn <- dimnames(res)[[2L]]
        namr <- rep(ax[["names"]], each = ncol(res)/length(x))
        dimnames(res) <- list(if(use.g.names) GRPnames(g) else NULL,
                              if(length(cn)) paste(namr, cn, sep = ".") else namr)
        return(res)
      }
    } else { # No expand wide (classical result)
      matl <- return == 3L
      isDTl <- !matl && return != 2L && inherits(x, "data.table")
      attributes(x) <- NULL
      if(return == 2L) ax <- list(names = ax[["names"]], row.names = ax[["row.names"]], class = "data.frame")
      if(use.g.names && !isDTl && length(gn <- GRPnames(g))) { # Using names...
        res <- vector("list", length(x))
        res1 <- lapply(gsplit(x[[1L]], g), FUN, ...)
        names(res1) <- gn
        res[[1L]] <- unlist(res1, FALSE, TRUE)
        namres1 <- names(res[[1L]])
        if(matl) dn <- list(namres1, ax[["names"]]) else
          if(length(namres1)) ax[["row.names"]] <- namres1 else
          if(length(res[[1L]]) != length(x[[1L]])) ax[["row.names"]] <- .set_row_names(length(res[[1L]]))
        if(length(namres1)) names(res[[1L]]) <- NULL
        if(matl) {
          if(length(res) > 1L) res[-1L] <- aplyfun(x[-1L], splaplfun, g, FUN, ...)
          res <- do.call(cbind, res)
          dimnames(res) <- dn
          return(res)
        } else {
          copyMostAttributes(res[[1L]], x[[1L]])
          if(length(res) > 1L) res[-1L] <- aplyfun(x[-1L], copysplaplfun, g, FUN, ...)
        }
      } else { # Not using names...
        if(matl) {
          res <- do.call(cbind, aplyfun(x, splaplfun, g, FUN, ...))
          sl <- isTRUE(g$ordered[2L]) && nrow(res) == length(x[[1L]])
          if(sl) rn1 <- ax[["row.names"]][1L]
          dimnames(res) <- list(if(sl && length(rn1) && is.character(rn1) && rn1 != "1")
            ax[["row.names"]] else NULL, ax[["names"]])
          return(res)
        } else {
          res <- aplyfun(x, copysplaplfun, g, FUN, ...)
          if(length(res[[1L]]) != length(x[[1L]]) || !isTRUE(g$ordered[2L]))
            ax[["row.names"]] <- .set_row_names(length(res[[1L]]))
        }
      }
    }
    return(condalcSA(res, ax, isDTl))
  }
  if(expand.wide) return(aplyfun(x, function(y) do.call(rbind, lapply(gsplit(y, g, use.g.names), FUN, ...))))
  return(aplyfun(x, function(y) lapply(gsplit(y, g, use.g.names), FUN, ...)))
}

BY.list <- function(x, ...) BY.data.frame(x, ...)

BY.matrix <- function(x, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                      expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                      return = c("same", "matrix", "data.frame", "list")) {
  if(!is.matrix(x)) stop("x needs to be a matrix")
  if(!(is.function(FUN) || is.character(FUN))) stop("FUN needs to be a function")
  aplyfun <- if(parallel) function(...) parallel::mclapply(..., mc.cores = mc.cores) else lapply
  return <- switch(return[1L], same = 3L, matrix = 2L, data.frame = 1L, list = 0L,
                   stop("Unknown return option!"))
  g <- GRP(g, return.groups = use.g.names, sort = sort, call = FALSE)
  if(return != 0L) {
    if(expand.wide) {
      if(return == 1L) { # Return data frame
        splitfun <- function(y) .Call(Cpp_mctl, do.call(rbind, lapply(gsplit(y, g), FUN, ...)), TRUE, 0L)
        res <- unlist(aplyfun(.Call(Cpp_mctl, x, TRUE, 0L), splitfun), recursive = FALSE, use.names = TRUE)
        ax <- list(names = names(res),
                   row.names = if(use.g.names && length(gn <- GRPnames(g))) gn else .set_row_names(length(res[[1L]])),
                   class = "data.frame")
      } else { # Return a matrix
        splitfun2 <- function(y) do.call(rbind, lapply(gsplit(y, g), FUN, ...))
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
    } else {
      if(use.g.names && length(gn <- GRPnames(g))) {
        res <- vector("list", ncol(x))
        res1 <- lapply(gsplit(`names<-`(x[, 1L], NULL), g), FUN, ...)
        names(res1) <- gn
        res[[1L]] <- unlist(res1, FALSE, TRUE)
        namres1 <- names(res[[1L]])
        if(length(namres1)) names(res[[1L]]) <- NULL
        if(length(res) > 1L) res[-1L] <- aplyfun(.Call(Cpp_mctl, x[, -1L, drop = FALSE], FALSE, 0L), splaplfun, g, FUN, ...)
        if(return > 1L) { # Return a matrix
          res <- do.call(cbind, res)
          dn <- list(namres1, dimnames(x)[[2L]])
          if(return == 2L) return(`dimnames<-`(res, dn))
          ax <- attributes(x)
          ax[["dim"]] <- dim(res)
          ax[["dimnames"]] <- dn
        } else { # Return a data frame
          ax <- list(names = dimnames(x)[[2L]],
                     row.names = if(length(namres1)) namres1 else .set_row_names(length(res[[1L]])),
                     class = "data.frame")
        }
      } else {
        res <- aplyfun(.Call(Cpp_mctl, x, TRUE, 0L), splaplfun, g, FUN, ...)
        if(return > 1L) { # Return a matrix
          res <- do.call(cbind, res)
          if(return == 2L) return(res)
          ax <- attributes(x)
          if(length(dimnames(x)[[1L]]) && !(isTRUE(g$ordered[2L]) && nrow(res) == nrow(x))) {
            ax[["dimnames"]][1L] <- list(NULL)
            ax[["dim"]] <- dim(res)
          }
        } else { # Return a data frame
          lr1 <- length(res[[1L]])
          ax <- list(names = names(res),
                     row.names = if(lr1 == nrow(x) && length(rn <- dimnames(x)[[1L]]) && isTRUE(g$ordered[2L])) rn else .set_row_names(lr1),
                     class = "data.frame")
        }
      }
    }
    return(setAttributes(res, ax))
  }
  if(expand.wide) return(aplyfun(.Call(Cpp_mctl, x, TRUE, 0L), function(y) do.call(rbind, lapply(gsplit(y, g, use.g.names), FUN, ...))))
  return(aplyfun(.Call(Cpp_mctl, x, TRUE, 0L), function(y) lapply(gsplit(y, g, use.g.names), FUN, ...)))
}

BY.grouped_df <- function(x, FUN, ..., keep.group_vars = TRUE, use.g.names = FALSE) {
  g <- GRP.grouped_df(x, call = FALSE)
  gn <- which(attr(x, "names") %in% g[[5L]])
  if(!length(gn)) {
   if(!isTRUE(g$ordered[2L])) return(BY.data.frame(fungroup(x), g, FUN, ..., use.g.names = use.g.names))
    res <- BY.data.frame(x, g, FUN, ..., use.g.names = use.g.names)
    if(!is.data.frame(res) || fnrow2(res) == fnrow2(x)) return(res) else return(fungroup(res))
  }
  res <- BY.data.frame(fcolsubset(x, -gn), g, FUN, ..., use.g.names = use.g.names)
  if(!is.data.frame(res)) return(res)
  nrr <- fnrow2(res)
  same_size <- nrr == fnrow2(x)
  if(!keep.group_vars) return(if(same_size && isTRUE(g$ordered[2L])) res else fungroup(res))
  if(!((same_size && isTRUE(g$ordered[2L])) || nrr == g[[1L]])) return(fungroup(res))
  if(same_size) {
    ar <- attributes(res)
    ar[["names"]] <- c(g[[5L]], ar[["names"]])
    return(condalcSA(c(.subset(x, gn), res), ar, any(ar$class == "data.table")))
  }
  ar <- attributes(fungroup(res))
  attributes(res) <- NULL
  ar[["names"]] <- c(g[[5L]], ar[["names"]])
  condalcSA(c(g[[4L]], res), ar, any(ar$class == "data.table"))
}
