

BY <- function(X, ...) UseMethod("BY") # , X

BY.default <- function(X, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                       expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                       return = c("same", "vector", "list")) { # what about ... in those other internal calls ?
  if(!is.atomic(X)) stop("X needs to be an atomic vector") # redundant ?
  if(is.matrix(X) && !inherits(X, "matrix")) return(UseMethod("BY", unclass(x)))
  if(!(is.function(FUN) || is.character(FUN))) stop("FUN needs to be a function")
  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply
  simplify <- switch(return[1L], same = 1L, vector = 2L, list = 3L, stop("BY.default only supports same, vector and list output!"))
  if(!is.factor(g)) g <- if(is.GRP(g)) as.factor_GRP(g) else if(is.list(g))
                         as.factor_GRP(GRP.default(g, sort = sort, call = FALSE)) else
                           qF(g, sort = sort, na.exclude = FALSE)
    res <- aplyfun(split.default(X, g), FUN, ...)
    if(simplify < 3L) {
      if(expand.wide) {
        res <- do.call(rbind, res)
        if(!use.g.names) dimnames(res) <- list(NULL, dimnames(res)[[2]])
      } else {
        if(use.g.names) {
          res <- unlist(res, recursive = FALSE)
          if(simplify == 1L) { #  && typeof(res) == typeof(X)   # length(res) == length(X) && # manually control it...
            ax <- attributes(X)
            if(length(ax)) {
              ax[["names"]] <- names(res)
              setattributes(res, ax) # attributes(res) <- ax
            }
          }
        } else {
          if(simplify == 1L) return(duplAttributes(unlist(res, recursive = FALSE, use.names = FALSE), X)) #  && typeof(res) == typeof(X) # length(res) == length(X) &&   what if X has names attribute ?
          ll <- length(res)
          nr1 <- names(res[[1L]]) # good solution ?
          res <- unlist(res, recursive = FALSE, use.names = FALSE)
          if(length(res) != ll) { # attributes(res) <- attributes(X)
            if(length(nr1) && length(res) == length(nr1)*ll) # additional check..
            names(res) <- rep(nr1, ll)
          }
        }
      }
    } else if(!use.g.names) names(res) <- NULL
  res
}

BY.data.frame <- function(X, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                          expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                          return = c("same", "matrix", "data.frame", "list")) {
  if(!is.list(X)) stop("X needs to be a list")
  if(!(is.function(FUN) || is.character(FUN))) stop("FUN needs to be a function")
  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply
  return <- switch(return[1L], same = 1L, matrix = 3L, data.frame = 2L, list = 0L,
                   stop("Unknown return option!"))
  if(!is.factor(g)) g <- if(is.GRP(g)) as.factor_GRP(g) else if(is.list(g))
                    as.factor_GRP(GRP.default(g, sort = sort, call = FALSE)) else
                      qF(g, sort = sort, na.exclude = FALSE)
  if(return != 0L) {
    ax <- attributes(X)
    if(expand.wide) {
      if(return < 3L) { # Return a data.frame
        splitfun <- function(x) .Call(Cpp_mctl, do.call(rbind, lapply(split.default(x, g), FUN, ...)), TRUE, 0L)
        res <- unlist(aplyfun(X, splitfun), recursive = FALSE, use.names = TRUE)
        if(return == 1L) {
          ax[["names"]] <- names(res)
          ax[["row.names"]] <- if(use.g.names && !inherits(X, "data.table")) attr(g, "levels") else .set_row_names(length(res[[1L]])) # faster than nlevels ?
        } else
          ax <- list(names = names(res), row.names = if(use.g.names) attr(g, "levels") else .set_row_names(length(res[[1L]])), class = "data.frame")
      } else { # Return a matrix
        attributes(X) <- NULL
        splitfun <- function(x) do.call(rbind, lapply(split.default(x, g), FUN, ...))
        res <- do.call(cbind, aplyfun(X, splitfun))
        dr <- dim(res)
        dn <- dimnames(res) # works, but what if dn[[2L]] is NULL ?
        if(!use.g.names) dn[1L] <- list(NULL) # character(0)
        dn[[2L]] <- paste(rep(ax[["names"]], each = dr[2L]/length(X)), dn[[2L]], sep = ".")
        ax <- list(dim = dr, dimnames = dn)
      }
    } else { # No expand wide (classical result)
      matl <- return == 3L
      if(return == 2L) ax <- list(names = ax[["names"]], row.names = ax[["row.names"]], class = "data.frame")
      if(use.g.names && (matl || !inherits(X, "data.table"))) { # using names...
        attributes(X) <- NULL
        res <- vector("list", length(X))
        res[[1L]] <- unlist(lapply(split.default(`names<-`(X[[1L]], ax[["row.names"]]), g), FUN, ...), FALSE, TRUE)
        if(matl) dn <- list(names(res[[1L]]), ax[["names"]]) else ax[["row.names"]] <- names(res[[1L]])
        setattr(res[[1L]], "names", NULL) # faster than  names(res[[1]]) <- NULL
        if(!matl && typeof(res[[1L]]) == typeof(X[[1L]])) { # length(res[[1]]) == nrow(X) &&   always safe ?
          setattr(X[[1L]], "names", NULL)
          duplattributes(res[[1L]], X[[1L]])
          splitfun <- function(x) duplAttributes(unlist(lapply(split.default(x, g), FUN, ...), FALSE, FALSE), x)
        } else splitfun <- function(x) unlist(lapply(split.default(x, g), FUN, ...), FALSE, FALSE)
        res[-1L] <- aplyfun(X[-1L], splitfun)
        if(matl) {
          res <- do.call(cbind, res)
          ax <- list(dim = lengths(dn, FALSE), dimnames = dn)
        }
      } else { # Not using generated rownames.
        attributes(X) <- NULL
        if(matl) {
          splitfun <- function(x) unlist(lapply(split.default(x, g), FUN, ...), FALSE, FALSE)
          res <- do.call(cbind, aplyfun(X, splitfun))
          dimr <- dim(res)
          ax <- list(dim = dimr, dimnames = list(if(length(X[[1L]]) == dimr[1L] && ax[["row.names"]][1L] != "1") ax[["row.names"]] else NULL, ax[["names"]]))
        } else {
          splitfun <- function(x) cond_duplAttributes(unlist(lapply(split.default(x, g), FUN, ...), FALSE, FALSE), x)
          res <- aplyfun(X, splitfun)
          if(length(res[[1L]]) != length(X[[1L]])) ax[["row.names"]] <- .set_row_names(length(res[[1L]]))
        }
      }
    }
    return(setAttributes(res, ax))
  }
  if(expand.wide) return(aplyfun(X, function(x) do.call(rbind, lapply(split.default(x, g), FUN, ...))))
  if(use.g.names) return(aplyfun(X, function(x) lapply(split.default(x, g), FUN, ...))) else
  return(aplyfun(X, function(x) `names<-`(lapply(split.default(x, g), FUN, ...), NULL)))
}

BY.list <- function(X, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                    expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                    return = c("same", "matrix", "data.frame", "list"))
  BY.data.frame(X, g, FUN, ..., use.g.names, sort, expand.wide, parallel, mc.cores, return)

BY.matrix <- function(X, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                      expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                      return = c("same", "matrix", "data.frame", "list")) {
  if(!is.matrix(X)) stop("X needs to be a matrix")
  if(!(is.function(FUN) || is.character(FUN))) stop("FUN needs to be a function")
  aplyfun <- if(parallel) function(...) parallel::mclapply(..., mc.cores = mc.cores) else lapply
  return <- switch(return[1L], same = 3L, matrix = 2L, data.frame = 1L, list = 0L,
                   stop("Unknown return option!"))
  if(!is.factor(g)) g <- if(is.GRP(g)) as.factor_GRP(g) else if(is.list(g))
                    as.factor_GRP(GRP.default(g, sort = sort, call = FALSE)) else
                    qF(g, sort = sort, na.exclude = FALSE)
  if(return != 0L) {
    ax <- attributes(X)
    if(expand.wide) {
      if(return == 1L) {
        splitfun <- function(x) .Call(Cpp_mctl, do.call(rbind, lapply(split.default(x, g), FUN, ...)), TRUE, 0L)
        res <- unlist(aplyfun(.Call(Cpp_mctl, X, TRUE, 0L), splitfun), recursive = FALSE, use.names = TRUE)
        ax <- list(names = names(res), row.names = if(use.g.names) attr(g, "levels") else .set_row_names(length(res[[1L]])), class = "data.frame")
      } else {
        splitfun2 <- function(x) do.call(rbind, lapply(split.default(x, g), FUN, ...))
        res <- do.call(cbind, aplyfun(.Call(Cpp_mctl, X, FALSE, 0L), splitfun2))
        dr <- dim(res)
        if(return == 2L) ax <- list(dim = dr, dimnames = NULL) else ax[["dim"]] <- dr
        ax[["dimnames"]] <- list(if(use.g.names) attr(g, "levels") else NULL, # works, but what if dn[[2L]] is NULL ?
                                 paste(rep(dimnames(X)[[2L]], each = dr[2L]/ncol(X)), dimnames(res)[[2L]], sep = "."))
      }
    } else {
      splitfun3 <- function(x, un = FALSE) unlist(lapply(split.default(x, g), FUN, ...), FALSE, un)
      if(return == 2L) ax <- list(dim = dim(X), dimnames = dimnames(X))
      if(use.g.names) {
        res <- vector("list", ncol(X))
        res[[1L]] <- splitfun3(X[, 1L], un = TRUE) # rewrite all in C++ ? # Note: X[, 1L] still keeps row.names, which are then interacted with group names.
        if(return > 1L) ax[["dimnames"]] <- list(names(res[[1L]]), dimnames(X)[[2L]]) else
          ax <- list(names = dimnames(X)[[2L]], row.names = names(res[[1L]]), class = "data.frame")
        setattr(res[[1L]], "names", NULL)
        res[-1L] <- aplyfun(.Call(Cpp_mctl, X[, -1L, drop = FALSE], FALSE, 0L), splitfun3)
        if(return > 1L) {
          res <- do.call(cbind, res)
          ax[["dim"]] <- dim(res)
        }
      } else {
        res <- aplyfun(.Call(Cpp_mctl, X, TRUE, 0L), splitfun3) # internal ?
        lr1 <- length(res[[1L]])
        if(return > 1L) {
          if(lr1 != nrow(X)) {
            if(length(dimnames(X))) ax[["dimnames"]][1L] <- list(NULL) # character(0) doesn't work !, the if check guards for error if no dimnames
            ax[["dim"]][1L] <- lr1
          }
          res <- do.call(cbind, res)
        } else {
          ax <- list(names = names(res), row.names = if(lr1 == nrow(X) && length(rn <- dimnames(X)[[1L]]))
            rn else .set_row_names(lr1), class = "data.frame")
        }
      }
    }
    return(setAttributes(res, ax))
  }
  if(expand.wide) return(aplyfun(.Call(Cpp_mctl, X, TRUE, 0L), function(x) do.call(rbind, lapply(split.default(x, g), FUN, ...))))
  if(use.g.names) return(aplyfun(.Call(Cpp_mctl, X, TRUE, 0L), function(x) lapply(split.default(x, g), FUN, ...))) else
  return(aplyfun(.Call(Cpp_mctl, X, TRUE, 0L), function(x) `names<-`(lapply(split.default(x, g), FUN, ...), NULL)))
}

BY.grouped_df <- function(X, FUN, ..., use.g.names = FALSE, keep.group_vars = TRUE,
                          expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                          return = c("same", "matrix", "data.frame", "list")) {
  g <- GRP.grouped_df(X, call = FALSE)
  groups <- g[[4L]]
  gnam <- g[[5L]]
  g <- as.factor_GRP(g)
  gn <- which(attr(X, "names") %in% gnam) # correct !
  if(length(gn)) {
    if(!keep.group_vars) return(BY.data.frame(X[-gn], g, FUN, ..., # colsubset(X, -gn) dont use colsubset -> doesn't drop group attachment !, for the other cases can always use ungroup !
                           use.g.names = use.g.names, sort = TRUE, expand.wide = expand.wide,
                           parallel = parallel, mc.cores = mc.cores, return = return))
      res <- BY.data.frame(fcolsubset(X, -gn), g, FUN, ...,
                           use.g.names = use.g.names, sort = TRUE, expand.wide = expand.wide,
                           parallel = parallel, mc.cores = mc.cores, return = return)
      if(is.data.frame(res)) {
        nrr <- fnrow2(res)
        same_size <- nrr == fnrow2(X)
        if(same_size || all(nrr == lengths(groups, FALSE))) {
          if(same_size) {
            ax <- attributes(X)
            attributes(res) <- NULL # faster removing attributes? (yes, a tiny bit!)  also set attributes of groups NULL ? -> Nah
            ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
            return(setAttributes(c(.subset(X, gn), res), ax))
          }
          ax <- attributes(res)
          attributes(res) <- NULL
          ax[["groups"]] <- NULL
          ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
          ax[["names"]] <- c(gnam, ax[["names"]])
          return(setAttributes(c(groups, res), ax))
        } else return(res)
      } else return(res)
  } else return(BY.data.frame(X, g, FUN, ..., use.g.names = use.g.names, sort = TRUE, expand.wide = expand.wide,
                              parallel = parallel, mc.cores = mc.cores, return = return))
}


# Notes / Experimental:

# What about split apply combining other data stricture i.e. factors, date and time ... -> try mode !!
# -> need fplit and unlist (original) to account for factors. Note that fsplit does not deal with date and time ... but unlist can't handle those either... but nobody aggregates dates anyway...

# fsplit <- function(x, f) {
#   if(is.null(attr(x, "class"))) .Call(Csplit, x, f)
#     # return(.Internal(split(x, f)))
#   lf <- levels(f)
#   y <- vector("list", length(lf))
#   names(y) <- lf
#   ind <- .Internal(split(seq_along(x), f))
#   for (k in lf) y[[k]] <- x[ind[[k]]]
#   y
# }

# fsplit <- split.default # slightly slower !! (Csplit not as fast as internal!! !!)
