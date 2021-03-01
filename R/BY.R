

BY <- function(x, ...) UseMethod("BY") # , x

BY.default <- function(x, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                       expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                       return = c("same", "vector", "list")) { # what about ... in those other internal calls ?
  if(!is.atomic(x)) stop("x needs to be an atomic vector") # redundant ?
  if(is.matrix(x) && !inherits(x, "matrix")) return(UseMethod("BY", unclass(x)))
  if(!(is.function(FUN) || is.character(FUN))) stop("FUN needs to be a function")
  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply
  simplify <- switch(return[1L], same = 1L, vector = 2L, list = 3L, stop("BY.default only supports same, vector and list output!"))
  if(!is.factor(g)) g <- if(is.GRP(g)) as.factor_GRP(g) else if(is.list(g))
                         as.factor_GRP(GRP.default(g, sort = sort, call = FALSE)) else
                           qF(g, sort = sort, na.exclude = FALSE)
    res <- aplyfun(split(x, g), FUN, ...)
    if(simplify < 3L) {
      if(expand.wide) {
        res <- do.call(rbind, res)
        if(!use.g.names) dimnames(res) <- list(NULL, dimnames(res)[[2]])
      } else {
        if(use.g.names) {
          res <- unlist(res, recursive = FALSE)
          if(simplify == 1L) { #  && typeof(res) == typeof(x)   # length(res) == length(x) && # manually control it...
            ax <- attributes(x)
            if(length(ax)) {
              ax[["names"]] <- names(res)
              setattributes(res, ax) # attributes(res) <- ax
            }
          }
        } else {
          if(simplify == 1L) return(duplAttributes(unlist(res, recursive = FALSE, use.names = FALSE), x)) #  && typeof(res) == typeof(x) # length(res) == length(x) &&   what if x has names attribute ?
          ll <- length(res)
          nr1 <- names(res[[1L]]) # good solution ?
          res <- unlist(res, recursive = FALSE, use.names = FALSE)
          if(length(res) != ll) { # attributes(res) <- attributes(x)
            if(length(nr1) && length(res) == length(nr1)*ll) # additional check..
            names(res) <- rep(nr1, ll)
          }
        }
      }
    } else if(!use.g.names) names(res) <- NULL
  res
}

BY.data.frame <- function(x, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                          expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                          return = c("same", "matrix", "data.frame", "list")) {
  if(!is.list(x)) stop("x needs to be a list")
  if(!(is.function(FUN) || is.character(FUN))) stop("FUN needs to be a function")
  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply
  return <- switch(return[1L], same = 1L, matrix = 3L, data.frame = 2L, list = 0L,
                   stop("Unknown return option!"))
  if(!is.factor(g)) g <- if(is.GRP(g)) as.factor_GRP(g) else if(is.list(g))
                    as.factor_GRP(GRP.default(g, sort = sort, call = FALSE)) else
                      qF(g, sort = sort, na.exclude = FALSE)
  if(return != 0L) {
    ax <- attributes(x)
    if(expand.wide) {
      if(return < 3L) { # Return a data.frame
        splitfun <- function(y) .Call(Cpp_mctl, do.call(rbind, lapply(split(y, g), FUN, ...)), TRUE, 0L)
        res <- unlist(aplyfun(x, splitfun), recursive = FALSE, use.names = TRUE)
        if(return == 1L) {
          ax[["names"]] <- names(res)
          ax[["row.names"]] <- if(use.g.names && !inherits(x, "data.table")) attr(g, "levels") else .set_row_names(length(res[[1L]])) # faster than nlevels ?
        } else
          ax <- list(names = names(res), row.names = if(use.g.names) attr(g, "levels") else .set_row_names(length(res[[1L]])), class = "data.frame")
      } else { # Return a matrix
        attributes(x) <- NULL
        splitfun <- function(y) do.call(rbind, lapply(split.default(y, g), FUN, ...))
        res <- do.call(cbind, aplyfun(x, splitfun))
        dr <- dim(res)
        dn <- dimnames(res) # works, but what if dn[[2L]] is NULL ?
        if(!use.g.names) dn[1L] <- list(NULL) # character(0)
        dn[[2L]] <- paste(rep(ax[["names"]], each = dr[2L]/length(x)), dn[[2L]], sep = ".")
        ax <- list(dim = dr, dimnames = dn)
      }
    } else { # No expand wide (classical result)
      matl <- return == 3L
      if(return == 2L) ax <- list(names = ax[["names"]], row.names = ax[["row.names"]], class = "data.frame")
      if(use.g.names && (matl || !inherits(x, "data.table"))) { # using names...
        attributes(x) <- NULL
        res <- vector("list", length(x))
        res[[1L]] <- unlist(lapply(split(`names<-`(x[[1L]], ax[["row.names"]]), g), FUN, ...), FALSE, TRUE)
        if(matl) dn <- list(names(res[[1L]]), ax[["names"]]) else ax[["row.names"]] <- names(res[[1L]])
        setattr(res[[1L]], "names", NULL) # faster than  names(res[[1]]) <- NULL
        if(!matl && typeof(res[[1L]]) == typeof(x[[1L]])) { # length(res[[1]]) == nrow(x) &&   always safe ?
          setattr(x[[1L]], "names", NULL)
          duplattributes(res[[1L]], x[[1L]])
          splitfun <- function(y) duplAttributes(unlist(lapply(split(y, g), FUN, ...), FALSE, FALSE), y)
        } else splitfun <- function(y) unlist(lapply(split(y, g), FUN, ...), FALSE, FALSE)
        res[-1L] <- aplyfun(x[-1L], splitfun)
        if(matl) {
          res <- do.call(cbind, res)
          ax <- list(dim = lengths(dn, FALSE), dimnames = dn)
        }
      } else { # Not using generated rownames.
        attributes(x) <- NULL
        if(matl) {
          splitfun <- function(y) unlist(lapply(split.default(y, g), FUN, ...), FALSE, FALSE)
          res <- do.call(cbind, aplyfun(x, splitfun))
          dimr <- dim(res)
          ax <- list(dim = dimr, dimnames = list(if(length(x[[1L]]) == dimr[1L] && ax[["row.names"]][1L] != "1") ax[["row.names"]] else NULL, ax[["names"]]))
        } else {
          splitfun <- function(y) cond_duplAttributes(unlist(lapply(split(y, g), FUN, ...), FALSE, FALSE), y)
          res <- aplyfun(x, splitfun)
          if(length(res[[1L]]) != length(x[[1L]])) ax[["row.names"]] <- .set_row_names(length(res[[1L]]))
        }
      }
    }
    return(setAttributes(res, ax))
  }
  if(expand.wide) return(aplyfun(x, function(y) do.call(rbind, lapply(split(y, g), FUN, ...))))
  if(use.g.names) return(aplyfun(x, function(y) lapply(split(y, g), FUN, ...))) else
  return(aplyfun(x, function(y) `names<-`(lapply(split(y, g), FUN, ...), NULL)))
}

BY.list <- function(x, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                    expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                    return = c("same", "matrix", "data.frame", "list"))
  BY.data.frame(x, g, FUN, ..., use.g.names, sort, expand.wide, parallel, mc.cores, return)

BY.matrix <- function(x, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                      expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                      return = c("same", "matrix", "data.frame", "list")) {
  if(!is.matrix(x)) stop("x needs to be a matrix")
  if(!(is.function(FUN) || is.character(FUN))) stop("FUN needs to be a function")
  aplyfun <- if(parallel) function(...) parallel::mclapply(..., mc.cores = mc.cores) else lapply
  return <- switch(return[1L], same = 3L, matrix = 2L, data.frame = 1L, list = 0L,
                   stop("Unknown return option!"))
  if(!is.factor(g)) g <- if(is.GRP(g)) as.factor_GRP(g) else if(is.list(g))
                    as.factor_GRP(GRP.default(g, sort = sort, call = FALSE)) else
                    qF(g, sort = sort, na.exclude = FALSE)
  if(return != 0L) {
    ax <- attributes(x)
    if(expand.wide) {
      if(return == 1L) {
        splitfun <- function(y) .Call(Cpp_mctl, do.call(rbind, lapply(split.default(y, g), FUN, ...)), TRUE, 0L)
        res <- unlist(aplyfun(.Call(Cpp_mctl, x, TRUE, 0L), splitfun), recursive = FALSE, use.names = TRUE)
        ax <- list(names = names(res), row.names = if(use.g.names) attr(g, "levels") else .set_row_names(length(res[[1L]])), class = "data.frame")
      } else {
        splitfun2 <- function(y) do.call(rbind, lapply(split.default(y, g), FUN, ...))
        res <- do.call(cbind, aplyfun(.Call(Cpp_mctl, x, FALSE, 0L), splitfun2))
        dr <- dim(res)
        if(return == 2L) ax <- list(dim = dr, dimnames = NULL) else ax[["dim"]] <- dr
        ax[["dimnames"]] <- list(if(use.g.names) attr(g, "levels") else NULL, # works, but what if dn[[2L]] is NULL ?
                                 paste(rep(dimnames(x)[[2L]], each = dr[2L]/ncol(x)), dimnames(res)[[2L]], sep = "."))
      }
    } else {
      splitfun3 <- function(y, un = FALSE) unlist(lapply(split.default(y, g), FUN, ...), FALSE, un)
      if(return == 2L) ax <- list(dim = dim(x), dimnames = dimnames(x))
      if(use.g.names) {
        res <- vector("list", ncol(x))
        res[[1L]] <- splitfun3(x[, 1L], un = TRUE) # rewrite all in C++ ? # Note: x[, 1L] still keeps row.names, which are then interacted with group names.
        if(return > 1L) ax[["dimnames"]] <- list(names(res[[1L]]), dimnames(x)[[2L]]) else
          ax <- list(names = dimnames(x)[[2L]], row.names = names(res[[1L]]), class = "data.frame")
        setattr(res[[1L]], "names", NULL)
        res[-1L] <- aplyfun(.Call(Cpp_mctl, x[, -1L, drop = FALSE], FALSE, 0L), splitfun3)
        if(return > 1L) {
          res <- do.call(cbind, res)
          ax[["dim"]] <- dim(res)
        }
      } else {
        res <- aplyfun(.Call(Cpp_mctl, x, TRUE, 0L), splitfun3) # internal ?
        lr1 <- length(res[[1L]])
        if(return > 1L) {
          if(lr1 != nrow(x)) {
            if(length(dimnames(x))) ax[["dimnames"]][1L] <- list(NULL) # character(0) doesn't work !, the if check guards for error if no dimnames
            ax[["dim"]][1L] <- lr1
          }
          res <- do.call(cbind, res)
        } else {
          ax <- list(names = names(res), row.names = if(lr1 == nrow(x) && length(rn <- dimnames(x)[[1L]]))
            rn else .set_row_names(lr1), class = "data.frame")
        }
      }
    }
    return(setAttributes(res, ax))
  }
  if(expand.wide) return(aplyfun(.Call(Cpp_mctl, x, TRUE, 0L), function(y) do.call(rbind, lapply(split.default(y, g), FUN, ...))))
  if(use.g.names) return(aplyfun(.Call(Cpp_mctl, x, TRUE, 0L), function(y) lapply(split.default(y, g), FUN, ...))) else
  return(aplyfun(.Call(Cpp_mctl, x, TRUE, 0L), function(y) `names<-`(lapply(split.default(y, g), FUN, ...), NULL)))
}

BY.grouped_df <- function(x, FUN, ..., use.g.names = FALSE, keep.group_vars = TRUE,
                          expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                          return = c("same", "matrix", "data.frame", "list")) {
  g <- GRP.grouped_df(x, call = FALSE)
  groups <- g[[4L]]
  gnam <- g[[5L]]
  g <- as.factor_GRP(g)
  gn <- which(attr(x, "names") %in% gnam) # correct !
  if(length(gn)) {
    if(!keep.group_vars) return(BY.data.frame(x[-gn], g, FUN, ..., # colsubset(x, -gn) dont use colsubset -> doesn't drop group attachment !, for the other cases can always use ungroup !
                           use.g.names = use.g.names, sort = TRUE, expand.wide = expand.wide,
                           parallel = parallel, mc.cores = mc.cores, return = return))
      res <- BY.data.frame(fcolsubset(x, -gn), g, FUN, ...,
                           use.g.names = use.g.names, sort = TRUE, expand.wide = expand.wide,
                           parallel = parallel, mc.cores = mc.cores, return = return)
      if(is.data.frame(res)) {
        nrr <- fnrow2(res)
        same_size <- nrr == fnrow2(x)
        if(same_size || all(nrr == lengths(groups, FALSE))) {
          if(same_size) {
            ax <- attributes(x)
            attributes(res) <- NULL # faster removing attributes? (yes, a tiny bit!)  also set attributes of groups NULL ? -> Nah
            ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
            return(setAttributes(c(.subset(x, gn), res), ax))
          }
          ax <- attributes(res)
          attributes(res) <- NULL
          ax[["groups"]] <- NULL
          ax[["class"]] <- fsetdiff(ax[["class"]], c("GRP_df", "grouped_df"))
          ax[["names"]] <- c(gnam, ax[["names"]])
          return(setAttributes(c(groups, res), ax))
        } else return(res)
      } else return(res)
  } else return(BY.data.frame(x, g, FUN, ..., use.g.names = use.g.names, sort = TRUE, expand.wide = expand.wide,
                              parallel = parallel, mc.cores = mc.cores, return = return))
}


# Notes / Experimental:

# What about split apply combining other data stricture i.e. factors, date and time ... -> try mode !!
# -> need fplit and unlist (original) to account for factors. Note that fsplit does not deal with date and time ... but unlist can't handle those either... but nobody aggregates dates anyway...

# fsplit <- function(y, f) {
#   if(is.null(attr(y, "class"))) .Call(Csplit, y, f)
#     # return(.Internal(split(y, f)))
#   lf <- levels(f)
#   z <- vector("list", length(lf))
#   names(z) <- lf
#   ind <- .Internal(split(seq_along(y), f))
#   for (k in lf) z[[k]] <- y[ind[[k]]]
#   z
# }

# fsplit <- split.default # slightly slower !! (Csplit not as fast as internal!! !!)
