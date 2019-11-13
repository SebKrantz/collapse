# todo:: remove attributes ??

BY.data.frame <- function(X, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
                          expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
                          return = c("same","matrix","data.frame","list")) {
  aplyfun <- if(parallel) function(...) parallel::mclapply(..., mc.cores = mc.cores) else lapply
  return <- switch(return[1L], same = 1L, matrix = 2L, data.frame = 1L, list = 0L,
                   stop("Unknown return option!"))
  if(!is.factor(g)) g <- if(is.GRP(g)) as.factor.GRP(g) else if(is.list(g))
                         as.factor.GRP(GRP(g, sort = sort)) else qF(g, ordered = sort)
  if(return != 0L) {
    ax <- attributes(X)
    if(expand.wide) {
      if(return == 1L) {
        splitfun <- function(x) mctl(do.call(rbind, lapply(fsplit(x, g), FUN, ...)), names = TRUE)
        res <- unlist(aplyfun(X, splitfun), recursive = FALSE, use.names = TRUE)
        ax[["row.names"]] <- if(use.g.names && !inherits(X, "data.table")) attr(g, "levels") else .set_row_names(length(res[[1L]])) # faster than nlevels ??
        ax[["names"]] <- names(res)
      } else {
        splitfun <- function(x) do.call(rbind, lapply(fsplit(x, g), FUN, ...))
        res <- aplyfun(X, splitfun)
        nam <- names(res)
        ax[["row.names"]] <- if(use.g.names && !inherits(X, "data.table")) attr(g, "levels") else .set_row_names(length(res[[1L]])) # faster than nlevels ??
        ax[["names"]] <- names(res)

      }
    } else {
      attributes(X) <- NULL
      if(use.g.names && !inherits(X, "data.table")) {
        res <- vector("list", length(X))
        res[[1L]] <- unlist(lapply(fsplit(X[[1L]], g), FUN, ...), FALSE, TRUE)
        ax[["row.names"]] <- names(res[[1L]])
        setattr_clp(res[[1L]], "names", NULL) # faster than  names(res[[1]]) <- NULL
        if(typeof(res[[1L]]) == typeof(X[[1L]])) { # length(res[[1]]) == nrow(X) &&   safe ??
          duplattributes(res[[1L]], X[[1L]])
          splitfun <- function(x) duplAttributes(unlist(lapply(fsplit(x, g), FUN, ...), FALSE, FALSE), x)
        } else splitfun <- function(x) unlist(lapply(fsplit(x, g), FUN, ...), FALSE, FALSE)
        res[-1L] <- aplyfun(X[-1L], splitfun)
      } else {
        splitfun <- function(x) cond_duplAttributes(unlist(lapply(fsplit(x, g), FUN, ...), FALSE, FALSE), x)
        res <- aplyfun(X, splitfun)
        lr1 <- length(res[[1L]])
        if(lr1 != nrow(X)) ax[["row.names"]] <- .set_row_names(lr1)
      }
    }
    return(setAttributes(res, ax))
  } else {
    if(expand.wide) return(aplyfun(X, function(x) do.call(rbind, lapply(fsplit(x, g), FUN, ...)))) else {
      if(use.g.names) return(aplyfun(X, function(x) lapply(fsplit(x, g), FUN, ...))) else
        return(aplyfun(X, function(x) `names<-`(lapply(fsplit(x, g), FUN, ...), NULL)))
    }
  }
}
