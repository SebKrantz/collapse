# library(Rcpp)
# sourceCpp("C:/Users/Sebastian Krantz/Documents/R/mrtl_type_dispatch.cpp")
# todo, Make setdapply(), or Xcols option??. Should be consistent with the rest. all set-functions should have copy option, or: if assigned to something make copy, else no!!
# definitely make add option for MARGIN = 1?? -> nah, rather setdapply, with Xcols and Add option
# -> Do just like B and W, both dapply and setdapply, and each has Xcols and add = 0,1,2 option.
# Note: in setdapply, 0 just means replace columns?? -> yes!!

# same as dapply 3 (compact), but takingdrop case before !! -> faster !! and also solving issue with row.names for matrices -> row and column names must be of same type !! as.matrix.data.frame converst row.names to character !!
dapply <- function(X, FUN, ..., MARGIN = 2, parallel = FALSE, # drop argument !!
                   mc.cores = 1L, return = c("same","matrix","data.frame"), drop = TRUE) {
  ax <- attributes(X)
  arl <- is.array(X)
  rowwl <- switch(MARGIN, `1` = TRUE, `2` = FALSE, stop("MARGIN only supports 2 - columns or 1 - rows"))
  retmatl <- switch(return[1L], same = arl, matrix = TRUE, data.frame = FALSE, stop("Unknown return option!"))
  aplyfun <- if(parallel) function(...) parallel::mclapply(..., mc.cores = mc.cores) else lapply
  if(arl) {
    dX <- dim(X)
    if(length(dX) > 2L) stop("dapply cannot handle higher-dimensional arrays")
    res <- if(rowwl) aplyfun(mrtl(X), FUN, ...) else aplyfun(mctl(X), FUN, ...)
    lx1 <- length(res[[1L]])
    if(lx1 == 1L && drop) return(setNames(unlist(res, use.names = FALSE), ax[["dimnames"]][[if(rowwl) 1L else 2L]]))
    if(!retmatl) {
      dn <- dimnames(X) # ax[["dimnames"]] # best ?? -> use res instead of reassigning X !!! -> no memory loss !!
      ax <- list(names = dn[[2L]], row.names = if(is.null(dn[[1L]])) .set_row_names(dX[1L]) else dn[[1L]],
                   class = "data.frame") # c( ... , ax[!(names(ax) %in% c("dim","dimnames","class"))]) # don't know why one would need this !!
    }
  } else {
    attributes(X) <- NULL
    dX <- c(length(X[[1L]]), length(X)) # much faster than dim(X) on a list !!
    res <- if(rowwl) aplyfun(mrtl(do.call(cbind, X)), FUN, ...) else aplyfun(X, FUN, ...) # do.call(cbind, X) is definitely faster than unlist(X, use.names = FALSE) and attaching dim attribute
    lx1 <- length(res[[1L]])
    if(lx1 == 1L && drop) return(setNames(unlist(res, use.names = FALSE), if(rowwl) charorNULL(ax[["row.names"]]) else ax[["names"]]))
    if(retmatl) ax <- list(dim = dX, dimnames = list(charorNULL(ax[["row.names"]]), ax[["names"]]))
                      # c(..., ax[!(names(ax) %in% c("names","row.names","class"))]) # don't know why one would need this !!
  }
  if(retmatl) {
    if(rowwl) {
      if(lx1 != dX[2L]) {
        ax[["dim"]][2L] <- lx1
        ax[["dimnames"]] <- list(ax[["dimnames"]][[1L]], if(!is.null(nx1 <- names(res[[1L]]))) nx1 else if(lx1 == 1L)
                    deparse(substitute(FUN)) else paste0(deparse(substitute(FUN)), seq_len(lx1)))
      }
      res <- matrix(unlist(res, use.names = FALSE), ncol = lx1, byrow = TRUE)
    } else {
      if(lx1 != dX[1L]) {
        ax[["dim"]][1L] <- lx1
        ax[["dimnames"]] <- list(if(!is.null(nx1 <- names(res[[1L]]))) nx1 else if(lx1 == 1L)
        deparse(substitute(FUN)) else paste0(deparse(substitute(FUN)), seq_len(lx1)), ax[["dimnames"]][[2L]])
      }
      res <- do.call(cbind, res)
    }
  } else {
    if(rowwl) {
      if(lx1 != dX[2L]) ax[["names"]] <- if(!is.null(nx1 <- names(res[[1L]]))) nx1 else if(lx1 == 1L)
        deparse(substitute(FUN)) else paste0(deparse(substitute(FUN)), seq_len(lx1))
      res <- mctl(matrix(unlist(res, use.names = FALSE), ncol = lx1, byrow = TRUE)) # definitely faster than do.call(rbind, X)
    } else if(lx1 != dX[1L])
      ax[["row.names"]] <- if(!is.null(nx1 <- names(res[[1L]]))) nx1 else .set_row_names(lx1) # could also make deparse(substitute(FUN)), but that is not so typical for data.frames !!
  }
  return(setAttributes(res, ax))
}
