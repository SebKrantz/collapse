# keep.w toggle w being kept even if passed externally ? -> Also not done with W, B , etc !! -> but they also don't keep by ..
collap <- function(X, by, FUN = fmean, catFUN = fmode, cols = NULL, w = NULL, wFUN = fsum, custom = NULL,
                   keep.by = TRUE, keep.w = TRUE, keep.col.order = TRUE, sort.row = TRUE,
                   parallel = FALSE, mc.cores = 1L,
                   return = c("wide","list","long","long_dupl"), give.names = "auto", ...) {

  return <- switch(return[1L], wide = 1L, list = 2L, long = 3L, long_dupl = 4L, stop("Unknown return output option"))
  widel <- return == 1L
  ncustoml <- is.null(custom)
  nwl <- is.null(w)
  if(!inherits(X, "data.frame")) X <- qDF(X)
  ax <- attributes(X)
  class(X) <- NULL
  nam <- names(X)
  # attributes(X) <- NULL
  # attr(X, "class") <- "data.frame" # class needed for method dispatch of fast functions, not for BY !

  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply

  # identifying by and cols
  vl <- TRUE
  bycalll <- is.call(by)
  if(bycalll) {
    if(length(by) == 3L) {
      v <- logical(length(X))
      v[ckmatch(all.vars(by[[2L]]), nam)] <- TRUE  # nam %in% all.vars(by[[2L]])
      namby <- all.vars(by[[3L]])
      numby <- ckmatch(namby, nam)
    } else {
      namby <- all.vars(by)
      numby <- ckmatch(namby, nam)
      if(ncustoml) {
        v <- if(is.null(cols)) !logical(length(X)) else cols2log(cols, X, nam)
        v[numby] <- FALSE
      }
    }
    by <- GRP.default(X, numby, sort = sort.row, return.groups = keep.by)
  } else if(is.atomic(by)) {
    namby <- l1orlst(as.character(substitute(by)))
    numby <- 1L
    if(ncustoml) if(is.null(cols)) vl <- FALSE else v <- cols2log(cols, X, nam)
    by <- GRP.default(`names<-`(list(by), namby), sort = sort.row, return.groups = keep.by)
  } else {
    if(ncustoml) if(is.null(cols)) vl <- FALSE else v <- cols2log(cols, X, nam)
    if(!is.GRP(by)) {
      numby <- seq_along(unclass(by))
      namby <- attr(by, "names") # faster if and only if by is a data.frame
      if(is.null(namby)) namby <- paste0("Group.", numby)
      by <- GRP.default(by, numby, sort = sort.row, return.groups = keep.by)
    } else {
      namby <- by[[5L]]
      if(is.null(namby)) namby <- paste0("Group.", seq_along(by[[4L]])) # necessary ?
      numby <- seq_along(namby)
    }
  }

  if(!nwl) {
    if(is.call(w)) {
      namw <- all.vars(w)
      numw <- ckmatch(namw, nam)
      if(vl && ncustoml) v[numw] <- FALSE
      w <- X[[numw]]
    } else if(keep.w) {
      numw <- 0L # length(X) + 1L
      namw <- l1orlst(as.character(substitute(w)))
    }
    if(keep.w) { # what about function name for give.names ?? What about give.names and custom ???
      if(is.function(wFUN)) {
        namwFUN <- l1orn(as.character(substitute(wFUN)), "wFUN")
      } else if(is.character(wFUN)) {
        namwFUN <- wFUN
        wFUN <- if(length(wFUN) > 1L) lapply(wFUN, match.fun, descend = FALSE) else match.fun(wFUN, descend = FALSE)
      } else if(is.list(wFUN)) {
        namwFUN <- names(wFUN)
        if(is.null(namwFUN)) namwFUN <- all.vars(substitute(wFUN))
      } else stop("wFUN needs to be a function, character vector of function names or list of functions!")
      if(!all(namwFUN %in% .FAST_STAT_FUN)) stop("wFUN needs to be fast statistical functions, see print(.FAST_STAT_FUN)")
      if(is.list(wFUN)) {
        namw <- paste(namwFUN, namw, sep = ".")
        by[[4L]] <- c(if(keep.by) by[[4L]], `names<-`(lapply(wFUN, function(f) f(w, by, ..., use.g.names = FALSE)), namw))
        if(keep.col.order) numby <- c(if(keep.by) numby, rep_len(numw, length(wFUN)))
      } else {
        if(isTRUE(give.names)) namw <- paste(namwFUN, namw, sep = ".")
        by[[4L]] <- c(if(keep.by) by[[4L]], `names<-`(list(wFUN(w, by, ..., use.g.names = FALSE)), namw))
        if(keep.col.order) numby <- c(if(keep.by) numby, numw)  # need to accommodate any option of keep.by, keep.w and keep.col.order
      }
      if(return == 2L) namby <- c(if(keep.by) namby, namw)
      keep.by <- TRUE
    }
  }

  if(ncustoml) {

    # Identifying data
    nu <- vapply(unattrib(X), is.numeric, TRUE)
    if(vl) {
      nnu <- which(!nu & v) # faster way ?
      nu <- which(nu & v)
    } else {
      nnu <- which(!nu)
      nu <- which(nu)
    }
    nul <- length(nu) > 0L
    nnul <- length(nnu) > 0L

    # Identifying FUN and catFUN:
    if(nul) if(is.function(FUN)) {
      namFUN <- l1orn(as.character(substitute(FUN)), "FUN")
    } else if(is.character(FUN)) {
      # FUN <- unlist(strsplit(FUN,",",fixed = TRUE), use.names = FALSE)
      namFUN <- FUN
      FUN <- if(length(FUN) > 1L) lapply(FUN, match.fun, descend = FALSE) else match.fun(FUN, descend = FALSE)
    } else if(is.list(FUN)) {
      namFUN <- names(FUN)
      if(is.null(namFUN)) namFUN <- all.vars(substitute(FUN))
    } else stop("FUN needs to be a function, character vector of function names or list of functions!")

    if(nnul) if(is.function(catFUN)) {
      namcatFUN <- l1orn(as.character(substitute(catFUN)), "catFUN")
    } else if(is.character(catFUN)) {
      # catFUN <- unlist(strsplit(catFUN,",",fixed = TRUE), use.names = FALSE)
      namcatFUN <- catFUN
      catFUN <- if(length(catFUN) > 1L) lapply(catFUN, match.fun, descend = FALSE) else match.fun(catFUN, descend = FALSE)
    } else if(is.list(catFUN)) {
      namcatFUN <- names(catFUN)
      if(is.null(namcatFUN)) namcatFUN <- all.vars(substitute(catFUN))
    } else stop("FUN needs to be a function, character vector of function names or list of functions!")

    if(give.names == "auto") give.names <- !widel || length(FUN) > 1L || length(catFUN) > 1L

    # Aggregator function # drop level of nesting i.e. make rest length(by)+length(FUN)+length(catFUN)  ?
    agg <- function(xnu, xnnu, ...) { # by, FUN, namFUN, catFUN, namcatFUN, drop.by
      lr <- nul + nnul + keep.by
      res <- vector("list", lr)
      if(keep.by) {
        res[[1L]] <- list(by[[4L]]) # could add later using "c" ?
        ind <- 2L
      } else ind <- 1L
      if(nul) {
        fFUN <- namFUN %in% .FAST_STAT_FUN
        if(is.list(FUN))
          res[[ind]] <- condsetn(aplyfun(seq_along(namFUN), function(i)
            if(fFUN[i]) FUN[[i]](xnu, by, ..., use.g.names = FALSE) else
              BY.data.frame(xnu, by, FUN[[i]], ..., use.g.names = FALSE)), namFUN, give.names) else
                res[[ind]] <- if(fFUN) condsetn(list(FUN(xnu, by, ..., use.g.names = FALSE)), namFUN, give.names) else # give.names || !widel
                  condsetn(list(BY.data.frame(xnu, by, FUN, ..., use.g.names = FALSE, parallel = parallel, mc.cores = mc.cores)), namFUN, give.names) # give.names || !widel
      }
      if(nnul) {
        fcatFUN <- namcatFUN %in% .FAST_STAT_FUN
        if(is.list(catFUN))
          res[[lr]] <- condsetn(aplyfun(seq_along(namcatFUN), function(i)
            if(fcatFUN[i]) catFUN[[i]](xnnu, by, ..., use.g.names = FALSE) else
              BY.data.frame(xnnu, by, catFUN[[i]], ..., use.g.names = FALSE)), namcatFUN, give.names) else
                res[[lr]] <- if(fcatFUN) condsetn(list(catFUN(xnnu, by, ..., use.g.names = FALSE)), namcatFUN, give.names) else # give.names || !widel
                  condsetn(list(BY.data.frame(xnnu, by, catFUN, ..., use.g.names = FALSE, parallel = parallel, mc.cores = mc.cores)), namcatFUN, give.names) # give.names || !widel
      }
      return(res)
    } # fastest using res list ?? or better combine at the end ??

    if(nwl) {
      res <- agg(if(nul) `oldClass<-`(X[nu], "data.frame") else NULL, if(nnul) `oldClass<-`(X[nnu], "data.frame") else NULL, ...)
    } else {
      res <- agg(if(nul) `oldClass<-`(X[nu], "data.frame") else NULL, if(nnul) `oldClass<-`(X[nnu], "data.frame") else NULL, w = w, ...)
    }


    if(keep.col.order && widel) o <- forder.int(c(if(!keep.by) NULL else if(!bycalll) rep(0L,length(numby)) else numby,
                                                  if(nul) rep(nu,length(FUN)) else NULL,
                                                  if(nnul) rep(nnu,length(catFUN)) else NULL))

  } else { # custom aggregation:

    if(give.names == "auto") give.names <- TRUE
    namFUN <- names(custom)
    if(!is.list(custom) || is.null(namFUN)) stop("custom needs to be a named list, see ?collap")
    fFUN <- namFUN %in% .FAST_STAT_FUN
    if(!keep.by) {
      res <- vector("list", 1L)
      ind <- 1L
    } else {
      res <- vector("list", 2L)
      res[[1L]] <- list(by[[4L]]) # could add later using "c" ?
      ind <- 2L
    }
    custom <- lapply(custom, cols2int, X, nam) # could integrate below, but then reorder doesn't work !
    #lx <- length(X)
    # custom <- lapply(custom, function(x) if(is.numeric(x) && max(abs(x)) <= lx)
    #                          x else if(is.character(x)) ckmatch(x, nam) else
    #                          stop("custom list content must be variable names or suitable column indices"))

    if(nwl) {
      res[[ind]] <- condsetn(aplyfun(seq_along(namFUN), function(i)
        if(fFUN[i]) match.fun(namFUN[i])(`oldClass<-`(X[custom[[i]]], "data.frame"), by, ..., use.g.names = FALSE) else
          BY.data.frame(X[custom[[i]]], by, namFUN[i], ..., use.g.names = FALSE)), namFUN, give.names)
    } else {
      if(!all(fFUN)) warning("collap can only perform weighted aggregations with the fast statistical functions (see .FAST_STAT_FUN): Ignoring weights argument to other functions")
      res[[ind]] <- condsetn(aplyfun(seq_along(namFUN), function(i)
        if(fFUN[i]) match.fun(namFUN[i])(`oldClass<-`(X[custom[[i]]], "data.frame"), by, w = w, ..., use.g.names = FALSE) else
          BY.data.frame(X[custom[[i]]], by, namFUN[i], ..., use.g.names = FALSE)), namFUN, give.names)
    }

    if(keep.col.order && return != 2L) { # && widel
      o <- unlist(custom, use.names = FALSE)
      o <- forder.int(c(if(!keep.by) NULL else if(!bycalll) rep(0L,length(numby)) else numby, if(widel) o else unique.default(o)))
    }
  }

  if(widel) res <- unlist(unlist(res, FALSE), FALSE) else {
    if(length(FUN) > 1L || length(catFUN) > 1L || length(custom) > 1L) {
      res <- unlist(res, FALSE)
      if(return == 2L) {
        ax[["row.names"]] <- .set_row_names(by[[1L]])
        if(!keep.by) return(lapply(res, function(e) {
          ax[["names"]] <- names(e)
          return(setAttributes(e, ax)) })) else
            return(lapply(res[-1L], function(e) {
              ax[["names"]] <- c(namby, names(e))
              setAttributes(c(res[[1L]], e), ax) }))
      } else {
        if(return != 4L) {
          res <- if(!keep.by) .Call(C_rbindlist, res, TRUE, TRUE, "Function") else # data.table:::Crbindlist
            .Call(C_rbindlist, lapply(res[-1L], function(e) c(res[[1L]], e)), TRUE, TRUE, "Function")
        } else {
          if(!ncustoml || !(nul && nnul)) stop("long_dupl is only meaningful for aggregations with both numeric and categorical data, and multiple functions used for only one of the two data types!")
          mFUN <- length(FUN) > 1L
          nid <- if(mFUN) length(res) else 2L-!keep.by
          if(!keep.by) {
            res <- if(mFUN) lapply(res[-nid], function(e) c(e, res[[nid]])) else
              lapply(res[-nid], function(e) c(res[[nid]], e))
          } else res <- if(mFUN) lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], e, res[[nid]])) else
            lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], res[[nid]], e))
          res <- .Call(C_rbindlist, res, FALSE, FALSE, "Function")
        }
        if(keep.col.order)  o <- if(ncustoml) forder.int(c(0L, if(!keep.by) NULL else if(!bycalll) rep(0L,length(numby)) else numby, nu, nnu)) else c(1L, o + 1L)
      }
    } else message("return options other than 'wide' are only meaningful if multiple functions are used!")
  }

  if(keep.col.order) .Call(C_setcolorder, res, o) # data.table:::Csetcolorder
  ax[["names"]] <- names(res)
  ax[["row.names"]] <- .set_row_names(length(res[[1L]]))
  return(setAttributes(res, ax))
}


# collapv: allows vector input to by and w
collapv <- function(X, by, FUN = fmean, catFUN = fmode, cols = NULL, w = NULL, wFUN = fsum, custom = NULL,
                    keep.by = TRUE, keep.w = TRUE, keep.col.order = TRUE, sort.row = TRUE,
                    parallel = FALSE, mc.cores = 1L,
                    return = c("wide","list","long","long_dupl"), give.names = "auto", ...) {

  return <- switch(return[1L], wide = 1L, list = 2L, long = 3L, long_dupl = 4L, stop("Unknown return output option"))
  widel <- return == 1L
  ncustoml <- is.null(custom)
  nwl <- is.null(w)
  if(!inherits(X, "data.frame")) X <- qDF(X)
  ax <- attributes(X)
  class(X) <- NULL
  nam <- names(X)

  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply

  # identifying by
  numby <- cols2int(by, X, nam)
  namby <- nam[numby]
  by <- GRP.default(X, numby, sort = sort.row, return.groups = keep.by)
  if(ncustoml) {
    v <- if(is.null(cols)) !logical(length(X)) else cols2log(cols, X, nam)
    v[numby] <- FALSE
  }

  if(!nwl) {
    if(length(w) == 1L) {
      numw <- cols2int(w, X, nam)
      namw <- nam[numw]
      if(ncustoml) v[numw] <- FALSE
      w <- X[[numw]]
    } else if(keep.w) {
      numw <- 0L
      namw <- l1orlst(as.character(substitute(w)))
    }
    if(keep.w) {
      if(is.function(wFUN)) {
        namwFUN <- l1orn(as.character(substitute(wFUN)), "wFUN")
      } else if(is.character(wFUN)) {
        namwFUN <- wFUN
        wFUN <- if(length(wFUN) > 1L) lapply(wFUN, match.fun, descend = FALSE) else match.fun(wFUN, descend = FALSE)
      } else if(is.list(wFUN)) {
        namwFUN <- names(wFUN)
        if(is.null(namwFUN)) namwFUN <- all.vars(substitute(wFUN))
      } else stop("wFUN needs to be a function, character vector of function names or list of functions!")
      if(!all(namwFUN %in% .FAST_STAT_FUN)) stop("wFUN needs to be fast statistical functions, see print(.FAST_STAT_FUN)")
      if(is.list(wFUN)) {
        namw <- paste(namwFUN, namw, sep = ".")
        by[[4L]] <- c(if(keep.by) by[[4L]], `names<-`(lapply(wFUN, function(f) f(w, by, ..., use.g.names = FALSE)), namw))
        if(keep.col.order) numby <- c(if(keep.by) numby, rep_len(numw, length(wFUN)))
      } else {
        if(isTRUE(give.names)) namw <- paste(namwFUN, namw, sep = ".")
        by[[4L]] <- c(if(keep.by) by[[4L]], `names<-`(list(wFUN(w, by, ..., use.g.names = FALSE)), namw))
        if(keep.col.order) numby <- c(if(keep.by) numby, numw)
      }
      if(return == 2L) namby <- c(if(keep.by) namby, namw)
      keep.by <- TRUE
    }
  }

  if(ncustoml) {

    # Identifying data
    nu <- vapply(unattrib(X), is.numeric, TRUE)
    nnu <- which(!nu & v) # faster way ?
    nu <- which(nu & v)
    nul <- length(nu) > 0L
    nnul <- length(nnu) > 0L

    # Identifying FUN and catFUN:
    if(nul) if(is.function(FUN)) {
      namFUN <- l1orn(as.character(substitute(FUN)), "FUN")
    } else if(is.character(FUN)) {
      namFUN <- FUN
      FUN <- if(length(FUN) > 1L) lapply(FUN, match.fun, descend = FALSE) else
        match.fun(FUN, descend = FALSE)
    } else if(is.list(FUN)) {
      namFUN <- names(FUN)
      if(is.null(namFUN)) namFUN <- all.vars(substitute(FUN))
    } else stop("FUN needs to be a function, character vector of function names or list of functions!")

    if(nnul) if(is.function(catFUN)) {
      namcatFUN <- l1orn(as.character(substitute(catFUN)), "catFUN")
    } else if(is.character(catFUN)) {
      namcatFUN <- catFUN
      catFUN <- if(length(catFUN) > 1L) lapply(catFUN, match.fun, descend = FALSE) else
        match.fun(catFUN, descend = FALSE)
    } else if(is.list(catFUN)) {
      namcatFUN <- names(catFUN)
      if(is.null(namcatFUN)) namcatFUN <- all.vars(substitute(catFUN))
    } else stop("FUN needs to be a function, character vector of function names or list of functions!")

    if(give.names == "auto") give.names <- !widel || length(FUN) > 1L || length(catFUN) > 1L

    # Aggregator function
    agg <- function(xnu, xnnu, ...) {
      lr <- nul + nnul + keep.by
      res <- vector("list", lr)
      if(keep.by) {
        res[[1L]] <- list(by[[4L]])
        ind <- 2L
      } else ind <- 1L
      if(nul) {
        fFUN <- namFUN %in% .FAST_STAT_FUN
        if(is.list(FUN))
          res[[ind]] <- condsetn(aplyfun(seq_along(namFUN), function(i)
            if(fFUN[i]) FUN[[i]](xnu, by, ..., use.g.names = FALSE) else
              BY.data.frame(xnu, by, FUN[[i]], ..., use.g.names = FALSE)), namFUN, give.names) else
                res[[ind]] <- if(fFUN) condsetn(list(FUN(xnu, by, ..., use.g.names = FALSE)), namFUN, give.names) else # give.names || !widel
                  condsetn(list(BY.data.frame(xnu, by, FUN, ..., use.g.names = FALSE, parallel = parallel, mc.cores = mc.cores)), namFUN, give.names) # give.names || !widel
      }
      if(nnul) {
        fcatFUN <- namcatFUN %in% .FAST_STAT_FUN
        if(is.list(catFUN))
          res[[lr]] <- condsetn(aplyfun(seq_along(namcatFUN), function(i)
            if(fcatFUN[i]) catFUN[[i]](xnnu, by, ..., use.g.names = FALSE) else
              BY.data.frame(xnnu, by, catFUN[[i]], ..., use.g.names = FALSE)), namcatFUN, give.names) else
                res[[lr]] <- if(fcatFUN) condsetn(list(catFUN(xnnu, by, ..., use.g.names = FALSE)), namcatFUN, give.names) else # give.names || !widel
                  condsetn(list(BY.data.frame(xnnu, by, catFUN, ..., use.g.names = FALSE, parallel = parallel, mc.cores = mc.cores)), namcatFUN, give.names) # give.names || !widel
      }
      return(res)
    }

    if(nwl) {
      res <- agg(if(nul) `oldClass<-`(X[nu], "data.frame") else NULL, if(nnul) `oldClass<-`(X[nnu], "data.frame") else NULL, ...)
    } else {
      res <- agg(if(nul) `oldClass<-`(X[nu], "data.frame") else NULL, if(nnul) `oldClass<-`(X[nnu], "data.frame") else NULL, w = w, ...)
    }

    if(keep.col.order && widel) o <- forder.int(c(if(!keep.by) NULL else numby,
                                                  if(nul) rep(nu,length(FUN)) else NULL,
                                                  if(nnul) rep(nnu,length(catFUN)) else NULL))

  } else { # custom aggregation:

    if(give.names == "auto") give.names <- TRUE
    namFUN <- names(custom)
    if(!is.list(custom) || is.null(namFUN)) stop("custom needs to be a named list, see ?collap")
    fFUN <- namFUN %in% .FAST_STAT_FUN
    if(!keep.by) {
      res <- vector("list", 1L)
      ind <- 1L
    } else {
      res <- vector("list", 2L)
      res[[1L]] <- list(by[[4L]])
      ind <- 2L
    }

    custom <- lapply(custom, cols2int, X, nam)

    if(nwl) {
      res[[ind]] <- condsetn(aplyfun(seq_along(namFUN), function(i)
        if(fFUN[i]) match.fun(namFUN[i])(`oldClass<-`(X[custom[[i]]], "data.frame"), by, ..., use.g.names = FALSE) else
          BY.data.frame(X[custom[[i]]], by, namFUN[i], ..., use.g.names = FALSE)), namFUN, give.names)
    } else {
      if(!all(fFUN)) warning("collapv can only perform weighted aggregations with the fast statistical functions (see .FAST_STAT_FUN): Ignoring weights argument to other functions")
      res[[ind]] <- condsetn(aplyfun(seq_along(namFUN), function(i)
        if(fFUN[i]) match.fun(namFUN[i])(`oldClass<-`(X[custom[[i]]], "data.frame"), by, w = w, ..., use.g.names = FALSE) else
          BY.data.frame(X[custom[[i]]], by, namFUN[i], ..., use.g.names = FALSE)), namFUN, give.names)
    }

    if(keep.col.order && return != 2L) {
      o <- unlist(custom, use.names = FALSE)
      o <- forder.int(c(if(!keep.by) NULL else numby, if(widel) o else unique.default(o)))
    }
  }

  if(widel) res <- unlist(unlist(res, FALSE), FALSE) else {
    if(length(FUN) > 1L || length(catFUN) > 1L || length(custom) > 1L) {
      res <- unlist(res, FALSE)
      if(return == 2L) {
        ax[["row.names"]] <- .set_row_names(by[[1L]])
        if(!keep.by) return(lapply(res, function(e) {
          ax[["names"]] <- names(e)
          return(setAttributes(e, ax)) })) else
            return(lapply(res[-1L], function(e) {
              ax[["names"]] <- c(namby, names(e))
              setAttributes(c(res[[1L]], e), ax) }))
      } else {
        if(return != 4L) {
          res <- if(!keep.by) .Call(C_rbindlist, res, TRUE, TRUE, "Function") else # data.table:::Crbindlist
            .Call(C_rbindlist, lapply(res[-1L], function(e) c(res[[1L]], e)), TRUE, TRUE, "Function")
        } else {
          if(!ncustoml || !(nul && nnul)) stop("long_dupl is only meaningful for aggregations with both numeric and categorical data, and multiple functions used for only one of the two data types!")
          mFUN <- length(FUN) > 1L
          nid <- if(mFUN) length(res) else 2L-!keep.by
          if(!keep.by) {
            res <- if(mFUN) lapply(res[-nid], function(e) c(e, res[[nid]])) else
              lapply(res[-nid], function(e) c(res[[nid]], e))
          } else res <- if(mFUN) lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], e, res[[nid]])) else
            lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], res[[nid]], e))
          res <- .Call(C_rbindlist, res, FALSE, FALSE, "Function")
        }
        if(keep.col.order)  o <- if(ncustoml) forder.int(c(0L, if(!keep.by) NULL else numby, nu, nnu)) else c(1L, o + 1L)
      }
    } else message("return options other than 'wide' are only meaningful if multiple functions are used!")
  }

  if(keep.col.order) .Call(C_setcolorder, res, o) # data.table:::Csetcolorder
  ax[["names"]] <- names(res)
  ax[["row.names"]] <- .set_row_names(length(res[[1L]]))
  return(setAttributes(res, ax))
}


# For dplyr integration: takes grouped_df as input
collapg <- function(X, FUN = fmean, catFUN = fmode, cols = NULL, w = NULL, wFUN = fsum, custom = NULL,
                    keep.group_vars = TRUE, keep.w = TRUE, keep.col.order = TRUE, sort.row = TRUE,
                    parallel = FALSE, mc.cores = 1L,
                    return = c("wide","list","long","long_dupl"), give.names = "auto", ...) {
  by <- GRP.grouped_df(X)
  if(is.null(custom)) ngn <- attr(X, "names") %!in% by[[5L]] # Note: this always leaves grouping columns on the left still !
  clx <- class(X)
  attr(X, "groups") <- NULL
  oldClass(X) <- clx[clx != "grouped_df"]
  if(length(wsym <- as.character(substitute(w))) == 1L) { # Non-standard evaluation of w argument
    if(any(windl <- wsym == attr(X, "names"))) {
      assign(wsym, .subset2(X, wsym))
      if(is.null(custom)) X <- fcolsubset(X, ngn & !windl) # else X <- X # Needed ?? -> nope !!
      return(eval(substitute(collap(X, by, FUN, catFUN, cols, w, wFUN, custom,
                                    keep.group_vars, keep.w, keep.col.order, sort.row, parallel,
                                    mc.cores, return, give.names, ...)), list(w = wsym)))
    }
  }
  if(is.null(custom)) X <- fcolsubset(X, ngn) # else X <- X # because of non-standard eval.. X is "."
  return(eval(substitute(collap(X, by, FUN, catFUN, cols, w, wFUN, custom,
              keep.group_vars, keep.w, keep.col.order, sort.row, parallel,
              mc.cores, return, give.names, ...))))
}
