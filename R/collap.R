# Need generic version for column-parallel apply and aggregating weights..
fsum_uw <- function(x, g, w, ...) fsum(x, g, ...)
fprod_uw <- function(x, g, w, ...) fprod(x, g, ...)
fmean_uw <- function(x, g, w, ...) fmean(x, g, ...)
fmedian_uw <- function(x, g, w, ...) fmedian(x, g, ...)
fvar_uw <- function(x, g, w, ...) fvar(x, g, ...)
fsd_uw <- function(x, g, w, ...) fsd(x, g, ...)
fmode_uw <- function(x, g, w, ...) fmode(x, g, ...)
fnth_uw <- function(x, n, g, w, ...) fnth(x, n, g, ...)

fmin_uw <- function(x, g, w, ...) fmin(x, g, ...)
fmax_uw <- function(x, g, w, ...) fmax(x, g, ...)
ffirst_uw <- function(x, g, w, ...) ffirst(x, g, ...)
flast_uw <- function(x, g, w, ...) flast(x, g, ...)
fnobs_uw <- function(x, g, w, ...) fnobs(x, g, ...)
fndistinct_uw <- function(x, g, w, ...) fndistinct(x, g, ...)
fNobs_uw <- function(x, g, w, ...) fnobs(x, g, ...)
fNdistinct_uw <- function(x, g, w, ...) fndistinct(x, g, ...)


mymatchfun <- function(FUN) {
  if(is.function(FUN)) return(FUN)
  switch(tochar(FUN),
         # cat(paste0(FSF, " = ", FSF, ",\n"))
         fmean = fmean,
         fmedian = fmedian,
         fmode = fmode,
         fsum = fsum,
         fprod = fprod,
         fsd = fsd,
         fvar = fvar,
         fmin = fmin,
         fmax = fmax,
         fnth = fnth,
         ffirst = ffirst,
         flast = flast,
         fnobs = fnobs,
         fndistinct = fndistinct,
         fNobs = fnobs,
         fNdistinct = fndistinct,
         # cat(paste0(paste0(FSF, "_uw"), " = ", paste0(FSF, "_uw"), ",\n"))
         fmean_uw = fmean_uw,
         fmedian_uw = fmedian_uw,
         fmode_uw = fmode_uw,
         fsum_uw = fsum_uw,
         fprod_uw = fprod_uw,
         fsd_uw = fsd_uw,
         fvar_uw = fvar_uw,
         fmin_uw = fmin_uw,
         fmax_uw = fmax_uw,
         fnth_uw = fnth_uw,
         ffirst_uw = ffirst_uw,
         flast_uw = flast_uw,
         fnobs_uw = fnobs_uw,
         fndistinct_uw = fndistinct_uw,
         fNobs_uw = fnobs_uw,
         fNdistinct_uw = fndistinct_uw,
         match.fun(FUN)) # get(FUN, mode = "function", envir = parent.frame(2)) -> no error message
}

# Column-level parallel implementation
applyfuns_internal <- function(data, by, FUN, fFUN, parallel, cores, ...) {
  oldClass(data) <- "data.frame" # Needed for correct method dispatch for fast functions...
  if(length(FUN) > 1L) {
    if(parallel) return(lapply(seq_along(FUN), function(i)
            if(fFUN[i]) mclapply(data, FUN[[i]], g = by, ..., use.g.names = FALSE, mc.cores = cores) else
              BY.data.frame(data, by, FUN[[i]], ..., use.g.names = FALSE, reorder = FALSE, return = "data.frame", parallel = parallel, mc.cores = cores))) # mclapply(data, copysplaplfun, by, FUN[[i]], ..., mc.cores = cores)

      return(lapply(seq_along(FUN), function(i)
              if(fFUN[i]) FUN[[i]](data, g = by, ..., use.g.names = FALSE) else
                   BY.data.frame(data, by, FUN[[i]], ..., use.g.names = FALSE, reorder = FALSE, return = "data.frame"))) # lapply(data, copysplaplfun, by, FUN[[i]], ...)
  }
  if(is.list(FUN)) FUN <- FUN[[1L]]

  if(parallel) if(fFUN) return(list(mclapply(data, FUN, g = by, ..., use.g.names = FALSE, mc.cores = cores))) else
      return(list(BY.data.frame(data, by, FUN, ..., use.g.names = FALSE, reorder = FALSE, return = "data.frame", parallel = parallel, mc.cores = cores))) # return(list(mclapply(data, copysplaplfun, by, FUN, ..., mc.cores = cores)))

  if(fFUN) return(list(FUN(data, g = by, ..., use.g.names = FALSE)))
  return(list(BY.data.frame(data, by, FUN, ..., use.g.names = FALSE, reorder = FALSE, return = "data.frame"))) # return(list(lapply(data, copysplaplfun, by, FUN, ...)))
}

rbindlist_factor <- function(l, idcol = "Function") {
  nam <- names(l)
  names(l) <- NULL
  res <- .Call(C_rbindlist, l, TRUE, TRUE, idcol)
  attr(res[[1L]], "levels") <- if (length(nam)) nam else as.character(seq_along(l))
  oldClass(res[[1L]]) <- "factor"
  res
}


# NOTE: CUSTOM SEPARATOR doesn't work because of unlist() !

# keep.w toggle w being kept even if passed externally ? -> Also not done with W, B , etc !! -> but they also don't keep by ..
collap <- function(X, by, FUN = fmean, catFUN = fmode, cols = NULL, w = NULL, wFUN = fsum, custom = NULL,
                   ...,
                   keep.by = TRUE, keep.w = TRUE, keep.col.order = TRUE, sort = .op[["sort"]], decreasing = FALSE,
                   na.last = TRUE, return.order = sort, method = "auto", parallel = FALSE, mc.cores = 2L,
                   return = c("wide","list","long","long_dupl"), give.names = "auto") {

  return <- switch(return[1L], wide = 1L, list = 2L, long = 3L, long_dupl = 4L, stop("Unknown return output option"))
  widel <- return == 1L
  ncustoml <- is.null(custom)
  autorn <- is.character(give.names) && give.names == "auto"
  nwl <- is.null(w)
  if(inherits(X, "data.frame")) DTl <- inherits(X, "data.table") else {
    X <- qDF(X)
    DTl <- FALSE
  }
  ax <- attributes(X)
  oldClass(X) <- NULL
  if(.Call(C_fnrow, X) == 0L) stop("data passed to collap() has 0 rows.") #160, 0 rows can cause segfault...
  nam <- names(X)
  # attributes(X) <- NULL
  # attr(X, "class") <- "data.frame" # class needed for method dispatch of fast functions, not for BY !

  # cl <- if(parallel) makeCluster(mc.cores) else NULL
  # aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply

  # identifying by and cols
  vl <- TRUE
  if(is.call(by)) {
    if(length(by) == 3L) {
      v <- ckmatch(all.vars(by[[2L]]), nam)
      numby <- ckmatch(all.vars(by[[3L]]), nam)
    } else {
      numby <- ckmatch(all.vars(by), nam)
      if(ncustoml) v <- if(is.null(cols)) seq_along(X)[-numby] else cols2int(cols, X, nam)
    }
    by <- GRP.default(X, numby, sort, decreasing, na.last, keep.by, return.order, method, call = FALSE)
  } else if(is.atomic(by)) {
    numby <- 0L
    if(ncustoml) if(is.null(cols)) vl <- FALSE else v <- cols2int(cols, X, nam)
    by <- GRP.default(`names<-`(list(by), l1orlst(as.character(substitute(by)))), NULL, sort, decreasing, na.last, keep.by, return.order, method, call = FALSE)
  } else {
    if(ncustoml) if(is.null(cols)) vl <- FALSE else v <- cols2int(cols, X, nam)
    if(!is_GRP(by)) by <- GRP.default(by, NULL, sort, decreasing, na.last, keep.by, return.order, method, call = FALSE)
    numby <- rep(0L, length(by[[5L]]))
    if(keep.by && !vl && any(m <- nam %in% by[[5L]])) {
      v <- whichv(m, FALSE)
      vl <- TRUE
    }
  }

  if(!nwl) {
    if(is.call(w)) {
      namw <- all.vars(w)
      numw <- ckmatch(namw, nam)
      if(ncustoml) if(vl) v <- fsetdiff(v, numw) else { # v[v != numw]
        v <- nam %!iin% namw; vl <- TRUE
      }
      w <- eval(w[[2L]], X, attr(w, ".Environment")) # w <- X[[numw]]
    } else if(keep.w) {
      numw <- 0L # length(X) + 1L
      namw <- l1orlst(as.character(substitute(w)))
    }
    if(keep.w) { # what about function name for give.names ? What about give.names and custom ?
      wFUN <- acr_get_funs(substitute(wFUN), wFUN, mymatchfun)
      namwFUN <- wFUN$namfun
      wFUN <- wFUN$funs
      if(!all(names(wFUN) %in% .FAST_STAT_FUN_EXT)) stop("wFUN needs to be fast statistical functions, see print(.FAST_STAT_FUN)")
      if(length(wFUN) > 1L) {
        namw <- paste(namwFUN, namw, sep = ".")
        by[[4L]] <- c(if(keep.by) by[[4L]], `names<-`(lapply(wFUN, function(f) f(w, g = by, ..., use.g.names = FALSE)), namw))
        if(keep.col.order) numby <- c(if(keep.by) numby, rep_len(numw, length(wFUN)))
      } else {
        wFUN <- wFUN[[1L]]
        if(isTRUE(give.names)) namw <- paste(namwFUN, namw, sep = ".")
        by[[4L]] <- c(if(keep.by) by[[4L]], `names<-`(list(wFUN(w, g = by, ..., use.g.names = FALSE)), namw))
        if(keep.col.order) numby <- c(if(keep.by) numby, numw)  # need to accommodate any option of keep.by, keep.w and keep.col.order
      }
      keep.by <- TRUE
    }
  }

  if(ncustoml) {

    # Identifying data
    nu <- .Call(C_vtypes, X, 1L) # vapply(unattrib(X), is.numeric, TRUE)
    if(vl) {
      temp <- nu[v]
      nnu <- v[!temp] # which(!nu & v) # faster way ?
      nu <- v[temp] # which(nu & v)
      rm(temp, v)
    } else {
      nnu <- whichv(nu, FALSE)
      nu <- which(nu)
    }
    nul <- length(nu) > 0L
    nnul <- length(nnu) > 0L

    # Identifying FUN and catFUN:
    if(nul) {
      FUN <- acr_get_funs(substitute(FUN), FUN, mymatchfun)
      namFUN <- FUN$namfun
      FUN <- FUN$funs
    }
    if(nnul) {
      catFUN <- acr_get_funs(substitute(catFUN), catFUN, mymatchfun)
      namcatFUN <- catFUN$namfun
      catFUN <- catFUN$funs
    }

    if(autorn) give.names <- !widel || length(FUN) > 1L || length(catFUN) > 1L

    # Aggregator function # drop level of nesting i.e. make rest length(by)+length(FUN)+length(catFUN)  ?
    agg <- function(xnu, xnnu, ...) { # by, FUN, namFUN, catFUN, namcatFUN, drop.by
      lr <- nul + nnul + keep.by
      res <- vector("list", lr)
      if(keep.by) {
        res[[1L]] <- list(by[[4L]]) # could add later using "c" ?
        ind <- 2L
      } else ind <- 1L
      if(nul) res[[ind]] <- condsetn(applyfuns_internal(xnu, by, FUN, names(FUN) %in% .FAST_STAT_FUN_EXT,
                                     parallel, mc.cores, ...), namFUN, give.names)
      if(nnul) res[[lr]] <- condsetn(applyfuns_internal(xnnu, by, catFUN, names(catFUN) %in% .FAST_STAT_FUN_EXT,
                                     parallel, mc.cores, ...), namcatFUN, give.names)
      return(res)
    } # fastest using res list ?? or better combine at the end ??

    # Fixes https://github.com/SebKrantz/collapse/issues/185
    if(widel && !give.names && ((length(nu) == 1L && !nnul && length(FUN) > 1L) || (length(nnu) == 1L && !nul && length(catFUN) > 1L))) {
       names(X) <- NULL
       give.names <- TRUE
    }

    if(nwl) {
      res <- agg(if(nul) X[nu] else NULL,
                 if(nnul) X[nnu] else NULL, ...)
    } else {
      res <- agg(if(nul) X[nu] else NULL,
                 if(nnul) X[nnu] else NULL, w = w, ...)
    }


    if(keep.col.order && widel) o <- forder.int(c(if(keep.by) numby else NULL,
                                                  if(nul) rep(nu,length(FUN)) else NULL,
                                                  if(nnul) rep(nnu,length(catFUN)) else NULL))

  } else { # custom aggregation:

    namFUN <- names(custom)
    if(!is.list(custom) || is.null(namFUN)) stop("custom needs to be a named list, see ?collap")
    fFUN <- namFUN %in% .FAST_STAT_FUN_EXT
    if(!keep.by) {
      res <- vector("list", 1L)
      ind <- 1L
    } else {
      res <- vector("list", 2L)
      res[[1L]] <- list(by[[4L]]) # could add later using "c" ?
      ind <- 2L
    }

    custom_names <- lapply(custom, names)
    custom <- lapply(custom, cols2int, X, nam) # could integrate below, but then reorder doesn't work !

    # if(autorn) give.names <- fanyDuplicated(funlist(custom))
    #lx <- length(X)
    # custom <- lapply(custom, function(x) if(is.numeric(x) && bmax(abs(x)) <= lx)
    #                          x else if(is.character(x)) ckmatch(x, nam) else
    #                          stop("custom list content must be variable names or suitable column indices"))

    if(nwl) {
      res[[ind]] <- lapply(seq_along(namFUN), function(i)
                      applyfuns_internal(setnck(X[custom[[i]]], custom_names[[i]]), by, mymatchfun(namFUN[i]),
                                         fFUN[i], parallel, mc.cores, ...)[[1L]])
    } else {
      if(!all(fFUN)) warning("collap can only perform weighted aggregations with the fast statistical functions (see .FAST_STAT_FUN): Ignoring weights argument to other functions")
      res[[ind]] <- lapply(seq_along(namFUN), function(i)
                      applyfuns_internal(setnck(X[custom[[i]]], custom_names[[i]]), by, mymatchfun(namFUN[i]),
                                         fFUN[i], parallel, mc.cores, w = w, ...)[[1L]])
    }
    # Better to do this check afterwards, because custom names may make column names unique...
    if(autorn && widel) give.names <- fanyDuplicated(funlist(lapply(res[[ind]], attr, "names")))
    if(!widel || give.names) names(res[[ind]]) <- namFUN

    if(keep.col.order && return != 2L) { # && widel
      o <- unlist(custom, use.names = FALSE)
      o <- forder.int(c(if(keep.by) numby else NULL, if(widel) o else unique.default(o)))
    }
  }

  # if(parallel) stopCluster(cl)
  if(widel) res <- unlist(unlist(res, FALSE), FALSE) else {
    # if(length(FUN) > 1L || length(catFUN) > 1L || length(custom) > 1L) {
      res <- unlist(res, FALSE)
      if(return == 2L) {
        ax[["row.names"]] <- .set_row_names(by[[1L]])
        if(!keep.by) return(lapply(res, function(e) {
                            ax[["names"]] <- names(e)
                            condalcSA(e, ax, DTl) }))
        namby <- attr(res[[1L]], "names") # always works ??
        return(lapply(res[-1L], function(e) {
          ax[["names"]] <- c(namby, names(e))
          condalcSA(c(res[[1L]], e), ax, DTl) }))
      } else {
        if(return != 4L) {
          if(keep.by) res <- lapply(res[-1L], function(e) c(res[[1L]], e))
        } else {
          if(!ncustoml || !(nul && nnul)) stop("long_dupl is only meaningful for aggregations with both numeric and categorical data, and multiple functions used for only one of the two data types!")
          mFUN <- length(FUN) > 1L
          nid <- if(mFUN) length(res) else 2L-!keep.by
          if(!keep.by) {
            res <- if(mFUN) lapply(res[-nid], function(e) c(e, res[[nid]])) else
              lapply(res[-nid], function(e) c(res[[nid]], e))
          } else res <- if(mFUN) lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], e, res[[nid]])) else
            lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], res[[nid]], e))
        }
        res <- rbindlist_factor(res)
        if(keep.col.order)  o <- if(ncustoml) forder.int(c(0L, if(keep.by) numby else NULL, nu, nnu)) else c(1L, o + 1L)
      }
    # } else message("return options other than 'wide' are only meaningful if multiple functions are used!")
  }

  if(keep.col.order) .Call(C_setcolorder, res, o) # data.table:::Csetcolorder
  ax[["names"]] <- names(res)
  ax[["row.names"]] <- .set_row_names(.Call(C_fnrow, res))
  return(condalcSA(res, ax, DTl))
}


# collapv: allows vector input to by and w
collapv <- function(X, by, FUN = fmean, catFUN = fmode, cols = NULL, w = NULL, wFUN = fsum, custom = NULL, ...,
                    keep.by = TRUE, keep.w = TRUE, keep.col.order = TRUE, sort = .op[["sort"]], decreasing = FALSE,
                    na.last = TRUE, return.order = sort, method = "auto", parallel = FALSE, mc.cores = 2L,
                    return = c("wide","list","long","long_dupl"), give.names = "auto") {

  return <- switch(return[1L], wide = 1L, list = 2L, long = 3L, long_dupl = 4L, stop("Unknown return output option"))
  widel <- return == 1L
  ncustoml <- is.null(custom)
  autorn <- is.character(give.names) && give.names == "auto"
  nwl <- is.null(w)
  if(inherits(X, "data.frame")) DTl <- inherits(X, "data.table") else {
    X <- qDF(X)
    DTl <- FALSE
  }
  ax <- attributes(X)
  oldClass(X) <- NULL
  if(.Call(C_fnrow, X) == 0L) stop("data passed to collapv() has 0 rows.") #160, 0 rows can cause segfault...
  nam <- names(X)

  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply

  # identifying by
  numby <- cols2int(by, X, nam)
  by <- GRP.default(X, numby, sort, decreasing, na.last, keep.by, return.order, method, call = FALSE)
  if(ncustoml) v <- if(is.null(cols)) seq_along(X)[-numby] else cols2int(cols, X, nam)


  if(!nwl) {
    if(length(w) == 1L) {
      numw <- cols2int(w, X, nam)
      namw <- nam[numw]
      if(ncustoml) v <- v[v != numw]
      w <- X[[numw]]
    } else if(keep.w) {
      numw <- 0L
      namw <- l1orlst(as.character(substitute(w)))
    }
    if(keep.w) {
      wFUN <- acr_get_funs(substitute(wFUN), wFUN, mymatchfun)
      namwFUN <- wFUN$namfun
      wFUN <- wFUN$funs
      if(!all(names(wFUN) %in% .FAST_STAT_FUN_EXT)) stop("wFUN needs to be fast statistical functions, see print(.FAST_STAT_FUN)")
      if(length(wFUN) > 1L) {
        namw <- paste(namwFUN, namw, sep = ".")
        by[[4L]] <- c(if(keep.by) by[[4L]], `names<-`(lapply(wFUN, function(f) f(w, g = by, ..., use.g.names = FALSE)), namw))
        if(keep.col.order) numby <- c(if(keep.by) numby, rep_len(numw, length(wFUN)))
      } else {
        wFUN <- wFUN[[1L]]
        if(isTRUE(give.names)) namw <- paste(namwFUN, namw, sep = ".")
        by[[4L]] <- c(if(keep.by) by[[4L]], `names<-`(list(wFUN(w, g = by, ..., use.g.names = FALSE)), namw))
        if(keep.col.order) numby <- c(if(keep.by) numby, numw)  # need to accommodate any option of keep.by, keep.w and keep.col.order
      }
      keep.by <- TRUE
    }
  }

  if(ncustoml) {

    # Identifying data
    nu <- .Call(C_vtypes, X, 1L) # vapply(unattrib(X), is.numeric, TRUE)
    temp <- nu[v]
    nnu <- v[!temp] # which(!nu & v) # faster way ?
    nu <- v[temp] # which(nu & v)
    rm(temp, v)
    nul <- length(nu) > 0L
    nnul <- length(nnu) > 0L

    # Identifying FUN and catFUN:
    if(nul) {
      FUN <- acr_get_funs(substitute(FUN), FUN, mymatchfun)
      namFUN <- FUN$namfun
      FUN <- FUN$funs
    }
    if(nnul) {
      catFUN <- acr_get_funs(substitute(catFUN), catFUN, mymatchfun)
      namcatFUN <- catFUN$namfun
      catFUN <- catFUN$funs
    }

    if(autorn) give.names <- !widel || length(FUN) > 1L || length(catFUN) > 1L

    # Aggregator function
    agg <- function(xnu, xnnu, ...) {
      lr <- nul + nnul + keep.by
      res <- vector("list", lr)
      if(keep.by) {
        res[[1L]] <- list(by[[4L]])
        ind <- 2L
      } else ind <- 1L
      if(nul) res[[ind]] <- condsetn(applyfuns_internal(xnu, by, FUN, names(FUN) %in% .FAST_STAT_FUN_EXT,
                                                        parallel, mc.cores, ...), namFUN, give.names)
      if(nnul) res[[lr]] <- condsetn(applyfuns_internal(xnnu, by, catFUN, names(catFUN) %in% .FAST_STAT_FUN_EXT,
                                                        parallel, mc.cores, ...), namcatFUN, give.names)
      return(res)
    }

    # Fixes https://github.com/SebKrantz/collapse/issues/185
    if(widel && !give.names && ((length(nu) == 1L && !nnul && length(FUN) > 1L) || (length(nnu) == 1L && !nul && length(catFUN) > 1L))) {
      names(X) <- NULL
      give.names <- TRUE
    }

    if(nwl) {
      res <- agg(if(nul) X[nu] else NULL,
                 if(nnul) X[nnu] else NULL, ...)
    } else {
      res <- agg(if(nul) X[nu] else NULL,
                 if(nnul) X[nnu] else NULL, w = w, ...)
    }

    if(keep.col.order && widel) o <- forder.int(c(if(keep.by) numby else NULL,
                                                  if(nul) rep(nu,length(FUN)) else NULL,
                                                  if(nnul) rep(nnu,length(catFUN)) else NULL))

  } else { # custom aggregation:

    namFUN <- names(custom)
    if(!is.list(custom) || is.null(namFUN)) stop("custom needs to be a named list, see ?collap")
    fFUN <- namFUN %in% .FAST_STAT_FUN_EXT
    if(!keep.by) {
      res <- vector("list", 1L)
      ind <- 1L
    } else {
      res <- vector("list", 2L)
      res[[1L]] <- list(by[[4L]])
      ind <- 2L
    }

    custom_names <- lapply(custom, names)
    custom <- lapply(custom, cols2int, X, nam)

    if(nwl) {
      res[[ind]] <- lapply(seq_along(namFUN), function(i)
        applyfuns_internal(setnck(X[custom[[i]]], custom_names[[i]]), by, mymatchfun(namFUN[i]),
                           fFUN[i], parallel, mc.cores, ...)[[1L]])
    } else {
      if(!all(fFUN)) warning("collap can only perform weighted aggregations with the fast statistical functions (see .FAST_STAT_FUN): Ignoring weights argument to other functions")
      res[[ind]] <- lapply(seq_along(namFUN), function(i)
        applyfuns_internal(setnck(X[custom[[i]]], custom_names[[i]]), by, mymatchfun(namFUN[i]),
                           fFUN[i], parallel, mc.cores, w = w, ...)[[1L]])
    }
    # Better to do this check afterwards, because custom names may make column names unique...
    if(autorn && widel) give.names <- fanyDuplicated(funlist(lapply(res[[ind]], attr, "names")))
    if(!widel || give.names) names(res[[ind]]) <- namFUN

    if(keep.col.order && return != 2L) {
      o <- unlist(custom, use.names = FALSE)
      o <- forder.int(c(if(keep.by) numby else NULL, if(widel) o else unique.default(o)))
    }
  }

  if(widel) res <- unlist(unlist(res, FALSE), FALSE) else {
    # if(length(FUN) > 1L || length(catFUN) > 1L || length(custom) > 1L) {
      res <- unlist(res, FALSE)
      if(return == 2L) {
        ax[["row.names"]] <- .set_row_names(by[[1L]])
        if(!keep.by) return(lapply(res, function(e) {
                            ax[["names"]] <- names(e)
                            condalcSA(e, ax, DTl) }))
        namby <- attr(res[[1L]], "names") # always works ??
        return(lapply(res[-1L], function(e) {
          ax[["names"]] <- c(namby, names(e))
          condalcSA(c(res[[1L]], e), ax, DTl) }))
      } else {
        if(return != 4L) {
          if(keep.by) res <- lapply(res[-1L], function(e) c(res[[1L]], e))
        } else {
          if(!ncustoml || !(nul && nnul)) stop("long_dupl is only meaningful for aggregations with both numeric and categorical data, and multiple functions used for only one of the two data types!")
          mFUN <- length(FUN) > 1L
          nid <- if(mFUN) length(res) else 2L-!keep.by
          if(!keep.by) {
            res <- if(mFUN) lapply(res[-nid], function(e) c(e, res[[nid]])) else
              lapply(res[-nid], function(e) c(res[[nid]], e))
          } else res <- if(mFUN) lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], e, res[[nid]])) else
            lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], res[[nid]], e))
        }
        res <- rbindlist_factor(res)
        if(keep.col.order)  o <- if(ncustoml) forder.int(c(0L, if(keep.by) numby else NULL, nu, nnu)) else c(1L, o + 1L)
      }
    # } else message("return options other than 'wide' are only meaningful if multiple functions are used!")
  }

  if(keep.col.order) .Call(C_setcolorder, res, o) # data.table:::Csetcolorder
  ax[["names"]] <- names(res)
  ax[["row.names"]] <- .set_row_names(.Call(C_fnrow, res))
  return(condalcSA(res, ax, DTl))
}


# For dplyr integration: takes grouped_df as input
collapg <- function(X, FUN = fmean, catFUN = fmode, cols = NULL, w = NULL, wFUN = fsum, custom = NULL,
                    keep.group_vars = TRUE, ...) {
  by <- GRP.grouped_df(X, return.groups = keep.group_vars, call = FALSE)
  if(is.null(by[[4L]])) keep.group_vars <- FALSE
  if(is.null(custom)) ngn <- attr(X, "names") %!in% by[[5L]] # Note: this always leaves grouping columns on the left still !
  # clx <- oldClass(X)
  attr(X, "groups") <- NULL
  oldClass(X) <- fsetdiff(oldClass(X), c("GRP_df", "grouped_df"))  # clx[clx != "grouped_df"]
  wsym <- substitute(w)
  if(!is.null(wsym)) { # Non-standard evaluation of w argument
    if(any(windl <- attr(X, "names") %in% all.vars(wsym))) {
      wchar <- if(length(wsym) == 1L) as.character(wsym) else deparse(wsym)
      assign(wchar,  eval(wsym, X, parent.frame())) # needs to be here !! (before subsetting!!)
      if(is.null(custom)) X <- fcolsubset(X, ngn & !windl) # else X <- X # Needed ?? -> nope !!
      expr <- substitute(collap(X, by, FUN, catFUN, cols, w, wFUN, custom, ...,
                                keep.by = keep.group_vars,
                                sort = TRUE,
                                decreasing = FALSE,
                                na.last = TRUE,
                                return.order = TRUE,
                                method = "auto"))
      expr[[7L]] <- as.symbol(wchar) # best solution !!
      return(eval(expr))
    }
  }
  if(is.null(custom)) X <- fcolsubset(X, ngn) # else X <- X # because of non-standard eval.. X is "."
  return(eval(substitute(collap(X, by, FUN, catFUN, cols, w, wFUN, custom, ...,
                                keep.by = keep.group_vars,
                                sort = TRUE,
                                decreasing = FALSE,
                                na.last = TRUE,
                                return.order = TRUE,
                                method = "auto"))))
}
