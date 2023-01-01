
collap <- function(X, by, FUN = fmean, catFUN = fmode, cols = NULL, custom = NULL,
                   keep.by = TRUE, keep.col.order = TRUE, sort.row = TRUE,
                   parallel = FALSE, mc.cores = 1L,
                   return = c("wide","list","long","long_dupl"), give.names = "auto", ...) {

  return <- switch(return[1L], wide = 1L, list = 2L, long = 3L, long_dupl = 4L, stop("Unknown return output option"))
  widel <- return == 1L
  customl <- !is.null(custom)
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
        if(!customl) {
          v <- if(is.null(cols)) !logical(length(X)) else cols2log(cols, X, nam)
          v[numby] <- FALSE
        }
      }
    by <- GRP.default(X, numby, sort = sort.row, return.groups = keep.by)
  } else if(is.atomic(by)) {
    namby <- deparse(substitute(by))
    numby <- 1L
    if(!customl) if(is.null(cols)) vl <- FALSE else v <- cols2log(cols, X, nam)
    by <- GRP.default(by, sort = sort.row, return.groups = keep.by)
  } else {
    if(!customl) if(is.null(cols)) vl <- FALSE else v <- cols2log(cols, X, nam)
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

  if(!customl) {

    # identifying data
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
    if(nul) if(is.character(FUN)) {
      # FUN <- unlist(strsplit(FUN,",",fixed = TRUE), use.names = FALSE)
      namFUN <- FUN
      FUN <- if(length(FUN) > 1L) lapply(FUN, match.fun, descend = FALSE) else
                                  match.fun(FUN, descend = FALSE)
    } else if(is.list(FUN)) {
      namFUN <- names(FUN)
      if(is.null(namFUN)) namFUN <- all.vars(substitute(FUN))
    } else namFUN <- l1orn(as.character(substitute(FUN)), "FUN") # Faster !

    if(nnul) if(is.character(catFUN)) {
      # catFUN <- unlist(strsplit(catFUN,",",fixed = TRUE), use.names = FALSE)
      namcatFUN <- catFUN
      catFUN <- if(length(catFUN) > 1L) lapply(catFUN, match.fun, descend = FALSE) else
                                        match.fun(catFUN, descend = FALSE)
    } else if(is.list(catFUN)) {
      namcatFUN <- names(catFUN)
      if(is.null(namcatFUN)) namcatFUN <- all.vars(substitute(catFUN))
    } else namcatFUN <- l1orn(as.character(substitute(catFUN)), "catFUN") # Faster !

    if(give.names == "auto") give.names <- !widel || length(FUN) > 1L || length(catFUN) > 1L

    # Aggregator function # drop level of nesting i.e. make rest length(by)+length(FUN)+length(catFUN)  ?
    agg <- function(xnu, xnnu) { # by, FUN, namFUN, catFUN, namcatFUN, drop.by
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
    } # fastest isung res list ?? or better combine at the end ??
    res <- agg(if(nul) `oldClass<-`(X[nu], "data.frame") else NULL, if(nnul) `oldClass<-`(X[nnu], "data.frame") else NULL)

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

    res[[ind]] <- condsetn(aplyfun(seq_along(namFUN), function(i)
            if(fFUN[i]) match.fun(namFUN[i])(`oldClass<-`(X[custom[[i]]], "data.frame"), by, ..., use.g.names = FALSE) else
            BY.data.frame(X[custom[[i]]], by, namFUN[i], ..., use.g.names = FALSE)), namFUN, give.names)
    if(keep.col.order && widel) {
      o <- unlist(custom, use.names = FALSE)
      o <- forder.int(c(if(!keep.by) NULL else if(!bycalll) rep(0L,length(numby)) else numby, o))
    }
  }
  if(widel) res <- unlist(unlist(res, FALSE), FALSE) else {
    if(length(FUN) > 1L || length(catFUN) > 1L || length(custom) > 1L) {
      res <- unlist(res, FALSE)
      if(return == 2L) {
        ax[["row.names"]] <- if(is.list(by)) .set_row_names(by[[1L]]) else .set_row_names(length(res[[1L]]))
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
          if(!(nul && nnul) || customl) stop("long_dupl is only meaningful for aggregations with both numeric and categorical data, and multiple functions used for only one of the two data types!")
          mFUN <- length(FUN) > 1L
          nid <- if(mFUN) length(res) else 2L-!keep.by
          if(!keep.by) {
            res <- if(mFUN) lapply(res[-nid], function(e) c(e, res[[nid]])) else
                            lapply(res[-nid], function(e) c(res[[nid]], e))
          } else res <- if(mFUN) lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], e, res[[nid]])) else
                                 lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], res[[nid]], e))
          res <- .Call(C_rbindlist, res, FALSE, FALSE, "Function")
        }
        if(keep.col.order)  o <- forder.int(c(0L, if(!keep.by) NULL else if(!bycalll) rep(0L,length(numby)) else numby, nu, nnu))
      }
    } else message("return options other than 'wide' are only meaningful if multiple functions are used!")
  }

  if(keep.col.order) .Call(C_setcolorder, res, o) # data.table:::Csetcolorder
  ax[["names"]] <- names(res)
  ax[["row.names"]] <- .set_row_names(length(res[[1L]]))
  return(setAttributes(res, ax))
}

# collapv: allows vector input to by !
collapv <- function(X, by, FUN = fmean, catFUN = fmode, cols = NULL, custom = NULL,
                   keep.by = TRUE, keep.col.order = TRUE, sort.row = TRUE,
                   parallel = FALSE, mc.cores = 1L,
                   return = c("wide","list","long","long_dupl"), give.names = "auto", ...) {

  return <- switch(return[1L], wide = 1L, list = 2L, long = 3L, long_dupl = 4L, stop("Unknown return output option"))
  widel <- return == 1L
  customl <- !is.null(custom)
  if(!inherits(X, "data.frame")) X <- qDF(X)
  ax <- attributes(X)
  class(X) <- NULL
  nam <- names(X)
  # attributes(X) <- NULL
  # attr(X, "class") <- "data.frame" # class needed for method dispatch of fast functions, not for BY !

  aplyfun <- if(parallel) function(...) mclapply(..., mc.cores = mc.cores) else lapply

  # identifying by
  numby <- cols2int(by, X, nam)
  namby <- nam[numby]
  by <- GRP.default(X, numby, sort = sort.row, return.groups = keep.by)

  if(!customl) {

    v <- if(is.null(cols)) !logical(length(X)) else cols2log(cols, X, nam)
    v[numby] <- FALSE

    # identifying data
    nu <- vapply(unattrib(X), is.numeric, TRUE)
    nnu <- which(!nu & v) # faster way ?
    nu <- which(nu & v)
    nul <- length(nu) > 0L
    nnul <- length(nnu) > 0L

    # Identifying FUN and catFUN:
    if(nul) if(is.character(FUN)) {
      # FUN <- unlist(strsplit(FUN,",",fixed = TRUE), use.names = FALSE)
      namFUN <- FUN
      FUN <- if(length(FUN) > 1L) lapply(FUN, match.fun, descend = FALSE) else
        match.fun(FUN, descend = FALSE)
    } else if(is.list(FUN)) {
      namFUN <- names(FUN)
      if(is.null(namFUN)) namFUN <- all.vars(substitute(FUN))
    } else namFUN <- l1orn(as.character(substitute(FUN)), "FUN") # Faster !

    if(nnul) if(is.character(catFUN)) {
      # catFUN <- unlist(strsplit(catFUN,",",fixed = TRUE), use.names = FALSE)
      namcatFUN <- catFUN
      catFUN <- if(length(catFUN) > 1L) lapply(catFUN, match.fun, descend = FALSE) else
        match.fun(catFUN, descend = FALSE)
    } else if(is.list(catFUN)) {
      namcatFUN <- names(catFUN)
      if(is.null(namcatFUN)) namcatFUN <- all.vars(substitute(catFUN))
    } else namcatFUN <- l1orn(as.character(substitute(catFUN)), "catFUN") # Faster !

    if(give.names == "auto") give.names <- !widel || length(FUN) > 1L || length(catFUN) > 1L

    # Aggregator function # drop level of nesting i.e. make rest length(by)+length(FUN)+length(catFUN)  ?
    agg <- function(xnu, xnnu) { # by, FUN, namFUN, catFUN, namcatFUN, drop.by
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
    } # fastest isung res list ?? or better combine at the end ??
    res <- agg(if(nul) `oldClass<-`(X[nu], "data.frame") else NULL, if(nnul) `oldClass<-`(X[nnu], "data.frame") else NULL)

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
      res[[1L]] <- list(by[[4L]]) # could add later using "c" ?
      ind <- 2L
    }
    custom <- lapply(custom, cols2int, X, nam) # could integrate below, but then reorder doesn't work !
    #lx <- length(X)
    # custom <- lapply(custom, function(x) if(is.numeric(x) && max(abs(x)) <= lx)
    #                          x else if(is.character(x)) ckmatch(x, nam) else
    #                          stop("custom list content must be variable names or suitable column indices"))

    res[[ind]] <- condsetn(aplyfun(seq_along(namFUN), function(i)
      if(fFUN[i]) match.fun(namFUN[i])(`oldClass<-`(X[custom[[i]]], "data.frame"), by, ..., use.g.names = FALSE) else
        BY.data.frame(X[custom[[i]]], by, namFUN[i], ..., use.g.names = FALSE)), namFUN, give.names)
    if(keep.col.order && widel) {
      o <- unlist(custom, use.names = FALSE)
      o <- forder.int(c(if(!keep.by) NULL else numby, o))
    }
  }
  if(widel) res <- unlist(unlist(res, FALSE), FALSE) else {
    if(length(FUN) > 1L || length(catFUN) > 1L || length(custom) > 1L) {
      res <- unlist(res, FALSE)
      if(return == 2L) {
        ax[["row.names"]] <- if(is.list(by)) .set_row_names(by[[1L]]) else .set_row_names(length(res[[1L]]))
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
          if(!(nul && nnul) || customl) stop("long_dupl is only meaningful for aggregations with both numeric and categorical data, and multiple functions used for only one of the two data types!")
          mFUN <- length(FUN) > 1L
          nid <- if(mFUN) length(res) else 2L-!keep.by
          if(!keep.by) {
            res <- if(mFUN) lapply(res[-nid], function(e) c(e, res[[nid]])) else
              lapply(res[-nid], function(e) c(res[[nid]], e))
          } else res <- if(mFUN) lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], e, res[[nid]])) else
            lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], res[[nid]], e))
          res <- .Call(C_rbindlist, res, FALSE, FALSE, "Function")
        }
        if(keep.col.order)  o <- forder.int(c(0L, if(!keep.by) NULL else numby, nu, nnu))
      }
    } else message("return options other than 'wide' are only meaningful if multiple functions are used!")
  }

  if(keep.col.order) .Call(C_setcolorder, res, o) # data.table:::Csetcolorder
  ax[["names"]] <- names(res)
  ax[["row.names"]] <- .set_row_names(length(res[[1L]]))
  return(setAttributes(res, ax))
}

# For dplyr integration: takes grouped_df as input
collapg <- function(X, FUN = fmean, catFUN = fmode, cols = NULL, custom = NULL,
                    keep.group_vars = TRUE, keep.col.order = TRUE, sort.row = TRUE,
                    parallel = FALSE, mc.cores = 1L,
                    return = c("wide","list","long","long_dupl"), give.names = "auto", ...) {
  by <- GRP.grouped_df(X)
  ngn <- attr(X, "names") %!in% by[[5L]] # Note: this always leaves grouping columns on the left still !
  clx <- class(X)
  attr(X, "groups") <- NULL
  oldClass(X) <- clx[clx != "grouped_df"]
  if(is.function(FUN)) FUN <- `names<-`(list(FUN), l1orn(as.character(substitute(FUN)), "FUN")) else
    if(is.list(FUN) && is.null(names(FUN))) names(FUN) <- all.vars(substitute(FUN))
  if(is.function(catFUN)) catFUN <- `names<-`(list(catFUN), l1orn(as.character(substitute(catFUN)), "catFUN")) else
    if(is.list(catFUN) && is.null(names(catFUN))) names(catFUN) <- all.vars(substitute(catFUN))

  return(collap(fcolsubset(X, ngn), by, FUN, catFUN, cols, custom,
                keep.group_vars, keep.col.order, sort.row, parallel,
                mc.cores, return, give.names, ...))
}
