# note: plyr has a function and class indexed_df...

# getpix <- function(x) switch(typeof(x), externalptr = .Call(C_geteptr, x), x)

findex <- function(x) {
  idx <- attr(x, "index_df")
  if(is.null(idx)) idx <- attr(x, "index")
  if(is.list(idx)) return(idx)
  .Call(C_geteptr, idx)
}

ix <- findex

# TODO use attr(ids, "optim_time") ? -> think about what is the smartest way to implement this. Also think in the long-term
# how further optimization (e.g. ordering vector) will take effect...
# also what about sorted data ?? If regular panel should be able to optimize... i.e. compute without time index.

to_plm <- function(x, row.names = FALSE) {
  index <- unclass(findex(x))
  if(is.null(index)) stop("Missing index!")
  if(length(index) < 2L) stop("plm compatible index must have at least 2 factors")
  if(length(index) > 2L) index <- c(list(id = finteraction(index[-length(index)], sort = FALSE)), index[length(index)])
  oldClass(index[[1L]]) <- "factor"
  oldClass(index[[2L]]) <- "factor"
  attr(index, "row.names") <- .set_row_names(length(index[[1L]]))
  oldClass(index) <- c("pindex", "data.frame")

  if(is.list(x) && inherits(x, "indexed_frame")) {
    res <- qDF(unindex(x), class = c("pdata.frame", "data.frame"))
    attr(res, "index") <- index
    if(row.names) attr(res, "row.names") <- do.call(paste, c(index, list(sep = "-")))
  } else if(inherits(x, "indexed_series")) {
    res <- unindex(x)
    attr(res, "index") <- index
    oldClass(res) <- c("pseries", class(res))
    if(row.names) names(res) <- do.call(paste, c(index, list(sep = "-")))
  } else stop("x must be 'indexed_frame' or 'indexed_series'")

  return(res)
}


# # fixest:
# time = unclass(wlddev$date)
# time_full = fixest:::quickUnclassFactor(time, addItem = TRUE, sorted = TRUE)
# time_unik = time_full$items
# all_steps = unique(diff(time_unik))
# my_step = fixest:::cpp_pgcd(all_steps)
# # we rescale time_unik
# time_unik_new = (time_unik - min(time_unik)) / my_step
# time = time_unik_new[time_full$x]

# TODO: also think of fixest's quf, checking if double is integer, break out of loop if fail...

timeid <- function(x, factor = FALSE, ordered = factor, extra = FALSE) {
  id <- .Call(C_group, x, TRUE, FALSE) # starts = TRUE, group.sizes = FALSE
  unik <- Csv(x, attr(id, "starts"))
  attributes(unik) <- NULL
  if(!is.numeric(unik)) stop("x needs to be numeric, otherwise use qF() or qG() instead of timeid()")
  is_dbl <- is.double(unik)
  o <- forder.int(unik, na.last = TRUE)
  ng <- length(o)
  unik_o <- if(attr(o, "sorted") && !extra) unik else Csv(unik, o) # !extra because of math by reference...
  if(is.na(unik_o[ng])) stop("Time variable may not contain missing values")
  r <- c(unik_o[1L], unik_o[ng])
  steps <- unik_o[-1L] %-=% unik_o[-ng] # tsibble uses abs(diff(unik_o)), but here we sort the values, so not necessary
  # if(is_dbl) steps <- round(steps, digits = 6) # This is pretty costly for long POSIXct sequences. Better not do it..
  gcd <- .Call(C_vecgcd, .Call(Cpp_sortunique, steps))
  if(is_dbl) {
    if(r[1L] != 1 || gcd != 1) unik %-=% (r[1L] - 1.4*gcd) # * 1.4 to make sure the as.integer conversion does proper rounding
    if(gcd != 1) unik %/=% gcd
    unik <- as.integer(unik)
  } else {
    if(r[1L] != 1L || gcd != 1L) unik %-=% (r[1L] - gcd)
    if(gcd != 1L) unik %/=% gcd
  }
  tid <- if(length(id) == ng) unik else Csv(unik, id)

  if(factor) {
    levnum <- if(is_dbl) seq.default(r[1L], r[2L]+0.4*gcd, gcd) else if(gcd == 1L) r[1L]:r[2L] else seq.int(r[1L], r[2L], gcd)
    if(is.object(x)) levnum <- copyMostAttrib(levnum, x)
    attr(tid, "levels") <- if(is.object(x) && is_date(x)) strftime(levnum, format = if(inherits(x, "Date")) "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S %Z") else as.character(levnum)
    oldClass(tid) <- c(if(ordered) "ordered", "factor", "na.included")
  } else {
    attr(tid, "N.groups") <- as.integer(if(is_dbl) (r[2L]+0.4*gcd-r[1L])/gcd else (r[2L]-r[1L])/gcd) + 1L
    oldClass(tid) <- c(if(ordered) "ordered", "qG", "na.included")
  }

  if(extra) {
    attr(tid, "unique_ints") <- unik
    attr(tid, "sort_unique_x") <- copyMostAttrib(unik_o, x)
    attr(tid, "range_x") <- copyMostAttrib(r, x)
    attr(tid, "step_x") <- gcd
  }
  tid
}

make_time_factor <- function(x) {
  if(inherits(x, c("factor", "qG"))) { # Make sure we handle irregularity correctly...
    if(is_qG(x)) return(as_factor_qG(x, na.exclude = FALSE))
    if(inherits(x, "na.included")) return(x)
    if(anyNA(x)) stop("Time variable may not contain missing values")
    oldClass(x) <- c("factor", "na.included")
    return(x)
  }
  if(is.numeric(x) && !is.object(x)) {
    idbl <- is.double(x)
    if(idbl) {
      # message("Time variable is of type double, but not a date/time object. It is therefore coerced to integer and assumed to represent unitary timesteps. If this is not desired pass timeid(t). To silence this message pass as.integer(t).")
      x <- as.integer(x)
    }
    r <- .Call(C_frange, x, FALSE) # na.rm = FALSE # Note that inside flag() and fgrowth() etc. we subtract the minimum within each group...
    if(anyNA(r)) stop("Time variable may not contain missing values")
    if(r[1L] != 1) {
      if(idbl) x %-=% (r[1L] - 1L)
      else x <- x - (r[1L] - 1L) # This is unfortunately quite a bit slower...
    }
    attr(x, "levels") <- as.character(r[1L]:r[2L])
    oldClass(x) <- c("ordered", "factor", "na.included")
    return(x)
  }
  if(is.numeric(unclass(x))) return(timeid(x, factor = TRUE, ordered = FALSE))
  qF(x, na.exclude = FALSE, sort = TRUE, method = "hash")
}

is_irregular <- function(x, any_id = TRUE) {
  if(is.object(x) && inherits(x, c("indexed_frame", "indexed_series"))) x <- findex(x)
  if(is.list(x) && inherits(x, "pindex")) {
    oldClass(x) <- NULL
    if(length(x) > 1L) {
      g <- if(length(x) <= 2L) x[[1L]] else if(any_id)
        group(x[-length(x)]) else finteraction(x[-length(x)], sort = FALSE)
      t <- x[[length(x)]]
      # if(!is.nmfactor(t)) stop("t must be a factor without any missing values")
      attributes(t) <- NULL
      rng_t <- fmax(t, g, use.g.names = !any_id)
      rng_t %-=% fmin(t, g, use.g.names = FALSE)
      rng_t %+=% 1L
      n_t <- fnobs(t, g, use.g.names = FALSE)
      if(any_id) return(!identical(rng_t, n_t))
      res <- rng_t != n_t
      names(res) <- names(rng_t)
      return(res)
    } else if(length(x) == 1L) {
      if(!isFALSE(attr(x, "single.id"))) stop("Index does not contain a time variable")
      t <- x[[1L]]
      # if(!is.nmfactor(t)) stop("t must be a factor without any missing values")
      return(fnlevels(t) != fndistinct(t))
    } else stop("Index has zero length")
  }
  if(!(is.atomic(x) && is.numeric(unclass(x)))) stop("x needs to be an 'indexed_frame', 'indexed_series' or 'pindex' object, or an atomic vector with storage type integer or double.")
  if(is.object(x)) {
    if(is.factor(x)) return(fnlevels(x) != fndistinct(x))
    if(is_qG(x)) return(attr(x, "N.groups") != fndistinct(as_factor_qG(x)))
  }
  attributes(x) <- NULL
  tid <- timeid(x, factor = FALSE, extra = TRUE)
  return(attr(tid, "N.groups") != length(attr(tid, "unique_ints")))
}

# Note: data returned as plain list with attributes !
index_series <- function(data, index, cl) {
  oldClass(data) <- NULL
  iptr <- .Call(C_createeptr, index)
  indexfun <- function(x) {
    attr(x, "index_df") <- iptr
    oldClass(x) <- unique.default(c("indexed_series", "pseries", class(x))) # Use OldClass??
    # class is better for methods such as as.data.frame.numeric (used inside plm) to apply..
    x
  }
  if(any(cl == "sf")) {
    geom <- whichv(names(data), attr(data, "sf_column"))
    data[-geom] <- lapply(data[-geom], indexfun)
    return(data)
  }
  data[] <- lapply(unattrib(data), indexfun) # dapply(data, indexfun)
  data
}

reindex <- function(x, index = findex(x), single = "auto") {
  n <- if(is.list(x)) fnrow(x) else NROW(x)
  if(is.atomic(index)) {
    if(length(index) != n) stop("index does not match data length")
    nam <- l1orlst(as.character(substitute(index)))
    idl <- switch(single, auto = anyDuplicated.default(index) > 0L, id = TRUE, time = FALSE, stop("'single' must be 'auto', 'id' or 'time'"))
    index <- list(if(idl) qF(index, sort = is.factor(index), na.exclude = FALSE) else make_time_factor(index))
    names(index) <- nam
    attr(index, "row.names") <- .set_row_names(n)
    attr(index, "single.id") <- idl
    oldClass(index) <- c("index_df", "pindex", "data.frame")
  } else {
    if(fnrow(index) != n) stop("index does not match data length")
    if(!inherits(index, "pindex")) {
      if(!all(.Call(C_vtypes, index, 2L))) stop("All variables in a valid index must be factors. Please prepare you data accordingly.")
      index <- qDF(index)
      if(fncol(index) == 1L)
        attr(index, "single.id") <- switch(single, auto = anyDuplicated.default(.subset2(index, 1L)) > 0L, id = TRUE, time = FALSE, stop("'single' must be 'auto', 'id' or 'time'"))
      oldClass(index) <- c("index_df", "pindex", "data.frame")
    }
  }
  if(is.list(x)) {
    clx <- oldClass(x)
    x <- index_series(x, index, clx) # x is list afterwards, so need to set class again
    attr(x, "index_df") <- index
    m <- match(c("indexed_frame", "pdata.frame", "data.frame"), clx, nomatch = 0L)
    oldClass(x) <- c("indexed_frame", if (length(mp <- m[m != 0L])) clx[-mp] else clx, "pdata.frame", if (m[3L]) "data.frame")
    if(any(clx == "data.table")) return(alc(x))
  } else {
    attr(x, "index_df") <- index
    oldClass(x) <- unique.default(c("indexed_series", "pseries", class(x)))
  }
  x
}

# TODO: group for integers use quf..
findex_by <- function(.X, ..., single = "auto", interact.ids = TRUE) { # pid = NULL, t
  clx <- oldClass(.X)
  oldClass(.X) <- NULL
  dots <- substitute(list(...))
  ids <- eval(dots, .X, parent.frame())
  nam <- names(ids)
  vars <- all.vars(dots, unique = FALSE)

  # If something else than NSE cols is supplied
  if(length(ids) == 1L && is.null(nam) && (length(vars) != 1L || !anyv(names(.X), vars))) { # !is.symbol(dots[[2L]]) || length(ids[[1L]]) != length(.X[[1L]]) || is.function(ids[[1L]]) # Fixes #320
    ids <- .X[cols2int(ids[[1L]], .X, names(.X), FALSE)]
  } else {
    if(length(nam)) {
      nonmiss <- nzchar(nam)
      if(!all(nonmiss)) names(ids) <- `[<-`(as.character(dots[-1L]), nonmiss, value = nam[nonmiss])
    } else names(ids) <- vars
  }

  # Single id
  if(length(ids) == 1L) {
    idl <- switch(single, auto = anyDuplicated.default(ids[[1L]]) > 0L, id = TRUE, time = FALSE, stop("'single' must be 'auto', 'id' or 'time'"))
    ids[[1L]] <- if(idl) qF(ids[[1L]], sort = is.factor(ids[[1L]]), na.exclude = FALSE) else make_time_factor(ids[[1L]])
    attr(ids, "single.id") <- idl
  } else {
    lids <- length(ids)
    if(lids > 2L) {
      if(interact.ids) {
        nam <- names(ids)
        ids <- c(`names<-`(list(finteraction(ids[-lids], sort = FALSE)), paste(nam[-lids], collapse = ".")), ids[lids])
        attr(ids, "nam") <- nam # This is a trick, fetched using attr(x, "nam"), before "names" attribute
      } else ids[-lids] <- lapply(ids[-lids], function(x) qF(x, sort = is.factor(x), na.exclude = FALSE))
    } else ids[[1L]] <- qF(ids[[1L]], sort = is.factor(ids[[1L]]), na.exclude = FALSE)
    ids[[length(ids)]] <- make_time_factor(ids[[length(ids)]])
  }
  attr(ids, "row.names") <- .set_row_names(length(ids[[1L]]))
  oldClass(ids) <- c("index_df", "pindex", "data.frame")
  m <- match(c("indexed_frame", "pdata.frame", "data.frame"), clx, nomatch = 0L)
  .X <- index_series(.X, ids, clx)
  attr(.X, "index_df") <- ids
  oldClass(.X) <- c("indexed_frame", if (length(mp <- m[m != 0L])) clx[-mp] else clx, "pdata.frame", if (m[3L]) "data.frame")
  if(any(clx == "data.table")) return(alc(.X))
  .X
}


iby <- findex_by


group_effect <- function(x, effect) {
  index <- findex(x)
  g <- if(length(effect) == 1L) .subset2(index, effect) else .subset(index, effect)
  if(is.factor(g)) return(g)
  g <- group(g)
  attr(g, "levels") <- seq_len(attr(g, "N.groups")) # This is just a trick for fnlevels..
  g
}

uncl2pix <- function(x, interact = FALSE) {
  ix <- unclass(findex(x))
  if(length(ix) == 2L) return(ix)
  if(length(ix) == 1L) {
    res <- if(isTRUE(attr(ix, "single.id"))) list(ix[[1L]], NULL) else list(0L, ix[[1L]])
  } else if(length(ix) > 2L) {
    if(interact) {
      g <- finteraction(ix[-length(ix)])
    } else {
      g <- group(ix[-length(ix)])
      attr(g, "levels") <- seq_len(attr(g, "N.groups"))
    }
    res <- list(g, ix[[length(ix)]])
  } else stop("invalid 'index' length")
  attr(res, "nam") <- names(ix)
  return(res)
}

plm_check_time <- function(x) {
  tlev <- attr(x, "levels")
  oldopts <- options(warn = -1L)
  on.exit(options(oldopts))
  if(is.finite(as.integer(tlev[1L]))) return(as.integer(tlev)[x])
  x
}

pseries_to_numeric <- function(x) {
  clx <- oldClass(x)
  m <- clx %in% c("integer", "logical", "complex", "raw")
  if(any(m)) oldClass(x) <- c(clx[!m], "numeric")
  x
}

unindex <- function(x) {
  attr(x, "index_df") <- NULL
  clx <- oldClass(x)
  if(is.list(x)) {
    oldClass(x) <- fsetdiff(clx, c("indexed_frame", "pdata.frame"))
    x <- fdapply(x, function(y) {
      attr(y, "index_df") <- NULL
      cly <- oldClass(y)
      oldClass(y) <- fsetdiff(cly, c("indexed_series", "pseries", if(length(cly) == 3L) class(unclass(y))))
      y
    })
    if(any(clx == "data.table")) return(alc(x))
  } else {
    oldClass(x) <- fsetdiff(clx, c("indexed_series", "pseries", if(length(clx) == 3L) class(unclass(x))))
  }
  x
}

unindex_light <- function(x) {
  clx <- oldClass(x)
  attr(x, "index_df") <- NULL
  oldClass(x) <- fsetdiff(clx, if(is.list(x)) c("indexed_frame", "pdata.frame") else c("indexed_series", "pseries", if(length(clx) == 3L) class(unclass(x))))
  x
}

index_stats <- function(index) {
  oldClass(index) <- NULL
  lix <- length(index)
  nam <- names(index)
  if(lix > 1L || isFALSE(attr(index, "single.id"))) {
    t <- index[[lix]]
    ndt <- fndistinct(t)
    tstat <- if(ndt == fnlevels(t)) paste0(nam[lix], " [", fnlevels(t), "]") else
      paste0(nam[lix], " [", ndt,  " (", fnlevels(t), ")]")
  } else tstat <- NULL
  if(lix > 1L || isTRUE(attr(index, "single.id"))) {
    if(lix <= 2L) {
      idstat <- paste0(nam[1L], " [", fnlevels(index[[1L]]), "]")
    } else {
      idstat <- paste(paste0(nam[-lix], " [", vapply(index[-lix], fnlevels, integer(1L)), "]"), collapse = " ")
    }
  } else idstat <- NULL
  return(paste(c(idstat, tstat), collapse = " | "))
}

print.indexed_series <- function(x, ...) {
  print(unindex_light(x), ...)
  # if(inherits(index, "pindex")) {
  cat("\nIndexed by: ", index_stats(findex(x)), "\n")
  # }
}

print.indexed_frame <- function(x, ...) {
  print(unindex(x), ...)
  # if(inherits(index, "pindex")) {
  cat("\nIndexed by: ", index_stats(findex(x)), "\n")
  # }
}


droplevels_index <- function(index, drop.index.levels = "id") {
  oi <- switch(drop.index.levels, none = 0L, id = 1L, time = 2L, all = 3L, stop("drop.index.levels must be one of 'all', 'id', 'time' or 'none'.") )
  if(oi == 0L) return(index)
  clix <- oldClass(index)
  oldClass(index) <- NULL
  if(oi == 1L) {
    if(length(index) > 2L) index[-length(index)] <- fdroplevels.list(index[-length(index)])
    else if(length(index) == 2L || isTRUE(attr(index, "single.id"))) index[[1L]] <- fdroplevels(index[[1L]])
  } else if(oi == 2L) {
    index[[length(index)]] <- fdroplevels(index[[length(index)]])
  } else {
    index <- fdroplevels.list(index)
  }
  oldClass(index) <- clix
  index
}

`[.index_df` <- function(x, i, j, drop = FALSE, drop.index.levels = "id") {
  res <- droplevels_index(ss(x, i, j), drop.index.levels)
  lr <- length(unclass(res))
  if(drop && lr == 1L) return(.subset(res, 1L))
  if(lr == 1L && length(unclass(x)) > 1L) {
      attr(res, "single.id") <- attr(res, "names") != l1orlst(attr(x, "names"))
  }
  res
}

print.index_df <- function(x, topn = 5, ...) {
   oldClass(x) <- "data.frame"
   if(fnrow(x) > 2*topn) {
     print(head(x, topn), ...)
     cat("---")
     print(`names<-`(tail(x, topn), NULL), ...)
   } else print(x, ...)
   cat("\n", index_stats(x), "\n", sep = "")
}

`[.indexed_frame` <- function(x, i, ..., drop.index.levels = "id") {

  clx <- oldClass(x)
  idDTl <- any(clx == "data.table")

  if(idDTl) {
    res <- unindex_light(x)
    # res <- NextMethod() # doesn't work with i
    ivsbl <- any(clx == "invisible")
    if(ivsbl) clx <- clx[clx != "invisible"] # for chaining...
    if(!missing(...)) {
      rem <- as.list(substitute(list(...))[-1L])
      cal <- as.call(c(list(quote(`[`), quote(res), substitute(i)), rem))
      rem <- as.character(rem)
      if(any(grepl(".SD", rem)) && !any(grepl("apply", rem)))
        warning("Found '.SD' in the call but no 'apply' function. Please note that .SD is not an indexed_frame but a plain data.table containing indexed_series. Thus indexed_frame / pdata.frame methods don't work on .SD! Consider using (m/l)apply(.SD, FUN) or reindex(.SD, ix(data)). If you are not performing indexed operations on .SD please ignore or suppress this warning.")
      if(any(grepl(":=", rem))) {
        res <- copyMostAttributes(eval(cal, list(res = alc(res)), parent.frame()), x)
        eval.parent(substitute(x <- res))
        oldClass(res) <- c("invisible", clx)
        return(res)
      }
    } else cal <- as.call(list(quote(`[`), quote(res), substitute(i)))
    res <- eval(cal, list(res = res), parent.frame())
    if(missing(i) && fnrow(res) != fnrow(x)) {
      if(ivsbl) oldClass(res) <- fsetdiff(oldClass(res), "invisible")
      return(unindex(res)) # data.table aggregation
    } else if(!missing(i)) i <- eval(substitute(i), x, parent.frame())
  } else res <- unindex(x)[i, ...] # does not respect data.table properties, but better for sf data frame and others which might check validity of "index_df" attribute

  index <- attr(x, "index_df")

  if(!missing(i) && (is.atomic(res) || fnrow(res) != fnrow(x) || length(i) == fnrow(x))) { # Problem: mtcars[1:10] selects columns, not rows!!
    index <- droplevels_index(ss(index, i), drop.index.levels)
    if(is.list(res)) {
      if(fnrow(res) != fnrow(index)) return(unindex(res)) # could be with data.table using i and also aggregating in j
      res <- index_series(res, index, clx)
    }
  } else if(!idDTl && is.list(res)) res <- index_series(res, index, clx)

  attr(res, "index_df") <- index

  if(is.atomic(res)) {
    oldClass(res) <- unique.default(c("indexed_series", "pseries", class(res)))
    return(res)
  }

  m <- match(c("indexed_frame", "pdata.frame", "data.frame"), clx, nomatch = 0L)
  oldClass(res) <- c("indexed_frame", if (length(mp <- m[m != 0L])) clx[-mp] else clx, "pdata.frame", if (m[3L]) "data.frame")

  if(any(clx == "data.table")) return(alc(res))
  res
}

`[.indexed_series` <- function(x, i, ..., drop.index.levels = "id") {

  res <- unindex_light(x)[i, ...]  # NextMethod("[", x, ...) # plm has no [.pseries method yet, but the drop.index.levels argument causes problems...
  if(length(res) <= 1L) return(res)

  if(!missing(i)) {
    attr(res, "index_df") <- droplevels_index(ss(findex(x), i), drop.index.levels)
  } else if(is.null(attr(res, "index_df"))) {
    attr(res, "index_df") <- findex(x)
  }

  oldClass(res) <- c("indexed_series", "pseries", class(res))
  res
}

`$.indexed_frame` <- function(x, name) {
  # res <- NextMethod() # don't use pdata.frame methods
  res <- .subset2(x, name, exact = FALSE) # as.character(substitute(name)) -> not necessary!
  if(is.null(res)) return(NULL)
  clr <- class(res)
  attr(res, "index_df") <- attr(x, "index_df")
  if(!any(clr == "indexed_series")) oldClass(res) <- c("indexed_series", "pseries", clr)
  res
}

`$<-.indexed_frame` <- function(x, name, value) {
  clx <- oldClass(x)
  oldClass(x) <- NULL
  if(is.null(value)) {
    x[[name]] <- NULL
    oldClass(x) <- clx
    if(any(clx == "data.table")) return(alc(x)) else return(x)
  }
  if(length(value) != .Call(C_fnrow, x)) {
    if(length(value) == 1L) value <- alloc(value, .Call(C_fnrow, x))
    else stop("length(value) must match nrow(x)")
  }
  attr(value, "index_df") <- .Call(C_createeptr, attr(x, "index_df"))
  oldClass(value) <- unique.default(c("indexed_series", "pseries", class(value)))
  x[[name]] <- value
  oldClass(x) <- clx
  if(any(clx == "data.table")) return(alc(x)) else return(x)
}

# What about i and j for data.frame?
`[[.indexed_frame` <- function(x, i, ...) {
  # res <- NextMethod() # don't use pdata.frame methods
  # oldClass(x) <- fsetdiff(oldClass(x), c("indexed_frame", "pdata.frame"))
  # res <- UseMethod("[[", x)
  res <- .subset2(x, i, ...)
  if(is.null(res)) return(NULL)
  clr <- class(res)
  attr(res, "index_df") <- attr(x, "index_df")
  if(!any(clr == "indexed_series")) oldClass(res) <- c("indexed_series", "pseries", clr)
  res
}

# No plm method... can use NextMethod? -> Yes, but this is faster, and I don't know of any other use cases (nobody uses df[[i, j]])
`[[<-.indexed_frame` <- function(x, i, value) {
  clx <- oldClass(x)
  oldClass(x) <- NULL
  if(is.null(value)) {
    x[[i]] <- NULL
    oldClass(x) <- clx
    if(any(clx == "data.table")) return(alc(x)) else return(x)
  }
  if(length(value) != .Call(C_fnrow, x)) {
    if(length(value) == 1L) value <- alloc(value, .Call(C_fnrow, x))
    else stop("length(value) must match nrow(x)")
  }
  attr(value, "index_df") <- .Call(C_createeptr, attr(x, "index_df"))
  oldClass(value) <- unique.default(c("indexed_series", "pseries", class(value)))
  x[[i]] <- value
  oldClass(x) <- clx
  if(any(clx == "data.table")) return(alc(x)) else return(x)
}

# no plm method... can use NextMethod!
`[<-.indexed_frame` <- function(x, i, j, value) {
  res <- NextMethod()
  if(missing(j)) return(res)
  if(!(missing(i) || missing(j)) && identical(attr(x, "names"), attr(res, "names"))) return(res)
  return(reindex(res))
}

# These are primarily needed to overwrite pseries methods when plm is attached...
# Note: could use reindex() instead of duplAttributes(), but the latter is more efficient,
# and I can't think of a single example where it would be undesirable.

Math.indexed_series <- function(x, ...) {
  duplAttributes(get(.Generic)(unindex_light(x), ...), x)
}

Ops.indexed_series <- function(e1, e2) {
  if(missing(e2)) { # unary operators (+, - and !)
    res <- get(.Generic)(unindex_light(e1))
    if(.Generic == "!") return(res)
    return(duplAttributes(res, e1))
  }
  res <- get(.Generic)(unindex_light(e1), unindex_light(e2))
  if(!any(.Generic == c("+", "-", "*", "/", "^", "%%", "%/%"))) return(res)
  if(inherits(e1, "indexed_series")) return(duplAttributes(res, e1))
  duplAttributes(res, e2)
}
